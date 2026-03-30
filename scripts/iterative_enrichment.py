"""
迭代富集管线 — 四大任务合一 (Iterative Enrichment Pipeline)
=============================================================

Task 1: LOTUS 3D SMILES 升维 (Stereo-SMILES Upgrade)
Task 2: 糖链二次解析 (Second-Pass Sugar Resolution)
Task 3: ChEMBL 活性二次挖掘 (Second-Pass ChEMBL Mining)
Task 4: PubChem API 兜底 (PubChem Fallback)

输入: reports/GlycoNP_Imputed_Merged.csv  (经 LOTUS 回填后的数据)
输出: reports/GlycoNP_Fully_Enriched.csv  (终极富集版本)

使用方法 / Usage:
  python scripts/iterative_enrichment.py
  python scripts/iterative_enrichment.py --skip-pubchem         # 跳过 Task 4
  python scripts/iterative_enrichment.py --pubchem-limit 500    # PubChem 查询限制
"""
import argparse
import os
import re
import sqlite3
import sys
import time
from collections import Counter
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
REPORT_DIR = os.path.join(BASE_DIR, "reports")
LOTUS_GZ = os.path.join(BASE_DIR, "data", "230106_frozen_metadata.csv.gz")
CHEMBL_DB = os.path.join(BASE_DIR, "data", "chembl_36",
                         "chembl_36_sqlite", "chembl_36.db")


# =====================================================================
# Task 1: LOTUS 3D SMILES 升维 (Stereo-SMILES Upgrade)
# =====================================================================

def runTask1_stereoUpgrade(df: pd.DataFrame) -> pd.DataFrame:
    """
    从 LOTUS 提取带立体化学的 SMILES, 替换我们的 2D 残次 SMILES。
    Extract stereo-rich SMILES from LOTUS to replace our flat 2D SMILES.

    替换条件:
      - 我们的 SMILES 不含 @ 符号 (无手性标记)
      - LOTUS 的 SMILES 包含 @ 符号 (有手性标记)
      - InChIKey Block-1 匹配 (同一分子骨架)

    安全: 原始 standard_inchi_key 不被覆盖。3D SMILES 写入新列,
    同时更新 canonical_smiles 列。原始 SMILES 保存在 original_smiles_2d。
    """
    print("\n" + "=" * 70)
    print("  Task 1: LOTUS 3D SMILES 升维 / Stereo-SMILES Upgrade")
    print("=" * 70)
    t0 = time.time()

    if not os.path.exists(LOTUS_GZ):
        print(f"  [SKIP] LOTUS file not found: {LOTUS_GZ}")
        return df

    # 加载 LOTUS 的 SMILES 列 (Load LOTUS SMILES columns)
    print("  [1/4] Loading LOTUS SMILES...")
    lotusDf = pd.read_csv(
        LOTUS_GZ, dtype=str, low_memory=False,
        usecols=["structure_inchikey", "structure_smiles", "structure_smiles_2D"])
    print(f"    LOTUS rows: {len(lotusDf):,}")

    # 生成 LOTUS Block-1 (Generate LOTUS Block-1)
    lotusDf["_block1"] = lotusDf["structure_inchikey"].astype(str).str[:14]
    lotusDf = lotusDf.dropna(subset=["_block1"])
    lotusDf = lotusDf[lotusDf["_block1"] != "nan"]

    # 优先使用 structure_smiles (3D), 回退到 structure_smiles_2D
    # Prefer structure_smiles (3D), fallback to structure_smiles_2D
    def pickBestSmiles(row):
        smi3d = str(row.get("structure_smiles", ""))
        smi2d = str(row.get("structure_smiles_2D", ""))
        if smi3d and smi3d != "nan" and "@" in smi3d:
            return smi3d
        if smi2d and smi2d != "nan" and "@" in smi2d:
            return smi2d
        return ""

    lotusDf["_best_smiles"] = lotusDf.apply(pickBestSmiles, axis=1)

    # 只保留有手性标记的 (Only keep those with stereo markers)
    lotusDf = lotusDf[lotusDf["_best_smiles"].str.contains("@", na=False)]
    lotusDf = lotusDf.drop_duplicates(subset=["_block1"])
    print(f"    LOTUS entries with stereo: {len(lotusDf):,}")

    # 构建查找表 (Build lookup table)
    print("  [2/4] Building Block-1 → 3D SMILES lookup...")
    lotusMap = dict(zip(lotusDf["_block1"], lotusDf["_best_smiles"]))
    print(f"    Lookup table size: {len(lotusMap):,}")

    # 识别需要升级的行 (Identify rows needing upgrade)
    print("  [3/4] Identifying upgrade candidates...")
    needs2dRescue = ~df["canonical_smiles"].astype(str).str.contains("@", na=False)
    hasBlock1 = df["non_isomeric_inchikey_block1"].notna()
    candidates = df[needs2dRescue & hasBlock1].index

    print(f"    SMILES without stereo: {needs2dRescue.sum():,}")
    print(f"    With Block-1 key: {hasBlock1.sum():,}")
    print(f"    Candidates: {len(candidates):,}")

    # 执行升维 (Execute upgrade)
    print("  [4/4] Upgrading to 3D SMILES...")
    df["original_smiles_2d"] = ""  # 备份列 (backup column)
    df["Is_3D_Rescued"] = False

    upgraded = 0
    for idx in candidates:
        block1 = str(df.at[idx, "non_isomeric_inchikey_block1"])
        lotusSmiles = lotusMap.get(block1)
        if lotusSmiles and "@" in lotusSmiles:
            # 保存原始 SMILES (备份)
            df.at[idx, "original_smiles_2d"] = df.at[idx, "canonical_smiles"]
            # 替换为 3D SMILES
            df.at[idx, "canonical_smiles"] = lotusSmiles
            df.at[idx, "Is_3D_Rescued"] = True
            upgraded += 1

    elapsed = time.time() - t0
    print(f"\n  Task 1 Results:")
    print(f"    3D Upgraded: {upgraded:,} / {len(candidates):,}")
    print(f"    Upgrade rate: {upgraded / len(candidates) * 100:.1f}%"
          if len(candidates) > 0 else "    N/A")
    print(f"    Time: {elapsed:.0f}s")

    return df


# =====================================================================
# Task 2: 糖链二次解析 (Second-Pass Sugar Resolution)
# =====================================================================

# 从 rescue_generic_sugars.py 复用的名称挖掘字典
# Reused from rescue_generic_sugars.py
TEXT_MINING_RULES = [
    (r'glucuronic\s*acid|glucuronide|glucuronosyl', 'D-GlcA'),
    (r'galacturonic\s*acid|galacturonide', 'D-GalA'),
    (r'N-acetylglucosamin|GlcNAc', 'D-GlcNAc'),
    (r'N-acetylgalactosamin|GalNAc', 'D-GalNAc'),
    (r'neuraminic|sialic|Neu5Ac|NeuAc', 'Neu5Ac'),
    (r'glucopyranosid|glucosid|glucosyl|\bgluco(?:se)?\b', 'D-Glc'),
    (r'galactopyranosid|galactosid|galactosyl|\bgalacto(?:se)?\b', 'D-Gal'),
    (r'mannopyranosid|mannosid|mannosyl|\bmanno(?:se)?\b', 'D-Man'),
    (r'xylopyranosid|xylosid|xylosyl|\bxylo(?:se)?\b', 'D-Xyl'),
    (r'arabinopyranosid|arabinosid|arabinosyl|\barabino(?:se)?\b', 'L-Ara'),
    (r'rhamnopyranosid|rhamnosid|rhamnosyl|\brhamno(?:se)?\b', 'L-Rha'),
    (r'fucopyranosid|fucosid|fucosyl|\bfuco(?:se)?\b', 'L-Fuc'),
    (r'ribopyranosid|ribosid|ribosyl|\bribo(?:se)?\b', 'D-Rib'),
    (r'fructopyranosid|fructosid|fructosyl|\bfructo(?:se)?\b', 'D-Fru'),
    (r'digitalosid|digitalose', 'D-Dig'),
    (r'oleandrosid|oleandrose', 'L-Ole'),
    (r'cymarosid|cymarose', 'D-Cym'),
    (r'thevetosid|thevetose', 'D-The'),
    (r'quinovo(?:se|sid)', 'D-Qui'),
]
COMPILED_TEXT_RULES = [
    (re.compile(p, re.IGNORECASE), s) for p, s in TEXT_MINING_RULES
]
GENERIC_PATTERN = re.compile(r'\b(Hex|Pen|dHex|HexA|Non|Oct|Hept)\b')


def _nameRescueSingleToken(
    token: str,
    namePool: str,
) -> str:
    """
    用化合物名称文本挖掘替换一个泛指 token。
    Text-mine a specific sugar name to replace a generic token.
    """
    if token not in ("Hex", "Pen", "dHex", "HexA", "Non", "Oct", "Hept"):
        return token
    for regex, sugar in COMPILED_TEXT_RULES:
        if regex.search(namePool):
            return sugar
    return token


def runTask2_secondPassSugar(df: pd.DataFrame) -> pd.DataFrame:
    """
    二次解析: 利用新 3D SMILES + 新 name 重新跑名称糖链匹配。
    Second pass: use upgraded 3D SMILES + backfilled names for sugar re-resolution.
    """
    print("\n" + "=" * 70)
    print("  Task 2: 糖链二次解析 / Second-Pass Sugar Resolution")
    print("=" * 70)
    t0 = time.time()

    # 统计当前泛指标签残留 (Count current generic labels)
    seqCol = "Sugar_Sequence"
    if "Rescued_Sugar_Sequence" in df.columns:
        seqCol = "Rescued_Sugar_Sequence"

    allSeqs = df[seqCol].dropna().astype(str)
    beforeTokens = []
    for seq in allSeqs:
        beforeTokens.extend(GENERIC_PATTERN.findall(seq))
    beforeCounter = Counter(beforeTokens)
    totalGenericBefore = sum(beforeCounter.values())

    print(f"  Before second pass:")
    print(f"    Total generic tokens: {totalGenericBefore:,}")
    for k, v in beforeCounter.most_common():
        print(f"      {k:>6s}: {v:,}")

    # 两段式名称修复 (Two-stage name rescue)
    # Stage A: 利用新回填的 name/synonyms 修复泛指标签
    # Stage B: 特别针对 Is_3D_Rescued=True 的记录重新尝试

    print(f"\n  Running name-based sugar rescue...")
    rescued = 0
    nameColCandidates = ["name", "iupac_name", "synonyms",
                         "structure_nameTraditional"]

    for idx in df.index:
        seq = str(df.at[idx, seqCol]) if pd.notna(df.at[idx, seqCol]) else ""
        if not GENERIC_PATTERN.search(seq):
            continue

        # 构建文本池 (Build text pool)
        namePool = ""
        for col in nameColCandidates:
            if col in df.columns and pd.notna(df.at[idx, col]):
                val = str(df.at[idx, col])
                if val not in ("", "nan", "None"):
                    namePool += " " + val

        if not namePool.strip():
            continue

        # 替换泛指 token (Replace generic tokens)
        def replaceMatch(m):
            return _nameRescueSingleToken(m.group(0), namePool)

        newSeq = GENERIC_PATTERN.sub(replaceMatch, seq)
        if newSeq != seq:
            df.at[idx, seqCol] = newSeq
            rescued += 1

    # 统计修复后泛指标签残留 (Count remaining generic labels)
    allSeqs = df[seqCol].dropna().astype(str)
    afterTokens = []
    for seq in allSeqs:
        afterTokens.extend(GENERIC_PATTERN.findall(seq))
    afterCounter = Counter(afterTokens)
    totalGenericAfter = sum(afterCounter.values())
    resolved = totalGenericBefore - totalGenericAfter

    elapsed = time.time() - t0
    print(f"\n  Task 2 Results:")
    print(f"    Rows modified: {rescued:,}")
    print(f"    Generic tokens resolved: {resolved:,}")
    print(f"    Remaining generic tokens: {totalGenericAfter:,}")
    print(f"    Resolution rate: "
          f"{resolved / totalGenericBefore * 100:.1f}%"
          if totalGenericBefore > 0 else "    N/A")
    print(f"\n    Remaining breakdown:")
    for k, v in afterCounter.most_common():
        print(f"      {k:>6s}: {v:,}")
    print(f"    Time: {elapsed:.0f}s")

    return df


# =====================================================================
# Task 3: ChEMBL 活性二次挖掘 (Second-Pass ChEMBL Mining)
# =====================================================================

def runTask3_secondPassChembl(df: pd.DataFrame) -> pd.DataFrame:
    """
    利用补全后的 InChIKey, 重新查询 ChEMBL 活性数据。
    Re-query ChEMBL activity data using enriched InChIKeys.
    """
    print("\n" + "=" * 70)
    print("  Task 3: ChEMBL 活性二次挖掘 / Second-Pass ChEMBL Mining")
    print("=" * 70)
    t0 = time.time()

    if not os.path.exists(CHEMBL_DB):
        print(f"  [SKIP] ChEMBL DB not found: {CHEMBL_DB}")
        return df

    # 统计首轮已有活性数据的化合物数 (Count first-pass hits)
    hasActivity = (
        df.get("ChEMBL_N_Activities", pd.Series(dtype=str)).notna() &
        (~df.get("ChEMBL_N_Activities", pd.Series(dtype=str))
         .astype(str).isin(["", "nan", "None", "0"]))
    )
    firstPassHits = hasActivity.sum()
    print(f"  First-pass compounds with activity: {firstPassHits:,}")

    # 找出首轮没有活性数据但有有效 InChIKey 的行
    # Find rows without activity but with valid InChIKey
    noActivity = ~hasActivity
    validIk = (
        df["standard_inchi_key"].notna() &
        (~df["standard_inchi_key"].astype(str).isin(["", "nan", "None"]))
    )
    candidates = df[noActivity & validIk]
    print(f"  Candidates for second pass: {len(candidates):,}")

    if len(candidates) == 0:
        print(f"  No new candidates!")
        return df

    # 收集 InChIKey (Collect InChIKeys)
    inchiKeys = candidates["standard_inchi_key"].astype(str).unique().tolist()
    print(f"  Unique InChIKeys to query: {len(inchiKeys):,}")

    # 批量查询 ChEMBL (Batch query ChEMBL)
    conn = sqlite3.connect(CHEMBL_DB, timeout=30)
    conn.execute("PRAGMA query_only = ON")

    # Step 1: InChIKey 碰撞 (InChIKey collision)
    print(f"  [1/2] InChIKey collision...")
    batchSize = 500
    ikToMolregno = {}
    for i in range(0, len(inchiKeys), batchSize):
        batch = inchiKeys[i:i + batchSize]
        placeholders = ",".join(["?"] * len(batch))
        sql = (f"SELECT standard_inchi_key, molregno "
               f"FROM compound_structures "
               f"WHERE standard_inchi_key IN ({placeholders})")
        cursor = conn.execute(sql, batch)
        for row in cursor:
            ikToMolregno[row[0]] = row[1]

    print(f"    ChEMBL hits: {len(ikToMolregno):,}")

    if not ikToMolregno:
        conn.close()
        print(f"  No new ChEMBL matches")
        return df

    # Step 2: 批量提取活性 (Batch extract activities)
    print(f"  [2/2] Extracting activities...")
    molregnos = list(ikToMolregno.values())
    allActivities = []

    for i in range(0, len(molregnos), batchSize):
        batch = molregnos[i:i + batchSize]
        placeholders = ",".join(["?"] * len(batch))
        sql = f"""
            SELECT
                md.molregno,
                td.pref_name AS target_name,
                act.standard_type,
                act.pchembl_value
            FROM activities act
            JOIN assays a ON act.assay_id = a.assay_id
            JOIN target_dictionary td ON a.tid = td.tid
            JOIN molecule_dictionary md ON act.molregno = md.molregno
            WHERE act.molregno IN ({placeholders})
              AND act.standard_type IN ('IC50','Ki','EC50','Kd')
              AND act.pchembl_value IS NOT NULL
        """
        cursor = conn.execute(sql, batch)
        allActivities.extend(cursor.fetchall())

    conn.close()
    print(f"    New activity rows: {len(allActivities):,}")

    if not allActivities:
        print(f"  No new activity data found")
        return df

    # 聚合活性数据 (Aggregate activity data)
    actDf = pd.DataFrame(allActivities,
                          columns=["molregno", "target_name",
                                   "standard_type", "pchembl_value"])
    actDf["pchembl_value"] = pd.to_numeric(actDf["pchembl_value"],
                                            errors="coerce")

    # 反向映射: molregno → inchikey
    molregnoToIk = {v: k for k, v in ikToMolregno.items()}

    newFilled = 0
    for ik, mregno in ikToMolregno.items():
        subDf = actDf[actDf["molregno"] == mregno]
        if subDf.empty:
            continue

        # 构建 summary (Build summary)
        topActs = subDf.nlargest(3, "pchembl_value")
        summaryParts = []
        for _, row in topActs.iterrows():
            summaryParts.append(
                f"{row['target_name']}({row['standard_type']}="
                f"{row['pchembl_value']:.1f})")
        summary = "; ".join(summaryParts)

        targets = "; ".join(subDf["target_name"].dropna().unique()[:5])
        actTypes = "; ".join(subDf["standard_type"].dropna().unique())
        bestPchembl = subDf["pchembl_value"].max()
        nActs = len(subDf)

        # 写入主数据集 (Write to main dataset)
        mask = df["standard_inchi_key"] == ik
        if mask.sum() > 0:
            df.loc[mask, "bioactivity_summary"] = summary
            df.loc[mask, "ChEMBL_Targets"] = targets
            df.loc[mask, "ChEMBL_Activity_Types"] = actTypes
            df.loc[mask, "ChEMBL_Best_pChEMBL"] = str(bestPchembl)
            df.loc[mask, "ChEMBL_N_Activities"] = str(nActs)
            newFilled += mask.sum()

    # 统计 (Statistics)
    hasActivityAfter = (
        df.get("ChEMBL_N_Activities", pd.Series(dtype=str)).notna() &
        (~df.get("ChEMBL_N_Activities", pd.Series(dtype=str))
         .astype(str).isin(["", "nan", "None", "0"]))
    )
    secondPassHits = hasActivityAfter.sum()
    newHits = secondPassHits - firstPassHits

    # pChEMBL>7 统计 (Count high-activity)
    pchemblVals = pd.to_numeric(
        df["ChEMBL_Best_pChEMBL"], errors="coerce")
    highActivity = (pchemblVals > 7).sum()

    elapsed = time.time() - t0
    print(f"\n  Task 3 Results:")
    print(f"    First-pass compounds with activity: {firstPassHits:,}")
    print(f"    After second-pass: {secondPassHits:,}")
    print(f"    New compounds with activity: +{newHits:,}")
    print(f"    Total high-activity (pChEMBL>7): {highActivity:,}")
    print(f"    Time: {elapsed:.0f}s")

    return df


# =====================================================================
# Task 4: PubChem API 兜底 (PubChem Fallback)
# =====================================================================

def runTask4_pubchemFallback(
    df: pd.DataFrame,
    limit: int = 300,
) -> pd.DataFrame:
    """
    PubChem API 兜底: 查询残留空白的 organisms/dois。
    PubChem API fallback for remaining empty organisms/dois.

    安全限制:
      - time.sleep(0.3) 每次请求
      - HTTP 429 退避重试
      - 严格限制查询数量
    """
    print("\n" + "=" * 70)
    print("  Task 4: PubChem API 兜底 / PubChem Fallback")
    print("=" * 70)
    t0 = time.time()

    try:
        import pubchempy as pcp
    except ImportError:
        print(f"  [ERROR] pubchempy not installed: pip install pubchempy")
        return df

    # 识别候选行: 有 name 但 organisms 仍空
    hasName = (
        df["name"].notna() &
        (~df["name"].astype(str).isin(["", "nan", "None"]))
    )
    needsOrg = (
        df["organisms"].isna() |
        df["organisms"].astype(str).isin(["", "nan", "None"])
    )
    candidates = df[hasName & needsOrg].head(limit)
    print(f"  Candidates (has name, needs organisms): {len(candidates):,}")
    print(f"  Query limit: {limit}")

    if len(candidates) == 0:
        print(f"  No candidates!")
        return df

    filled = 0
    errors = 0
    retries = 0

    for i, idx in enumerate(candidates.index):
        nameStr = str(df.at[idx, "name"]).strip()
        if not nameStr or nameStr == "nan":
            continue

        # 安全限流 (Rate limiting)
        time.sleep(0.3)

        try:
            results = pcp.get_compounds(nameStr, "name")
            if results:
                compound = results[0]
                cid = str(compound.cid) if compound.cid else ""
                if cid:
                    df.at[idx, "PubChem_CID"] = cid
                    filled += 1
            else:
                errors += 1
        except Exception as e:
            errStr = str(e)
            if "429" in errStr or "Too Many" in errStr:
                # HTTP 429 退避重试 (Backoff retry)
                retries += 1
                print(f"    [429] Rate limited at query {i+1}, "
                      f"backing off 10s...")
                time.sleep(10)
                try:
                    results = pcp.get_compounds(nameStr, "name")
                    if results and results[0].cid:
                        df.at[idx, "PubChem_CID"] = str(results[0].cid)
                        filled += 1
                    else:
                        errors += 1
                except Exception:
                    errors += 1
            else:
                errors += 1

        if (i + 1) % 50 == 0:
            print(f"    Progress: {i+1}/{len(candidates)} "
                  f"(filled={filled}, errors={errors}, retries={retries})")

    elapsed = time.time() - t0
    print(f"\n  Task 4 Results:")
    print(f"    Queried: {len(candidates):,}")
    print(f"    CID filled: {filled:,}")
    print(f"    Errors: {errors:,}")
    print(f"    429 retries: {retries:,}")
    print(f"    Time: {elapsed:.0f}s")

    return df


# =====================================================================
# 终极审查报告 (Final Audit Report)
# =====================================================================

def printFinalAudit(df: pd.DataFrame, total: int):
    """输出全库核心字段完整度。 / Print final audit of all key fields."""
    print("\n" + "=" * 70)
    print("  终极审查报告 / Final Audit Report")
    print("=" * 70)

    def pct(n):
        return f"{n / total * 100:.1f}%"

    def countValid(col):
        if col not in df.columns:
            return 0
        return (df[col].notna() &
                ~df[col].astype(str).isin(["", "nan", "None"])).sum()

    fields = [
        ("canonical_smiles", "SMILES"),
        ("standard_inchi_key", "InChIKey"),
        ("non_isomeric_inchikey_block1", "2D Block-1"),
        ("name", "Name"),
        ("organisms", "Organisms"),
        ("dois", "DOI"),
        ("Sugar_Sequence", "Sugar_Sequence"),
        ("Rescued_Sugar_Sequence", "Rescued_Sugar_Seq"),
        ("bioactivity_summary", "Bioactivity"),
        ("ChEMBL_N_Activities", "ChEMBL_Activities"),
        ("Is_3D_Rescued", "3D_Rescued"),
        ("LOTUS_kingdom", "LOTUS_kingdom"),
        ("LOTUS_family", "LOTUS_family"),
        ("PubChem_CID", "PubChem_CID"),
    ]

    print(f"\n  {'Field':<22} {'Valid':>10} {'%':>8}")
    print(f"  {'-'*22} {'-'*10} {'-'*8}")

    for col, label in fields:
        if col not in df.columns:
            continue
        if col == "Is_3D_Rescued":
            n = (df[col].astype(str) == "True").sum()
        else:
            n = countValid(col)
        print(f"  {label:<22} {n:>10,} {pct(n):>8}")

    # 糖链泛指标签残留统计 (Generic sugar label residual stats)
    seqCol = "Rescued_Sugar_Sequence" if "Rescued_Sugar_Sequence" in df.columns else "Sugar_Sequence"
    if seqCol in df.columns:
        allTokens = []
        for seq in df[seqCol].dropna().astype(str):
            allTokens.extend(GENERIC_PATTERN.findall(seq))
        tokenCounter = Counter(allTokens)
        totalGeneric = sum(tokenCounter.values())
        print(f"\n  Generic sugar labels remaining: {totalGeneric:,}")
        for k, v in tokenCounter.most_common():
            print(f"    {k:>6s}: {v:,}")

    print(f"\n  Total rows: {total:,}")
    print("=" * 70)


# =====================================================================
# Main
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Iterative Enrichment Pipeline")
    parser.add_argument("--input", type=str, default=None)
    parser.add_argument("--output", type=str, default=None)
    parser.add_argument("--skip-pubchem", action="store_true")
    parser.add_argument("--pubchem-limit", type=int, default=300)
    args = parser.parse_args()

    inputPath = args.input or os.path.join(
        REPORT_DIR, "GlycoNP_Imputed_Merged.csv")
    outputPath = args.output or os.path.join(
        REPORT_DIR, "GlycoNP_Fully_Enriched.csv")

    print("=" * 70)
    print("  GlycoNP Iterative Enrichment Pipeline")
    print("  糖缀合物迭代富集管线")
    print("=" * 70)
    print(f"  Input:  {inputPath}")
    print(f"  Output: {outputPath}")

    if not os.path.exists(inputPath):
        print(f"\n  [ERROR] Input not found: {inputPath}")
        return

    t0 = time.time()
    df = pd.read_csv(inputPath, dtype=str, low_memory=False,
                     encoding="utf-8-sig")
    total = len(df)
    print(f"  Loaded: {total:,} rows")

    # ============ Task 1: LOTUS 3D SMILES ============
    df = runTask1_stereoUpgrade(df)

    # ============ Task 2: 糖链二次解析 ============
    df = runTask2_secondPassSugar(df)

    # ============ Task 3: ChEMBL 二次挖掘 ============
    df = runTask3_secondPassChembl(df)

    # ============ Task 4: PubChem 兜底 ============
    if not args.skip_pubchem:
        df = runTask4_pubchemFallback(df, limit=args.pubchem_limit)
    else:
        print(f"\n  [SKIP] Task 4: PubChem fallback disabled (--skip-pubchem)")

    # ============ 保存 (Save) ============
    df.to_csv(outputPath, index=False, encoding="utf-8-sig")
    elapsed = time.time() - t0

    # ============ 终极审查 ============
    printFinalAudit(df, total)

    print(f"\n  Output: {outputPath}")
    print(f"  Total pipeline time: {elapsed:.0f}s")


if __name__ == "__main__":
    main()
