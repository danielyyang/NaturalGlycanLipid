"""
本地数据降维补全引擎 — 无手性 InChIKey + 离线知识库回填 + PubChem 兜底
Dimensionality Reduction Imputation Engine

在进行 ChEMBL 活性挖掘之前, 先解决数据集中的三大缺口:
  - DOI 缺失 58.1% (54,720 条)
  - Organisms 缺失 48.4% (45,622 条)
  - Name 缺失 62.2% (58,573 条)

三大任务 / Three Tasks:
  Task 1: 生成无手性二维 InChIKey Block-1 (万能二维钥匙)
  Task 2: 对接离线知识库 (LOTUS/COCONUT/NPAtlas) 进行 DOI/Name/Organisms 回填
  Task 3: PubChem CID 名称反向兜底 (pubchempy, 可选离线映射表)

使用方法 / Usage:
  # 完整运行 (全部三个任务)
  python scripts/impute_2d_key_and_backfill.py

  # 仅生成 2D key (跳过离线合并和 PubChem)
  python scripts/impute_2d_key_and_backfill.py --task1-only

  # 指定离线数据库路径
  python scripts/impute_2d_key_and_backfill.py --lotus data/external/lotus_frozen.csv.gz

  # 启用 PubChem 兜底 (需要网络, --pubchem-limit 控制数量)
  python scripts/impute_2d_key_and_backfill.py --enable-pubchem --pubchem-limit 200
"""
import argparse
import os
import re
import sys
import time
from collections import Counter
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd
import numpy as np

from rdkit import Chem
from rdkit.Chem.inchi import MolToInchi, InchiToInchiKey

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


# =====================================================================
# Task 1: 生成万能二维钥匙 (Robust 2D Identifier Generation)
# =====================================================================

def generateNonIsomericInchikeyBlock1(smiles: str) -> Optional[str]:
    """
    从 SMILES 生成无手性的二维 InChIKey Block-1 (前 14 字符)。
    Generate non-isomeric 2D InChIKey Block-1 from SMILES.

    流程 / Pipeline:
      1. RDKit 解析 SMILES
      2. RemoveStereochemistry() — 剥离所有立体化学 (E/Z, R/S)
      3. MolToInchi() → InchiToInchiKey() → 取前 14 字符
      (前 14 字符 = connectivity layer, 与 3D 构型无关)

    设计意图: 不同手性的同一骨架在二维 InChIKey Block-1 上完全一致。
    这允许跨数据库关联那些因手性标注差异而匹配失败的分子。
    """
    if not smiles or str(smiles) in ("nan", "", "None"):
        return None

    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return None

        # 关键: 强制剥离所有立体化学信息和同位素标记
        # Critical: force-strip all stereochemistry and isotope labels
        Chem.RemoveStereochemistry(mol)
        # 同时清除同位素 (Also clear isotope labels)
        for atom in mol.GetAtoms():
            atom.SetIsotope(0)

        # 生成非异构 InChI → InChIKey
        inchi = MolToInchi(mol)
        if not inchi:
            return None

        inchiKey = InchiToInchiKey(inchi)
        if not inchiKey or len(inchiKey) < 14:
            return None

        # 返回 Block-1: connectivity hash (前 14 字符)
        return inchiKey[:14]

    except Exception:
        return None


def runTask1(df: pd.DataFrame) -> pd.DataFrame:
    """
    Task 1 主函数: 为全数据集生成 non_isomeric_inchikey_block1 列。
    Task 1 main: generate non_isomeric_inchikey_block1 column for all rows.
    """
    print("\n" + "=" * 70)
    print("  Task 1: 生成万能二维钥匙 / Generate 2D InChIKey Block-1")
    print("=" * 70)
    t0 = time.time()

    # 批量生成 (Batch generation)
    smilesList = df["canonical_smiles"].astype(str).tolist()
    results = []
    failCount = 0

    for i, smi in enumerate(smilesList):
        block1 = generateNonIsomericInchikeyBlock1(smi)
        results.append(block1)
        if block1 is None:
            failCount += 1

        if (i + 1) % 10000 == 0:
            print(f"    Processed {i+1:,}/{len(smilesList):,}...")

    df["non_isomeric_inchikey_block1"] = results

    # 也生成原始 InChIKey 的 Block-1 (方便对比)
    # Also extract original InChIKey Block-1 for comparison
    df["original_inchikey_block1"] = (
        df["standard_inchi_key"].astype(str).str[:14]
    )

    # 统计 (Statistics)
    validCount = df["non_isomeric_inchikey_block1"].notna().sum()
    uniqueBlock1 = df["non_isomeric_inchikey_block1"].nunique()

    # 测量降维效果: 手性异构体合并数
    # Measure dimensionality reduction: isomer merging
    originalUnique = df["original_inchikey_block1"].nunique()
    reduction = originalUnique - uniqueBlock1

    elapsed = time.time() - t0
    print(f"\n  Results:")
    print(f"    Valid 2D keys: {validCount:,} / {len(df):,} ({validCount/len(df)*100:.1f}%)")
    print(f"    Failed: {failCount:,}")
    print(f"    Original unique Block-1: {originalUnique:,}")
    print(f"    Non-isomeric unique Block-1: {uniqueBlock1:,}")
    print(f"    Dimensionality reduction: {reduction:,} isomers merged")
    print(f"    Time: {elapsed:.1f}s")

    return df


# =====================================================================
# Task 2: 离线知识库 Merge & 回填 (Offline KB Merge & Backfill)
# =====================================================================

def runTask2(
    df: pd.DataFrame,
    lotusPath: Optional[str] = None,
    coconutPath: Optional[str] = None,
    npatlasPath: Optional[str] = None,
) -> pd.DataFrame:
    """
    Task 2 主函数: 使用 non_isomeric_inchikey_block1 对接离线知识库。
    Task 2 main: merge with offline KBs using 2D InChIKey Block-1.

    回填策略 (Backfill strategy):
      - 仅填充 NaN/空值, 不覆盖已有数据 (.fillna() 语义)
      - 优先级: LOTUS > COCONUT > NPAtlas
    """
    print("\n" + "=" * 70)
    print("  Task 2: 离线知识库 Merge & 回填 / Offline KB Merge & Backfill")
    print("=" * 70)

    totalFilled = Counter()

    def _mergeAndBackfill(
        mainDf: pd.DataFrame,
        extPath: str,
        dbName: str,
        extKeyCol: str,
        fieldMapping: Dict[str, str],
        sep: str = ",",
    ) -> pd.DataFrame:
        """
        通用合并+回填函数。
        Generic merge+backfill function.

        Args:
            mainDf: 主数据集
            extPath: 外部文件路径
            dbName: 数据库名 (用于日志)
            extKeyCol: 外部文件中的 InChIKey/Block-1 列名
            fieldMapping: {外部列名: 主数据集目标列名}
            sep: 文件分隔符
        """
        if not os.path.exists(extPath):
            print(f"\n  [{dbName}] SKIP — File not found: {extPath}")
            return mainDf

        print(f"\n  [{dbName}] Loading: {extPath}")
        t0 = time.time()

        try:
            # 大文件分块读取 (Chunked reading for large files)
            if extPath.endswith(".gz") or os.path.getsize(extPath) > 500_000_000:
                extDf = pd.read_csv(extPath, sep=sep, dtype=str,
                                    low_memory=False, compression="infer")
            else:
                extDf = pd.read_csv(extPath, sep=sep, dtype=str,
                                    low_memory=False)
        except Exception as e:
            print(f"    [ERROR] Failed to read: {e}")
            return mainDf

        print(f"    External rows: {len(extDf):,}")

        if extKeyCol not in extDf.columns:
            # 尝试自动检测 InChIKey 列 (Auto-detect InChIKey column)
            inchikeyLikelyCols = [c for c in extDf.columns
                                  if "inchi" in c.lower() and "key" in c.lower()]
            if inchikeyLikelyCols:
                extKeyCol = inchikeyLikelyCols[0]
                print(f"    Auto-detected key column: '{extKeyCol}'")
            else:
                print(f"    [ERROR] Key column '{extKeyCol}' not found")
                print(f"    Available: {list(extDf.columns[:10])}")
                return mainDf

        # 生成外部数据的 Block-1 (Generate Block-1 for external data)
        extDf["_ext_block1"] = extDf[extKeyCol].astype(str).str[:14]
        extDf = extDf.dropna(subset=["_ext_block1"])
        extDf = extDf[extDf["_ext_block1"] != "nan"]

        # 去重: 每个 Block-1 只保留一条 (Dedup: keep first per Block-1)
        extDf = extDf.drop_duplicates(subset=["_ext_block1"])

        # Merge on Block-1 (2D match)
        merged = mainDf.merge(
            extDf[["_ext_block1"] + list(fieldMapping.keys())],
            left_on="non_isomeric_inchikey_block1",
            right_on="_ext_block1",
            how="left",
            suffixes=("", f"_{dbName}"),
        )

        # 回填 (Backfill) — 仅填充空值
        for extCol, mainCol in fieldMapping.items():
            srcCol = extCol if extCol in merged.columns else f"{extCol}_{dbName}"
            if srcCol not in merged.columns:
                continue

            if mainCol not in merged.columns:
                merged[mainCol] = ""

            # 识别空值 (Identify nulls in main column)
            isBlank = (
                merged[mainCol].isna() |
                (merged[mainCol].astype(str).isin(["", "nan", "None"]))
            )
            hasFill = (
                merged[srcCol].notna() &
                (~merged[srcCol].astype(str).isin(["", "nan", "None"]))
            )
            fillMask = isBlank & hasFill

            nFilled = fillMask.sum()
            if nFilled > 0:
                merged.loc[fillMask, mainCol] = merged.loc[fillMask, srcCol]
                totalFilled[mainCol] += nFilled
                print(f"    Backfilled {mainCol}: {nFilled:,} rows")

        # 清理临时列 (Clean temp columns)
        dropCols = [c for c in merged.columns
                    if c.startswith("_ext_") or c.endswith(f"_{dbName}")]
        merged.drop(columns=dropCols, inplace=True, errors="ignore")

        elapsed = time.time() - t0
        print(f"    [{dbName}] Time: {elapsed:.1f}s")

        return merged

    # ---- LOTUS ----
    if not lotusPath:
        lotusPath = os.path.join(
            os.path.dirname(__file__), "..", "data", "external",
            "lotus_frozen.csv.gz")
    df = _mergeAndBackfill(
        df, lotusPath, "LOTUS",
        extKeyCol="structure_inchikey",
        fieldMapping={
            "organism_name": "organisms",
            "reference_doi": "dois",
            "structure_nameTraditional": "name",
        },
    )

    # ---- COCONUT ----
    if not coconutPath:
        coconutPath = os.path.join(
            os.path.dirname(__file__), "..", "data", "external",
            "COCONUT_DB.csv")
    df = _mergeAndBackfill(
        df, coconutPath, "COCONUT",
        extKeyCol="standard_inchi_key",
        fieldMapping={
            "name": "name",
            "synonyms": "synonyms",
            "organisms": "organisms",
            "dois": "dois",
        },
    )

    # ---- NPAtlas ----
    if not npatlasPath:
        npatlasPath = os.path.join(
            os.path.dirname(__file__), "..", "data", "external",
            "npatlas_download.tsv")
    df = _mergeAndBackfill(
        df, npatlasPath, "NPAtlas",
        extKeyCol="compound_inchikey",
        fieldMapping={
            "compound_names": "name",
            "genus": "organisms",
            "doi": "dois",
        },
        sep="\t",
    )

    # 汇总 (Summary)
    print(f"\n  Task 2 Total Backfill:")
    for col, count in totalFilled.most_common():
        print(f"    {col}: {count:,} rows filled")
    if not totalFilled:
        print(f"    No offline databases found. Download and place in data/external/")

    return df


# =====================================================================
# Task 3: PubChem CID 名称反向兜底 (PubChem Name Fallback)
# =====================================================================

def runTask3(
    df: pd.DataFrame,
    limit: int = 200,
    enablePubchem: bool = False,
    localSynonymPath: Optional[str] = None,
) -> pd.DataFrame:
    """
    Task 3: 通过名称从 PubChem 兜底获取缺失的 SMILES/DOI。
    Task 3: Fallback PubChem lookup by name for missing data.

    两条路径:
      A) 本地离线: PubChem CID-Synonym 映射表 (如果用户已下载)
      B) 在线查询: pubchempy API (需要网络, 限流 200 条)

    安全原则: 仅在 enablePubchem=True 时启用在线查询。
    """
    print("\n" + "=" * 70)
    print("  Task 3: PubChem 名称兜底 / PubChem Name Fallback")
    print("=" * 70)

    # 识别需要兜底的行: 有 name 但 DOI 仍缺失
    # Identify rows needing fallback: has name but DOI still missing
    hasName = (
        df["name"].notna() &
        (~df["name"].astype(str).isin(["", "nan", "None"]))
    )
    needsDoi = (
        df["dois"].isna() |
        df["dois"].astype(str).isin(["", "nan", "None"])
    )
    candidates = df[hasName & needsDoi]
    print(f"  Candidates (has name, needs DOI): {len(candidates):,}")

    if len(candidates) == 0:
        print(f"  All DOIs already filled!")
        return df

    # ---- Path A: 本地离线 CID-Synonym 映射表 ----
    if localSynonymPath and os.path.exists(localSynonymPath):
        print(f"\n  [Path A] Loading local PubChem synonym table...")
        print(f"    File: {localSynonymPath}")
        try:
            synDf = pd.read_csv(localSynonymPath, sep="\t", dtype=str,
                                low_memory=False)
            print(f"    Loaded: {len(synDf):,} rows")

            # 构建 name → CID 映射 (Build name → CID mapping)
            if "Synonym" in synDf.columns and "CID" in synDf.columns:
                synDf["_name_lower"] = synDf["Synonym"].str.lower().str.strip()
                nameMap = synDf.set_index("_name_lower")["CID"].to_dict()

                filled = 0
                for idx in candidates.index:
                    nameStr = str(df.at[idx, "name"]).lower().strip()
                    if nameStr in nameMap:
                        cid = nameMap[nameStr]
                        df.at[idx, "PubChem_CID"] = cid
                        filled += 1

                print(f"    Local synonym match: {filled:,} / {len(candidates):,}")
        except Exception as e:
            print(f"    [ERROR] Failed: {e}")

    # ---- Path B: pubchempy 在线查询 ----
    if enablePubchem:
        print(f"\n  [Path B] PubChem online lookup via pubchempy...")
        print(f"    Limit: {limit} queries")

        try:
            import pubchempy as pcp
        except ImportError:
            print(f"    [ERROR] pubchempy not installed: pip install pubchempy")
            return df

        # 重新计算候选行 (refresh after Path A)
        needsDoi = (
            df["dois"].isna() |
            df["dois"].astype(str).isin(["", "nan", "None"])
        )
        candidates = df[hasName & needsDoi]

        # 抽样 (Sample to respect limit)
        sampleIdxs = candidates.head(limit).index

        filled = 0
        errors = 0
        t0 = time.time()

        for i, idx in enumerate(sampleIdxs):
            nameStr = str(df.at[idx, "name"]).strip()
            if not nameStr or nameStr == "nan":
                continue

            try:
                results = pcp.get_compounds(nameStr, "name")
                if results:
                    compound = results[0]
                    # 回填 CID
                    df.at[idx, "PubChem_CID"] = str(compound.cid)

                    # 回填 InChIKey (如果更权威)
                    if compound.inchikey:
                        currentIk = str(df.at[idx, "standard_inchi_key"])
                        if currentIk in ("", "nan", "None"):
                            df.at[idx, "standard_inchi_key"] = compound.inchikey

                    filled += 1
                else:
                    errors += 1

            except Exception:
                errors += 1

            # 限流: 每次查询间隔 200ms (Rate limit)
            time.sleep(0.2)

            if (i + 1) % 50 == 0:
                print(f"    Progress: {i+1}/{len(sampleIdxs)} "
                      f"(filled={filled}, errors={errors})")

        elapsed = time.time() - t0
        print(f"\n  [PubChem] Results:")
        print(f"    Queried: {len(sampleIdxs):,}")
        print(f"    Filled: {filled:,}")
        print(f"    Errors: {errors:,}")
        print(f"    Time: {elapsed:.0f}s")
    else:
        print(f"\n  [PubChem online] DISABLED (use --enable-pubchem to activate)")
        print(f"  [提示] 建议先下载 PubChem 离线 CID-Synonym 文件:")
        print(f"    https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-Synonym-filtered.gz")
        print(f"    解压后放入 data/external/, 使用 --local-synonyms 指定路径")

    return df


# =====================================================================
# 统计报告 (Final Statistics Report)
# =====================================================================

def printBeforeAfterStats(
    beforeStats: Dict[str, int],
    afterDf: pd.DataFrame,
    total: int,
):
    """打印补全前后对比。 / Print before/after comparison."""
    print("\n" + "=" * 70)
    print("  补全效果对比 / Imputation Results Comparison")
    print("=" * 70)

    fields = ["standard_inchi_key", "name", "dois", "organisms"]
    fieldLabels = ["InChIKey", "Name", "DOI", "Organisms"]

    print(f"\n  {'Field':<15} {'Before':>12} {'After':>12} {'Gained':>10} {'Rate':>8}")
    print(f"  {'-'*15} {'-'*12} {'-'*12} {'-'*10} {'-'*8}")

    for field, label in zip(fields, fieldLabels):
        before = beforeStats.get(field, 0)

        afterValid = afterDf[field].notna() & (
            ~afterDf[field].astype(str).isin(["", "nan", "None"])
        )
        after = afterValid.sum()
        gained = after - before
        rate = after / total * 100

        gainStr = f"+{gained:,}" if gained >= 0 else f"{gained:,}"
        print(f"  {label:<15} {before:>12,} {after:>12,} {gainStr:>10} {rate:>7.1f}%")

    # 新列统计 (New column stats)
    newCols = ["non_isomeric_inchikey_block1", "PubChem_CID"]
    for col in newCols:
        if col in afterDf.columns:
            valid = afterDf[col].notna() & (
                ~afterDf[col].astype(str).isin(["", "nan", "None"])
            )
            print(f"  {col:<15} {'(new)':>12} {valid.sum():>12,} {'':>10} "
                  f"{valid.sum()/total*100:>7.1f}%")


# =====================================================================
# Main
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Dimensionality Reduction Imputation Engine")
    parser.add_argument("--input", type=str, default=None)
    parser.add_argument("--output", type=str, default=None)
    parser.add_argument("--task1-only", action="store_true",
                        help="Only generate 2D keys, skip merge and PubChem")
    parser.add_argument("--lotus", type=str, default=None)
    parser.add_argument("--coconut", type=str, default=None)
    parser.add_argument("--npatlas", type=str, default=None)
    parser.add_argument("--enable-pubchem", action="store_true",
                        help="Enable PubChem online lookup")
    parser.add_argument("--pubchem-limit", type=int, default=200,
                        help="Max PubChem online queries")
    parser.add_argument("--local-synonyms", type=str, default=None,
                        help="Path to PubChem CID-Synonym local file")
    args = parser.parse_args()

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")

    inputPath = args.input or os.path.join(
        reportDir, "GlycoNP_Pipeline_Full_Cleaned.csv")
    outputPath = args.output or os.path.join(
        reportDir, "GlycoNP_Imputed.csv")

    print("=" * 70)
    print("  GlycoNP Dimensionality Reduction Imputation Engine")
    print("  糖缀合物本地数据降维补全引擎")
    print("=" * 70)
    print(f"  Input:  {inputPath}")
    print(f"  Output: {outputPath}")

    if not os.path.exists(inputPath):
        print(f"\n  [ERROR] Input not found: {inputPath}")
        return

    t0 = time.time()
    df = pd.read_csv(inputPath, low_memory=False, dtype=str,
                     encoding="utf-8-sig")
    total = len(df)
    print(f"  Loaded: {total:,} rows")

    # 记录补全前的统计 (Record pre-imputation stats)
    beforeStats = {}
    for field in ["standard_inchi_key", "name", "dois", "organisms"]:
        validMask = df[field].notna() & (
            ~df[field].astype(str).isin(["", "nan", "None"])
        )
        beforeStats[field] = validMask.sum()

    # ============ Task 1 ============
    df = runTask1(df)

    if args.task1_only:
        df.to_csv(outputPath, index=False, encoding="utf-8-sig")
        printBeforeAfterStats(beforeStats, df, total)
        print(f"\n  Output: {outputPath}")
        return

    # ============ Task 2 ============
    df = runTask2(
        df,
        lotusPath=args.lotus,
        coconutPath=args.coconut,
        npatlasPath=args.npatlas,
    )

    # ============ Task 3 ============
    df = runTask3(
        df,
        limit=args.pubchem_limit,
        enablePubchem=args.enable_pubchem,
        localSynonymPath=args.local_synonyms,
    )

    # ---- 保存 (Save) ----
    df.to_csv(outputPath, index=False, encoding="utf-8-sig")
    elapsed = time.time() - t0

    # ---- 最终统计 (Final stats) ----
    printBeforeAfterStats(beforeStats, df, total)

    print(f"\n  Output: {outputPath}")
    print(f"  Total time: {elapsed:.0f}s")
    print("=" * 70)
    print(f"\n  下一步 / Next Steps:")
    print(f"  1. 下载离线数据库放入 data/external/, 重跑本脚本")
    print(f"  2. 用补全后的数据运行: python scripts/mine_chembl_sqlite.py --input {outputPath}")
    print(f"  3. 补全与活性数据融合后, 运行: python scripts/aggregate_preference.py")


if __name__ == "__main__":
    main()
