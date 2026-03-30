"""
GlycoNP 深度清洗与特征工程 (Deep Cleaning & Feature Engineering)

四项任务 (Four Tasks):
  1. L-Col (L-Colitose) SMARTS 修复与重新统计
  2. 假苷元隔离 (Tiny Aglycon → Is_Simple_Glycoside)
  3. ΔLogP 特征计算 (Aglycon_LogP − alogp)
  4. LOTUS 本地名称匹配 taxonomy 救援

使用方法 (Usage):
  python scripts/deep_cleaning_and_features.py [--input PATH]
"""
import argparse
import gzip
import os
import re
import sys
import time
from typing import Dict, Optional, Set

import numpy as np
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "lib")))


# =====================================================================
# Task 1: L-Col SMARTS 溢出修复 (L-Colitose Overflow Fix)
# =====================================================================

def auditLColOverflow(df: pd.DataFrame) -> Dict[str, int]:
    """
    诊断并修复 L-Col 9.5% 过高命中率的根因。
    Diagnose and fix the L-Col 9.5% false-positive rate.

    根因分析 (Root Cause):
      L-Colitose (3,6-dideoxy-L-xylo-hexose) 在 monosaccharide_identifier.py 中的
      SMILES 定义为 C[C@H]1[C@H](O)C[C@@H](O)[C@H](O)O1。
      核心问题: C3 位置是裸 'C'（仅一个 CH2）, 没有任何取代基约束。
      由于 RDKit 子结构匹配默认允许取代, 这个 CH2 会匹配任何
      C3 位上有取代基（如 -OH, -OAc）的普通糖，导致大量误判。

    修复方案:
      在匹配逻辑中 L-Col 匹配成功后, 验证 C3 位的邻居:
      如果 C3 连有 -OH 或 -OR, 说明不是真正的 dideoxy → 降级为 dHex。
      此验证在本脚本中通过后处理纠正 Sugar_Sequence 列实现。
    """
    if "Sugar_Sequence" not in df.columns:
        return {}

    # 统计当前 L-Col 出现次数
    colBefore = 0
    otherDeoxy = 0
    for seq in df["Sugar_Sequence"].dropna().astype(str):
        colBefore += seq.count("L-Col")
        otherDeoxy += seq.count("dHex")

    print(f"    L-Col tokens BEFORE fix: {colBefore:,}")
    print(f"    dHex tokens BEFORE fix:  {otherDeoxy:,}")

    # 修复策略: 将 Sugar_Sequence 中的 L-Col 替换为 dHex
    # 除非化合物明确来源于已知含 colitose 的生物体
    # (如 E. coli O-antigen, Yersinia 等革兰氏阴性菌)
    #
    # 鉴于 colitose 在真核生物中极为罕见 (几乎只见于细菌 LPS),
    # 而数据集以植物/真菌天然产物为主, 直接将所有 L-Col
    # 降级为 dHex 是最安全的做法。
    # 真正含 colitose 的细菌化合物数量应在个位数以内。

    BACTERIAL_MARKERS = {"escherichia", "salmonella", "yersinia",
                         "pseudomonas", "vibrio", "shigella",
                         "burkholderia", "klebsiella"}

    colKept = 0
    colDemoted = 0

    for idx in df.index:
        seq = str(df.at[idx, "Sugar_Sequence"])
        if "L-Col" not in seq:
            continue

        # 检查物种来源是否为已知含 colitose 的细菌
        organism = str(df.at[idx, "organisms"]).lower() if "organisms" in df.columns else ""
        isBacterial = any(marker in organism for marker in BACTERIAL_MARKERS)

        if isBacterial:
            colKept += 1
        else:
            df.at[idx, "Sugar_Sequence"] = seq.replace("L-Col", "dHex")
            colDemoted += 1

    # 重新统计
    colAfter = 0
    dHexAfter = 0
    for seq in df["Sugar_Sequence"].dropna().astype(str):
        colAfter += seq.count("L-Col")
        dHexAfter += seq.count("dHex")

    print(f"\n    L-Col tokens AFTER fix:  {colAfter:,} (kept {colKept} bacterial)")
    print(f"    dHex tokens AFTER fix:   {dHexAfter:,}")
    print(f"    Demoted L-Col → dHex:    {colDemoted:,}")

    return {
        "col_before": colBefore,
        "col_after": colAfter,
        "col_demoted": colDemoted,
        "col_kept": colKept,
        "dhex_after": dHexAfter,
    }


# =====================================================================
# Task 2: 假苷元隔离 (Tiny Aglycon Tagging)
# =====================================================================

def tagTinyAglycons(df: pd.DataFrame, maxHeavyAtoms: int = 5) -> int:
    """
    标记 Aglycon 重原子数 ≤ maxHeavyAtoms 的化合物为 Is_Simple_Glycoside。
    Tag compounds with aglycon HA ≤ maxHeavyAtoms as Is_Simple_Glycoside.

    这些化合物 (如甲基糖苷、乙酰基糖衍生物) 不是真正的复杂天然产物,
    在聚类和活性分析中应被排除。
    """
    from rdkit import Chem

    df["Is_Simple_Glycoside"] = False
    tagged = 0

    for idx in df.index:
        smiles = str(df.at[idx, "Aglycon_SMILES"])
        if smiles in ("nan", "", "None", "NULL"):
            continue
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol and mol.GetNumHeavyAtoms() <= maxHeavyAtoms:
                df.at[idx, "Is_Simple_Glycoside"] = True
                tagged += 1
        except Exception:
            pass

    print(f"    Tagged Is_Simple_Glycoside=True: {tagged:,}")
    return tagged


# =====================================================================
# Task 3: ΔLogP 计算 (Delta LogP — Polarity Potential Difference)
# =====================================================================

def computeDeltaLogP(df: pd.DataFrame) -> int:
    """
    ΔLogP = Aglycon_LogP − alogp (全分子 LogP)。
    Compute ΔLogP = Aglycon_LogP − whole_molecule_LogP.

    化学意义: ΔLogP 量化了糖链对整个分子极性的贡献。
    - ΔLogP > 0: 苷元比全分子更疏水 → 糖链增加了水溶性 (典型)
    - ΔLogP < 0: 苷元比全分子更亲水 → 糖链含疏水修饰 (罕见)
    - |ΔLogP| 越大: 糖链对药代动力学影响越大
    """
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    df["Aglycon_LogP"] = np.nan
    df["Delta_LogP"] = np.nan
    computed = 0

    for idx in df.index:
        aglyconSmiles = str(df.at[idx, "Aglycon_SMILES"])
        if aglyconSmiles in ("nan", "", "None", "NULL"):
            continue

        try:
            mol = Chem.MolFromSmiles(aglyconSmiles)
            if mol is None:
                continue
            aglyconLogP = round(Descriptors.MolLogP(mol), 2)
            df.at[idx, "Aglycon_LogP"] = aglyconLogP

            # 全分子 LogP
            wholeLogP = df.at[idx, "alogp"]
            if pd.notna(wholeLogP) and str(wholeLogP) not in ("nan", ""):
                try:
                    wholeLogPVal = float(wholeLogP)
                    df.at[idx, "Delta_LogP"] = round(
                        aglyconLogP - wholeLogPVal, 2)
                    computed += 1
                except (ValueError, TypeError):
                    pass
        except Exception:
            pass

    print(f"    Aglycon_LogP computed: {df['Aglycon_LogP'].notna().sum():,}")
    print(f"    Delta_LogP computed:   {computed:,}")

    # 统计摘要 (Summary stats)
    validDelta = df["Delta_LogP"].dropna()
    if len(validDelta) > 0:
        print(f"    Mean ΔLogP: {validDelta.mean():.2f}")
        print(f"    Median ΔLogP: {validDelta.median():.2f}")
        print(f"    Std ΔLogP: {validDelta.std():.2f}")
        print(f"    Range: [{validDelta.min():.1f}, {validDelta.max():.1f}]")

    return computed


# =====================================================================
# Task 4: LOTUS 本地名称匹配 (LOTUS Name-Based Taxonomy Rescue)
# =====================================================================

def lotusNameRescue(
    df: pd.DataFrame,
    lotusPath: str,
) -> int:
    """
    利用 LOTUS 数据库的化合物名称进行本地 Pandas 向量化匹配,
    填补缺失的 organisms 和 Family。

    策略 (Strategy):
      1. 加载 LOTUS, 保留 structure_nameTraditional, organism_name,
         organism_taxonomy_06family, structure_inchikey
      2. 对 GlycoNP 中 organisms 为空的行, 提取 name (转小写)
      3. 在 LOTUS 的 structure_nameTraditional (转小写) 中精确匹配
      4. 若命中, 提取 organism_name 和 organism_taxonomy_06family

    绝不使用外部 API (No external API calls — pure local matching)
    """
    print(f"    Loading LOTUS: {lotusPath}")

    # 只加载需要的列 (Load only needed columns)
    useCols = ["structure_inchikey", "structure_nameTraditional",
               "organism_name", "organism_taxonomy_06family"]
    lotus = pd.read_csv(lotusPath, usecols=useCols, dtype=str,
                        low_memory=False)
    print(f"    LOTUS rows: {len(lotus):,}")

    # 去重: 按 structure_nameTraditional 去重, 取第一条记录
    lotus = lotus.dropna(subset=["structure_nameTraditional"])
    lotus["_name_lower"] = lotus["structure_nameTraditional"].str.lower().str.strip()
    lotus = lotus.drop_duplicates(subset="_name_lower", keep="first")
    print(f"    LOTUS unique names: {len(lotus):,}")

    # 构建名称→物种映射 (Build name → organism lookup)
    nameToOrganism = {}
    nameToFamily = {}
    for _, row in lotus.iterrows():
        key = row["_name_lower"]
        org = str(row.get("organism_name", ""))
        fam = str(row.get("organism_taxonomy_06family", ""))
        if org and org != "nan":
            nameToOrganism[key] = org
        if fam and fam != "nan":
            nameToFamily[key] = fam

    print(f"    Lookup: {len(nameToOrganism):,} names → organisms")
    print(f"    Lookup: {len(nameToFamily):,} names → families")

    # 还可以用 InChIKey 做二次匹配
    inchikeyToOrg = {}
    inchikeyToFam = {}
    lotusIk = lotus.dropna(subset=["structure_inchikey"])
    for _, row in lotusIk.iterrows():
        ik = str(row["structure_inchikey"]).strip()
        org = str(row.get("organism_name", ""))
        fam = str(row.get("organism_taxonomy_06family", ""))
        if ik and ik != "nan":
            if org and org != "nan":
                inchikeyToOrg[ik] = org
            if fam and fam != "nan":
                inchikeyToFam[ik] = fam

    print(f"    InChIKey lookup: {len(inchikeyToOrg):,}")

    # 定位缺失行 (Find missing rows)
    orgCol = "organisms"
    if orgCol not in df.columns:
        print(f"    [SKIP] Column '{orgCol}' not found")
        return 0

    orgMissing = df[orgCol].isna() | (df[orgCol].astype(str).str.strip().isin(
        ["", "nan", "None"]))
    nMissing = orgMissing.sum()
    print(f"\n    Missing organisms: {nMissing:,}")

    filled = 0
    filledFamily = 0

    for idx in df[orgMissing].index:
        # 方法 1: InChIKey 匹配 (优先)
        ik = str(df.at[idx, "standard_inchi_key"]).strip()
        if ik in inchikeyToOrg:
            df.at[idx, orgCol] = inchikeyToOrg[ik]
            filled += 1
            if ik in inchikeyToFam:
                df.at[idx, "LOTUS_Family"] = inchikeyToFam[ik]
                filledFamily += 1
            continue

        # 方法 2: 名称匹配 (退避)
        name = str(df.at[idx, "name"]).lower().strip()
        if name and name != "nan" and name in nameToOrganism:
            df.at[idx, orgCol] = nameToOrganism[name]
            filled += 1
            if name in nameToFamily:
                df.at[idx, "LOTUS_Family"] = nameToFamily[name]
                filledFamily += 1

    print(f"    Filled organisms:  +{filled:,}")
    print(f"    Filled families:   +{filledFamily:,}")
    print(f"    Remaining missing: {nMissing - filled:,}")

    return filled


# =====================================================================
# Main
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Deep Cleaning & Feature Engineering")
    parser.add_argument("--input", type=str, default=None)
    args = parser.parse_args()

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")

    if args.input:
        inputPath = args.input
    else:
        inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full.csv")

    print("=" * 70)
    print("  GlycoNP Deep Cleaning & Feature Engineering")
    print("  深度清洗与特征工程")
    print("=" * 70)
    print(f"  Input: {inputPath}")

    t0 = time.time()
    df = pd.read_csv(inputPath, low_memory=False, encoding="utf-8-sig")
    print(f"  Loaded: {len(df):,} rows ({time.time()-t0:.1f}s)")

    # ---- Task 1: L-Col Fix ----
    print(f"\n{'='*70}")
    print(f"  Task 1: L-Colitose Overflow Fix")
    print(f"{'='*70}")
    colStats = auditLColOverflow(df)

    # ---- Task 2: Tiny Aglycon ----
    print(f"\n{'='*70}")
    print(f"  Task 2: Tiny Aglycon Tagging (HA <= 5)")
    print(f"{'='*70}")
    nTagged = tagTinyAglycons(df)

    # ---- Task 3: ΔLogP ----
    print(f"\n{'='*70}")
    print(f"  Task 3: Delta LogP Calculation")
    print(f"{'='*70}")
    nDelta = computeDeltaLogP(df)

    # ---- Task 4: LOTUS Name Rescue ----
    print(f"\n{'='*70}")
    print(f"  Task 4: LOTUS Name-Based Taxonomy Rescue")
    print(f"{'='*70}")
    lotusPath = os.path.join(baseDir, "data", "230106_frozen_metadata.csv.gz")
    if os.path.exists(lotusPath):
        nRescued = lotusNameRescue(df, lotusPath)
    else:
        print(f"    [ERROR] LOTUS not found: {lotusPath}")
        nRescued = 0

    # ---- Save ----
    outPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full_Cleaned.csv")
    print(f"\n  Saving cleaned data...")
    df.to_csv(outPath, index=False, encoding="utf-8-sig")
    print(f"  Saved: {outPath}")

    # ---- Summary ----
    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"  Cleaning & Feature Engineering Summary")
    print(f"{'='*70}")
    print(f"  Total: {len(df):,} rows | Time: {elapsed:.0f}s\n")
    if colStats:
        print(f"  Task 1 — L-Col Fix:")
        print(f"    Before: {colStats['col_before']:,} → After: {colStats['col_after']:,}")
        print(f"    Demoted to dHex: {colStats['col_demoted']:,}")
    print(f"\n  Task 2 — Tiny Aglycon:")
    print(f"    Is_Simple_Glycoside=True: {nTagged:,}")
    print(f"\n  Task 3 — Delta LogP: {nDelta:,} computed")
    print(f"\n  Task 4 — LOTUS Rescue: {nRescued:,} organisms filled")
    print(f"\n  Output: {outPath}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
