"""
立体化学完整度比对与局部修补 (Stereo-Completeness Patch)
=======================================================

重构 Task 1: 用 RDKit FindMolChiralCenters 比较手性中心绝对数量,
而非简单判断 SMILES 是否含 @ 符号。

规则: LOTUS_Chiral_Count > Original_Chiral_Count → 无条件替换

随后立即重跑 Task 2 名称糖链修复, 看残余 Hex/Non/Pen 是否进一步消灭。

Usage:
  python scripts/stereo_completeness_patch.py
"""
import os
import re
import sys
import time
from collections import Counter
from typing import Optional

import pandas as pd
from rdkit import Chem

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
REPORT_DIR = os.path.join(BASE_DIR, "reports")
LOTUS_GZ = os.path.join(BASE_DIR, "data", "230106_frozen_metadata.csv.gz")
INPUT_CSV = os.path.join(REPORT_DIR, "GlycoNP_Fully_Enriched.csv")


# =====================================================================
# 辅助函数: 手性中心绝对计数 (Chiral Center Counting)
# =====================================================================

def countDefinedChiralCenters(smiles: str) -> int:
    """
    计算 SMILES 中已明确定义的四面体手性中心数量。
    Count the number of explicitly defined tetrahedral chiral centers.

    使用 FindMolChiralCenters(includeUnassigned=False) — 只计入
    有明确 R/S 标注的中心, 忽略未定义的潜在手性位点。

    设计意图: 区分 "完全无手性" 和 "局部手性缺失" 两种情况。
    一个分子可能有 12 个手性中心但 SMILES 只定义了 2 个 — 这就是
    "partial stereochemistry", 我们需要用 LOTUS 的完全版本替换它。
    """
    if not smiles or str(smiles) in ("nan", "", "None"):
        return 0
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return 0
        # includeUnassigned=False: 只计已定义的手性中心
        # force=True: 确保重新计算
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        chiralCenters = Chem.FindMolChiralCenters(
            mol, includeUnassigned=False, useLegacyImplementation=False)
        return len(chiralCenters)
    except Exception:
        return 0


def countTotalPotentialChiralCenters(smiles: str) -> int:
    """
    计算分子中所有潜在的手性中心数量 (包括未定义的)。
    Count total potential chiral centers (including unassigned).

    用于计算 "手性完整度" = defined / total
    """
    if not smiles or str(smiles) in ("nan", "", "None"):
        return 0
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return 0
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        allCenters = Chem.FindMolChiralCenters(
            mol, includeUnassigned=True, useLegacyImplementation=False)
        return len(allCenters)
    except Exception:
        return 0


# =====================================================================
# Task 1 重构: 立体结构完整度比对 (Stereo-Completeness Check)
# =====================================================================

def runStereoCompletnessPatch(df: pd.DataFrame) -> pd.DataFrame:
    """
    核心逻辑:
      1. 加载 LOTUS Block-1 → structure_smiles 映射
      2. 对每个 Block-1 匹配的分子对:
         - 计算 original chiral count
         - 计算 LOTUS chiral count
         - 若 LOTUS > original → 替换
    """
    print("\n" + "=" * 70)
    print("  Stereo-Completeness Patch")
    print("  立体化学完整度比对与局部修补")
    print("=" * 70)
    t0 = time.time()

    if not os.path.exists(LOTUS_GZ):
        print(f"  [ERROR] LOTUS not found: {LOTUS_GZ}")
        return df

    # 加载 LOTUS SMILES (Load LOTUS SMILES)
    print("  [1/4] Loading LOTUS SMILES...")
    lotusDf = pd.read_csv(
        LOTUS_GZ, dtype=str, low_memory=False,
        usecols=["structure_inchikey", "structure_smiles"])
    lotusDf["_block1"] = lotusDf["structure_inchikey"].astype(str).str[:14]
    lotusDf = lotusDf.dropna(subset=["_block1", "structure_smiles"])
    lotusDf = lotusDf[
        (lotusDf["_block1"] != "nan") &
        (lotusDf["structure_smiles"] != "nan") &
        (lotusDf["structure_smiles"] != "")
    ]
    # 去重: 每个 Block-1 保留 SMILES 最长的记录 (通常信息最丰富)
    # Dedup: keep longest SMILES per Block-1 (usually most informative)
    lotusDf["_smi_len"] = lotusDf["structure_smiles"].str.len()
    lotusDf = lotusDf.sort_values("_smi_len", ascending=False)
    lotusDf = lotusDf.drop_duplicates(subset=["_block1"])
    lotusMap = dict(zip(lotusDf["_block1"], lotusDf["structure_smiles"]))
    print(f"    LOTUS lookup table: {len(lotusMap):,} entries")

    # 构建候选集 (Build candidate set)
    print("  [2/4] Finding matching Block-1 candidates...")
    hasBlock1 = df["non_isomeric_inchikey_block1"].notna()
    matchedIdxs = []
    for idx in df.index[hasBlock1]:
        block1 = str(df.at[idx, "non_isomeric_inchikey_block1"])
        if block1 in lotusMap:
            matchedIdxs.append(idx)
    print(f"    Block-1 matches: {len(matchedIdxs):,}")

    # 逐条比对手性中心数量 (Compare chiral center counts)
    print("  [3/4] Comparing chiral center counts...")
    upgraded = 0
    partialCount = 0   # 有部分手性但不完整
    alreadyFull = 0    # 已经完整, 无需替换
    lotusWorse = 0     # LOTUS 竟然更差
    sameCounts = 0     # 相等

    if "Is_Stereo_Upgraded" not in df.columns:
        df["Is_Stereo_Upgraded"] = False

    gainDistribution = Counter()  # 统计手性中心增量分布 (gain distribution)

    for i, idx in enumerate(matchedIdxs):
        origSmiles = str(df.at[idx, "canonical_smiles"])
        block1 = str(df.at[idx, "non_isomeric_inchikey_block1"])
        lotusSmiles = lotusMap[block1]

        origChiral = countDefinedChiralCenters(origSmiles)
        lotusChiral = countDefinedChiralCenters(lotusSmiles)

        if lotusChiral > origChiral:
            # LOTUS 有更多已定义手性中心 → 无条件替换
            if "original_smiles_2d" not in df.columns:
                df["original_smiles_2d"] = ""
            if not df.at[idx, "original_smiles_2d"] or \
               str(df.at[idx, "original_smiles_2d"]) in ("", "nan", "None"):
                df.at[idx, "original_smiles_2d"] = origSmiles
            df.at[idx, "canonical_smiles"] = lotusSmiles
            df.at[idx, "Is_Stereo_Upgraded"] = True

            gain = lotusChiral - origChiral
            gainDistribution[gain] += 1
            upgraded += 1

            if origChiral > 0:
                partialCount += 1
        elif lotusChiral == origChiral:
            sameCounts += 1
            if origChiral > 0:
                alreadyFull += 1
        else:
            lotusWorse += 1

        if (i + 1) % 10000 == 0:
            print(f"    Compared {i+1:,}/{len(matchedIdxs):,}... "
                  f"(upgraded={upgraded:,})")

    elapsed = time.time() - t0
    print(f"\n  [4/4] Stereo-Completeness Patch Results:")
    print(f"    Total Block-1 matches compared: {len(matchedIdxs):,}")
    print(f"    [OK] Upgraded (LOTUS > Original): {upgraded:,}")
    print(f"       - Of which partial stereo fixed: {partialCount:,}")
    print(f"    [==] Same (equal counts): {sameCounts:,}")
    print(f"       - Already fully defined: {alreadyFull:,}")
    print(f"    [--] LOTUS worse: {lotusWorse:,}")
    print(f"    Time: {elapsed:.0f}s")

    if gainDistribution:
        print(f"\n    Chiral center gain distribution:")
        for gain, count in sorted(gainDistribution.items()):
            print(f"      +{gain} centers: {count:,} compounds")

    return df


# =====================================================================
# Task 2 快速重跑: 名称糖链修复 (Quick Sugar Re-parse)
# =====================================================================

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


def quickSugarReparse(df: pd.DataFrame) -> pd.DataFrame:
    """快速重跑名称糖链修复, 专门针对残余的 31 个泛指标签。"""
    print("\n" + "=" * 70)
    print("  Quick Sugar Re-parse (Post Stereo Patch)")
    print("  糖链快速修复 (立体结构补全后)")
    print("=" * 70)

    seqCol = "Sugar_Sequence"
    nameColCandidates = ["name", "iupac_name", "synonyms"]

    # 统计当前残留 (Count current residuals)
    allSeqs = df[seqCol].dropna().astype(str)
    beforeTokens = []
    for seq in allSeqs:
        beforeTokens.extend(GENERIC_PATTERN.findall(seq))
    beforeCounter = Counter(beforeTokens)
    totalBefore = sum(beforeCounter.values())

    print(f"  Before re-parse:")
    for k, v in beforeCounter.most_common():
        print(f"    {k:>6s}: {v:,}")
    print(f"    Total: {totalBefore:,}")

    # 执行修复 (Execute rescue)
    rescued = 0
    for idx in df.index:
        seq = str(df.at[idx, seqCol]) if pd.notna(df.at[idx, seqCol]) else ""
        if not GENERIC_PATTERN.search(seq):
            continue

        namePool = ""
        for col in nameColCandidates:
            if col in df.columns and pd.notna(df.at[idx, col]):
                val = str(df.at[idx, col])
                if val not in ("", "nan", "None"):
                    namePool += " " + val

        if not namePool.strip():
            continue

        def replaceMatch(m):
            token = m.group(0)
            if token not in ("Hex", "Pen", "dHex", "HexA", "Non", "Oct", "Hept"):
                return token
            for regex, sugar in COMPILED_TEXT_RULES:
                if regex.search(namePool):
                    return sugar
            return token

        newSeq = GENERIC_PATTERN.sub(replaceMatch, seq)
        if newSeq != seq:
            df.at[idx, seqCol] = newSeq
            rescued += 1

    # 统计修复后残留 (Count after)
    allSeqs = df[seqCol].dropna().astype(str)
    afterTokens = []
    for seq in allSeqs:
        afterTokens.extend(GENERIC_PATTERN.findall(seq))
    afterCounter = Counter(afterTokens)
    totalAfter = sum(afterCounter.values())
    resolved = totalBefore - totalAfter

    print(f"\n  After re-parse:")
    if afterCounter:
        for k, v in afterCounter.most_common():
            print(f"    {k:>6s}: {v:,}")
        print(f"    Total: {totalAfter:,}")
    else:
        print(f"    >>> ALL GENERIC LABELS ELIMINATED! <<<")
    print(f"  Resolved this pass: {resolved:,}")
    print(f"  Rows modified: {rescued:,}")

    return df


# =====================================================================
# Main
# =====================================================================

def main():
    print("=" * 70)
    print("  Stereo-Completeness Patch + Sugar Re-parse")
    print("  立体化学完整度比对 + 糖链快速修复")
    print("=" * 70)

    if not os.path.exists(INPUT_CSV):
        print(f"  [ERROR] Input not found: {INPUT_CSV}")
        return

    t0 = time.time()
    df = pd.read_csv(INPUT_CSV, dtype=str, low_memory=False,
                     encoding="utf-8-sig")
    total = len(df)
    print(f"  Loaded: {total:,} rows")

    # Step 1: Stereo patch
    df = runStereoCompletnessPatch(df)

    # Step 2: Sugar re-parse
    df = quickSugarReparse(df)

    # Save
    df.to_csv(INPUT_CSV, index=False, encoding="utf-8-sig")
    elapsed = time.time() - t0
    print(f"\n  Saved: {INPUT_CSV}")
    print(f"  Total time: {elapsed:.0f}s")
    print("=" * 70)


if __name__ == "__main__":
    main()
