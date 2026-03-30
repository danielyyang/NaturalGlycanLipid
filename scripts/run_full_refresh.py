"""
GlycoNP 全量糖序列刷新管线 (Full Sugar Sequence Refresh Pipeline)
=================================================================

在现有 GlycoNP_Deep_Enriched_v2.csv 基础上:
  Step 1: 利用双向 O-Gate 重新匹配 Sugar_Sequence (不重跑 Phase 2-6)
  Step 2: Rescue Phase 1 (策略 A/C/D) — 泛指糖精确修复
  Step 3: Rescue Phase 2 (策略 E/F/G) — 深度修复
  Step 4: PMC 抓取 + NLP 文本挖掘 + Safe Merge
  Step 5: 统计汇总 + v4 HTML 报告重新生成

On existing GlycoNP_Deep_Enriched_v2.csv:
  Step 1: Re-match Sugar_Sequence using bidirectional O-Gate
  Step 2: Rescue Phase 1 (strategies A/C/D)
  Step 3: Rescue Phase 2 (strategies E/F/G)
  Step 4: PMC fetch + NLP mining + Safe Merge
  Step 5: Statistics + v4 HTML report regeneration

Usage:
  python scripts/run_full_refresh.py [--skip-pmc] [--pmc-limit N]
"""
import argparse
import os
import re
import sys
import time
from collections import Counter

import pandas as pd
from rdkit import Chem
from tqdm import tqdm

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from lib.monosaccharide_identifier import generate_refined_sequence

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
REPORT_DIR = os.path.join(BASE_DIR, "reports")
INPUT_CSV = os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched_v2.csv")
OUTPUT_CSV = os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched_v3.csv")

GENERIC_PAT = re.compile(r'\b(Hex|Pen|dHex|HexA|Non|Oct|Hept)\b')


# =====================================================================
# Step 1: 糖序列重新匹配 (Sugar Sequence Re-match)
# 设计意图: 仅重新运行 sugar_sequence 的 identify 逻辑,
# 不触碰 Glycan_SMILES/Aglycon_SMILES (保留现有切分数据)
# Design: Only re-run sugar sequence identification logic,
# preserve existing Glycan_SMILES/Aglycon_SMILES splits
# =====================================================================

def step1_rematch(df: pd.DataFrame) -> pd.DataFrame:
    """利用双向 O-Gate 重新匹配所有行的 Sugar_Sequence。
    Re-match Sugar_Sequence for all rows using bidirectional O-Gate.
    """
    print("\n" + "=" * 70)
    print("  Step 1: Sugar Sequence Re-match (Bidirectional O-Gate)")
    print("  糖序列重新匹配 (双向氧门控)")
    print("=" * 70)
    t0 = time.time()

    # 保存旧序列用于对比 (Save old sequences for comparison)
    oldSeqCol = "Sugar_Sequence_Before_Refresh"
    df[oldSeqCol] = df["Sugar_Sequence"].copy()

    # 使用 Glycan_SMILES 重新生成 Sugar_Sequence
    # Use Glycan_SMILES to regenerate Sugar_Sequence
    glycanCol = "Glycan_SMILES"
    seqChanged = 0
    lcolBefore = df["Sugar_Sequence"].astype(str).str.contains(
        "L-Col", na=False).sum()

    for idx in tqdm(df.index, desc="  Re-matching", ncols=80):
        glycanSmi = str(df.at[idx, glycanCol]) if pd.notna(
            df.at[idx, glycanCol]) else ""
        if not glycanSmi or glycanSmi in ("nan", "None", "", "NULL"):
            continue
        try:
            glycanMol = Chem.MolFromSmiles(glycanSmi)
            if glycanMol is None:
                continue
            newSeq, newMods = generate_refined_sequence(glycanMol)
            if newSeq:
                oldSeq = str(df.at[idx, "Sugar_Sequence"])
                df.at[idx, "Sugar_Sequence"] = newSeq
                if newSeq != oldSeq:
                    seqChanged += 1
                # 同步更新修饰列 (Sync modification column)
                if newMods and "Glycan_Modifications" in df.columns:
                    df.at[idx, "Glycan_Modifications"] = newMods
        except Exception:
            pass

    lcolAfter = df["Sugar_Sequence"].astype(str).str.contains(
        "L-Col", na=False).sum()

    elapsed = time.time() - t0
    print(f"\n  Re-match Results ({elapsed:.0f}s):")
    print(f"    Sequences changed: {seqChanged:,}")
    print(f"    L-Col: {lcolBefore:,} → {lcolAfter:,} "
          f"(Δ = {lcolAfter - lcolBefore:+,})")

    return df


# =====================================================================
# Step 2: Rescue Phase 1 (A/C/D)
# =====================================================================

def step2_rescuePhase1(df: pd.DataFrame) -> pd.DataFrame:
    """执行 Rescue Phase 1: 文本挖掘(A) + 骨架规则(C) + 统计推断(D)。
    Execute Rescue Phase 1: text mining(A) + scaffold rules(C) + stats(D).
    """
    print("\n" + "=" * 70)
    print("  Step 2: Rescue Phase 1 (Strategies A/C/D)")
    print("=" * 70)

    # 导入 rescue 函数 (Import rescue functions)
    from scripts.rescue_generic_sugars import (
        rescueSequence, buildStatisticalPrior,
    )

    t0 = time.time()
    targetMask = df["Sugar_Sequence"].str.contains(
        GENERIC_PAT, na=False)
    targetCount = targetMask.sum()
    print(f"  Rows with generic sugars: {targetCount:,}")

    if targetCount == 0:
        print("  No generic sugars found, skipping.")
        return df

    scaffoldPrior, classPrior = buildStatisticalPrior(df)
    logCounter = {"A": 0, "C": 0, "D": 0, "MISS": 0}
    rescued = 0

    for idx in tqdm(df.index[targetMask], desc="  Rescue P1", ncols=80):
        oldSeq = str(df.at[idx, "Sugar_Sequence"])
        newSeq = rescueSequence(
            oldSeq,
            df.at[idx, "name"] if pd.notna(df.at[idx, "name"]) else None,
            df.at[idx, "iupac_name"] if "iupac_name" in df.columns and pd.notna(df.at[idx, "iupac_name"]) else None,
            df.at[idx, "synonyms"] if "synonyms" in df.columns and pd.notna(df.at[idx, "synonyms"]) else None,
            df.at[idx, "Superclass"] if "Superclass" in df.columns and pd.notna(df.at[idx, "Superclass"]) else None,
            df.at[idx, "Murcko_Scaffold"] if "Murcko_Scaffold" in df.columns and pd.notna(df.at[idx, "Murcko_Scaffold"]) else None,
            scaffoldPrior, classPrior, logCounter,
        )
        if newSeq != oldSeq:
            df.at[idx, "Sugar_Sequence"] = newSeq
            rescued += 1

    print(f"\n  Phase 1 Results ({time.time()-t0:.0f}s):")
    print(f"    Modified: {rescued:,} / {targetCount:,}")
    print(f"    A(Text):  {logCounter['A']:>6,}")
    print(f"    C(Rule):  {logCounter['C']:>6,}")
    print(f"    D(Stats): {logCounter['D']:>6,}")
    print(f"    MISS:     {logCounter['MISS']:>6,}")

    return df


# =====================================================================
# Step 3: Rescue Phase 2 (E/F/G)
# =====================================================================

def step3_rescuePhase2(df: pd.DataFrame) -> pd.DataFrame:
    """执行 Rescue Phase 2: IUPAC手性(E) + 生物合成(F) + InChIKey(G)。
    Execute Rescue Phase 2: IUPAC chiral(E) + biological(F) + InChIKey(G).
    """
    print("\n" + "=" * 70)
    print("  Step 3: Rescue Phase 2 (Strategies E/F/G)")
    print("=" * 70)

    from scripts.rescue_generic_sugars_phase2 import (
        rescueSequencePhase2, buildInchiKeyPrior, GENERIC_PATTERN,
    )

    t0 = time.time()
    targetMask = df["Sugar_Sequence"].str.contains(
        GENERIC_PAT, na=False)
    targetCount = targetMask.sum()
    print(f"  Remaining generic: {targetCount:,}")

    if targetCount == 0:
        print("  All resolved, skipping.")
        return df

    inchiKeyPrior = buildInchiKeyPrior(df)
    logCounter = {"E": 0, "F": 0, "G": 0, "MISS": 0}
    rescued = 0

    for idx in tqdm(df.index[targetMask], desc="  Rescue P2", ncols=80):
        oldSeq = str(df.at[idx, "Sugar_Sequence"])
        newSeq = rescueSequencePhase2(
            oldSeq,
            str(df.at[idx, "iupac_name"]) if "iupac_name" in df.columns and pd.notna(df.at[idx, "iupac_name"]) else None,
            str(df.at[idx, "Superclass"]) if "Superclass" in df.columns and pd.notna(df.at[idx, "Superclass"]) else None,
            str(df.at[idx, "canonical_smiles"]) if pd.notna(df.at[idx, "canonical_smiles"]) else None,
            str(df.at[idx, "standard_inchi_key"]) if "standard_inchi_key" in df.columns and pd.notna(df.at[idx, "standard_inchi_key"]) else None,
            inchiKeyPrior, logCounter,
        )
        if newSeq != oldSeq:
            df.at[idx, "Sugar_Sequence"] = newSeq
            rescued += 1

    print(f"\n  Phase 2 Results ({time.time()-t0:.0f}s):")
    print(f"    Modified: {rescued:,} / {targetCount:,}")
    print(f"    E(IUPAC):   {logCounter['E']:>6,}")
    print(f"    F(Bio):     {logCounter['F']:>6,}")
    print(f"    G(InChI):   {logCounter['G']:>6,}")
    print(f"    MISS:       {logCounter['MISS']:>6,}")

    return df


# =====================================================================
# Step 4: PMC 抓取 + NLP 挖掘 + Safe Merge
# =====================================================================

def step4_nlpMining(df: pd.DataFrame, skipPmc: bool = False,
                    pmcLimit: int = 0) -> pd.DataFrame:
    """PMC 增量抓取 + NLP 糖名/活性挖掘 + 零污染安全合并。
    PMC incremental fetch + NLP sugar/bioactivity mining + safe merge.
    """
    print("\n" + "=" * 70)
    print("  Step 4: PMC Fetch + NLP Mining + Safe Merge")
    print("=" * 70)

    from scripts.nlp_literature_mining import (
        runPmcFetch, runNlpMining, runSafeMerge,
    )

    # 选择目标: 有 DOI 且序列含泛指糖的行
    # Select targets: rows with DOI and generic sugar in sequence
    hasDoi = (df["dois"].notna() &
              ~df["dois"].astype(str).isin(["", "nan", "None"]))
    hasGeneric = df["Sugar_Sequence"].str.contains(
        GENERIC_PAT, na=False)

    # 全部有 DOI 的行都作为 NLP 挖掘目标
    # All rows with DOI are NLP mining targets
    targetIdxs = df.index[hasDoi].tolist()

    if pmcLimit > 0:
        targetIdxs = targetIdxs[:pmcLimit]

    print(f"  NLP targets: {len(targetIdxs):,} (with DOI)")
    print(f"  Still generic: {hasGeneric.sum():,}")

    # PMC 抓取 (PMC Fetch)
    if not skipPmc:
        df = runPmcFetch(df, targetIdxs)
    else:
        print("  [SKIP] PMC fetch (--skip-pmc)")

    # NLP 挖掘 (NLP Mining)
    df = runNlpMining(df, targetIdxs)

    # 安全合并 (Safe Merge)
    df = runSafeMerge(df)

    return df


# =====================================================================
# Step 5: 统计汇总 (Statistics Summary)
# =====================================================================

def step5_statistics(df: pd.DataFrame):
    """打印修复前后的全量统计。
    Print before/after statistics for the full refresh.
    """
    print("\n" + "=" * 70)
    print("  Step 5: Refresh Statistics Summary")
    print("  刷新统计汇总")
    print("=" * 70)

    oldCol = "Sugar_Sequence_Before_Refresh"
    seqCol = "Sugar_Sequence"
    consensusCol = "Consensus_Sugar_Sequence"

    # L-Col 统计 (L-Col stats)
    if oldCol in df.columns:
        lcolBefore = df[oldCol].astype(str).str.contains(
            "L-Col", na=False).sum()
    else:
        lcolBefore = "N/A"
    lcolAfter = df[seqCol].astype(str).str.contains(
        "L-Col", na=False).sum()
    lcolConsensus = (df[consensusCol].astype(str).str.contains(
        "L-Col", na=False).sum() if consensusCol in df.columns else "N/A")

    print(f"\n  L-Col Records:")
    print(f"    Before refresh: {lcolBefore}")
    print(f"    After re-match: {lcolAfter}")
    print(f"    In Consensus:   {lcolConsensus}")

    # 泛指糖统计 (Generic sugar stats)
    genericPat = r'\bHex\b|\bPen\b|\bdHex\b|\bHexA\b|\bNon\b|\bOct\b|\bHept\b'
    if oldCol in df.columns:
        genericBefore = df[oldCol].str.contains(
            genericPat, na=False, regex=True).sum()
    else:
        genericBefore = "N/A"
    genericAfter = df[seqCol].str.contains(
        genericPat, na=False, regex=True).sum()
    genericConsensus = (df[consensusCol].str.contains(
        genericPat, na=False, regex=True).sum()
        if consensusCol in df.columns else "N/A")

    print(f"\n  Generic Sugar Rows:")
    print(f"    Before: {genericBefore}")
    print(f"    After:  {genericAfter}")
    print(f"    Consensus: {genericConsensus}")

    # 糖名分布 Top 15 (Sugar name distribution top 15)
    sugarPattern = (
        r'Neu5Ac|Neu5Gc|KDO|'
        r'[DL]-[A-Z][a-z]+[A-Z]?[a-z]*(?:\(NLP\))?|'
        r'Hex|dHex|Pen|HexA|HexN|HexNAc|Non|Oct|Hept'
    )
    targetCol = consensusCol if consensusCol in df.columns else seqCol
    allTokens = []
    for seq in df[targetCol].dropna():
        allTokens.extend(re.findall(sugarPattern, str(seq)))
    tokenCounter = Counter(allTokens)

    print(f"\n  Sugar Token Distribution (Top 20):")
    for token, count in tokenCounter.most_common(20):
        print(f"    {token:<25s} {count:>7,}")

    print(f"\n  Total unique tokens: {len(tokenCounter):,}")
    print(f"  Total token count:   {sum(tokenCounter.values()):,}")


# =====================================================================
# Main
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        description="GlycoNP Full Sugar Sequence Refresh")
    parser.add_argument("--skip-pmc", action="store_true",
                        help="跳过 PMC 抓取 / Skip PMC fetch")
    parser.add_argument("--pmc-limit", type=int, default=0,
                        help="PMC 抓取行数上限 / Max rows for PMC fetch (0=all)")
    args = parser.parse_args()

    print("=" * 70)
    print("  GlycoNP Full Sugar Sequence Refresh Pipeline")
    print("  全量糖序列刷新管线 (双向 O-Gate + 10 策略)")
    print("=" * 70)
    t0 = time.time()

    # 加载数据 (Load data)
    print(f"\n  Loading: {INPUT_CSV}")
    df = pd.read_csv(INPUT_CSV, low_memory=False, dtype=str,
                     encoding="utf-8-sig")
    print(f"  Rows: {len(df):,}, Columns: {len(df.columns)}")

    # Step 1: 糖序列重新匹配 (Sugar sequence re-match)
    df = step1_rematch(df)

    # Step 2: Rescue Phase 1 (A/C/D)
    df = step2_rescuePhase1(df)

    # Step 3: Rescue Phase 2 (E/F/G)
    df = step3_rescuePhase2(df)

    # Step 4: PMC + NLP + Safe Merge
    df = step4_nlpMining(df, skipPmc=args.skip_pmc,
                         pmcLimit=args.pmc_limit)

    # Step 5: 统计汇总 (Statistics)
    step5_statistics(df)

    # 清理临时列 (Clean temp columns)
    if "Sugar_Sequence_Before_Refresh" in df.columns:
        df.drop(columns=["Sugar_Sequence_Before_Refresh"], inplace=True)

    # 保存 (Save)
    print(f"\n  Saving to: {OUTPUT_CSV}")
    df.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")
    print(f"  Saved: {len(df):,} rows, {len(df.columns)} columns")

    totalTime = time.time() - t0
    print(f"\n  Total time: {totalTime:.0f}s ({totalTime/60:.1f}min)")
    print("=" * 70)
    print("  [OK] Full refresh complete!")
    print("  Next: python scripts/dict_recalc_v3.py")
    print("    (reads GlycoNP_Deep_Enriched_v3.csv, generates v4 HTML)")
    print("=" * 70)


if __name__ == "__main__":
    main()

