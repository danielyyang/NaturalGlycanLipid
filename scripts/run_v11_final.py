"""
V11 终极全库重刷 (V11 Full Pipeline Final Rerun)
================================================

用 KEGG 完美字典 + v10 工业级管线对全库 94,242 条数据进行终极重刷。

Pipeline:
  Step 1: v10 引擎全量 Sugar_Sequence 重匹配
          (虚拟水解剥离修饰 -> 裸糖匹配 -> KEGG 字典)
  Step 2: Rescue Phase 1 (A/C/D) — 泛指糖精确修复
  Step 3: Rescue Phase 2 (E/F/G) — 深度修复
  Step 4: 终极糖名分布榜单 + 保存 FINAL.csv

V11 Full Pipeline Final Rerun with KEGG-validated dictionary.
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
INPUT_CSV = os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched_v8.csv")
OUTPUT_CSV = os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched_V12.csv")

GENERIC_PAT = re.compile(r'\b(Hex|Pen|dHex|HexA|Non|Oct|Hept)\b')

# =====================================================================
# 任务 1 & 2: v10 全量重匹配 (v10 Full Re-match)
# =====================================================================
def step1_v10Rematch(df: pd.DataFrame) -> pd.DataFrame:
    """用 v10 管线对全量数据重匹配 Sugar_Sequence。
    Re-match Sugar_Sequence using v10 pipeline with KEGG dictionary."""
    print("\n" + "=" * 70)
    print("  Step 1: V10 Industrial Pipeline Full Re-match")
    print("  V10 工业级管线全量重匹配 (KEGG 完美字典)")
    print("=" * 70)
    t0 = time.time()

    # 保存旧序列 (Save old sequences for comparison)
    oldSeqCol = "Sugar_Sequence_v8"
    df[oldSeqCol] = df["Sugar_Sequence"].copy()

    glycanCol = "Glycan_SMILES"
    seqChanged = 0
    errorCount = 0
    processed = 0

    totalRows = len(df)
    for idx in tqdm(df.index, desc="  V10 Re-match", ncols=80):
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
                if newMods and "Glycan_Modifications" in df.columns:
                    df.at[idx, "Glycan_Modifications"] = newMods
            processed += 1
        except Exception:
            errorCount += 1

    elapsed = time.time() - t0
    print(f"\n  V10 Re-match Results ({elapsed:.0f}s):")
    print(f"    Processed:  {processed:,} / {totalRows:,}")
    print(f"    Changed:    {seqChanged:,}")
    print(f"    Errors:     {errorCount:,}")

    return df


# =====================================================================
# 任务 3: Rescue Phase 1 (A/C/D)
# =====================================================================
def step2_rescuePhase1(df: pd.DataFrame) -> pd.DataFrame:
    """Rescue Phase 1: 文本挖掘(A) + 骨架规则(C) + 统计推断(D)。"""
    print("\n" + "=" * 70)
    print("  Step 2: Rescue Phase 1 (Strategies A/C/D)")
    print("=" * 70)

    try:
        from scripts.rescue_generic_sugars import (
            rescueSequence, buildStatisticalPrior,
        )
    except ImportError:
        print("  [SKIP] rescue_generic_sugars not found")
        return df

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

    nameCol = "name" if "name" in df.columns else "COMPOUND_NAME"

    for idx in tqdm(df.index[targetMask], desc="  Rescue P1", ncols=80):
        oldSeq = str(df.at[idx, "Sugar_Sequence"])
        try:
            newSeq = rescueSequence(
                oldSeq,
                df.at[idx, nameCol] if nameCol in df.columns and pd.notna(df.at[idx, nameCol]) else None,
                df.at[idx, "iupac_name"] if "iupac_name" in df.columns and pd.notna(df.at[idx, "iupac_name"]) else None,
                df.at[idx, "synonyms"] if "synonyms" in df.columns and pd.notna(df.at[idx, "synonyms"]) else None,
                df.at[idx, "Superclass"] if "Superclass" in df.columns and pd.notna(df.at[idx, "Superclass"]) else None,
                df.at[idx, "Murcko_Scaffold"] if "Murcko_Scaffold" in df.columns and pd.notna(df.at[idx, "Murcko_Scaffold"]) else None,
                scaffoldPrior, classPrior, logCounter,
            )
            if newSeq != oldSeq:
                df.at[idx, "Sugar_Sequence"] = newSeq
                rescued += 1
        except Exception:
            pass

    print(f"\n  Phase 1 Results ({time.time()-t0:.0f}s):")
    print(f"    Modified: {rescued:,} / {targetCount:,}")
    print(f"    A(Text):  {logCounter['A']:>6,}")
    print(f"    C(Rule):  {logCounter['C']:>6,}")
    print(f"    D(Stats): {logCounter['D']:>6,}")
    print(f"    MISS:     {logCounter['MISS']:>6,}")

    return df


# =====================================================================
# Rescue Phase 2 (E/F/G)
# =====================================================================
def step3_rescuePhase2(df: pd.DataFrame) -> pd.DataFrame:
    """Rescue Phase 2: IUPAC(E) + Biological(F) + InChIKey(G)."""
    print("\n" + "=" * 70)
    print("  Step 3: Rescue Phase 2 (Strategies E/F/G)")
    print("=" * 70)

    try:
        from scripts.rescue_generic_sugars_phase2 import (
            rescueSequencePhase2, buildInchiKeyPrior,
        )
    except ImportError:
        print("  [SKIP] rescue_generic_sugars_phase2 not found")
        return df

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
        try:
            newSeq = rescueSequencePhase2(
                oldSeq,
                str(df.at[idx, "iupac_name"]) if "iupac_name" in df.columns and pd.notna(df.at[idx, "iupac_name"]) else None,
                str(df.at[idx, "Superclass"]) if "Superclass" in df.columns and pd.notna(df.at[idx, "Superclass"]) else None,
                str(df.at[idx, "canonical_smiles"]) if "canonical_smiles" in df.columns and pd.notna(df.at[idx, "canonical_smiles"]) else None,
                str(df.at[idx, "standard_inchi_key"]) if "standard_inchi_key" in df.columns and pd.notna(df.at[idx, "standard_inchi_key"]) else None,
                inchiKeyPrior, logCounter,
            )
            if newSeq != oldSeq:
                df.at[idx, "Sugar_Sequence"] = newSeq
                rescued += 1
        except Exception:
            pass

    print(f"\n  Phase 2 Results ({time.time()-t0:.0f}s):")
    print(f"    Modified: {rescued:,} / {targetCount:,}")
    print(f"    E(IUPAC):   {logCounter['E']:>6,}")
    print(f"    F(Bio):     {logCounter['F']:>6,}")
    print(f"    G(InChI):   {logCounter['G']:>6,}")
    print(f"    MISS:       {logCounter['MISS']:>6,}")

    return df


# =====================================================================
# Step 4: V11 终极糖名分布榜单
# =====================================================================
def step4_finalReport(df: pd.DataFrame):
    """V11 终极糖名分布榜单。"""
    print("\n" + "=" * 70)
    print("  Step 4: V11 Final Sugar Distribution Report")
    print("  V11 终极糖名分布榜单")
    print("=" * 70)

    seqCol = "Sugar_Sequence"
    oldCol = "Sugar_Sequence_v8"

    # v8 分布 (Before)
    sugarPattern = (
        r'Neu5Ac|Neu5Gc|KDO|'
        r'[DL]-[A-Z][a-z]+[A-Z]?[a-z]*(?:\([^\)]*\))?|'
        r'Hex|dHex|Pen|HexA|HexN|HexNAc|Non|Oct|Hept'
    )

    v8Tokens = Counter()
    for seq in df[oldCol].dropna():
        v8Tokens.update(re.findall(sugarPattern, str(seq)))

    v11Tokens = Counter()
    for seq in df[seqCol].dropna():
        v11Tokens.update(re.findall(sugarPattern, str(seq)))

    print(f"\n  {'Rank':>4s}  {'Sugar':^30s}  {'v8':>8s}  {'V11':>8s}  {'Delta':>8s}")
    print("  " + "-" * 65)

    # 合并所有出现过的糖名
    allSugars = set(list(v8Tokens.keys()) + list(v11Tokens.keys()))
    v11Sorted = sorted(allSugars, key=lambda x: -v11Tokens.get(x, 0))

    for rank, sugar in enumerate(v11Sorted[:30], 1):
        v8c = v8Tokens.get(sugar, 0)
        v11c = v11Tokens.get(sugar, 0)
        delta = v11c - v8c
        marker = ""
        if delta > 100:
            marker = " ++"
        elif delta < -100:
            marker = " --"
        print(f"  {rank:4d}  {sugar:<30s}  {v8c:8,}  {v11c:8,}  {delta:+8,}{marker}")

    print(f"\n  Total unique sugar names: v8={len(v8Tokens):,}, V11={len(v11Tokens):,}")
    print(f"  Total sugar tokens: v8={sum(v8Tokens.values()):,}, V11={sum(v11Tokens.values()):,}")

    # 关键指标
    hexV8 = v8Tokens.get("Hex", 0)
    hexV11 = v11Tokens.get("Hex", 0)
    talV8 = v8Tokens.get("D-Tal", 0)
    talV11 = v11Tokens.get("D-Tal", 0)
    glcV11 = v11Tokens.get("D-Glc", 0)
    galV11 = v11Tokens.get("D-Gal", 0)
    rhaV11 = v11Tokens.get("L-Rha", 0)

    print(f"\n  KEY METRICS:")
    print(f"    Hex:   {hexV8:,} -> {hexV11:,} ({hexV11-hexV8:+,})")
    print(f"    D-Tal: {talV8:,} -> {talV11:,} ({talV11-talV8:+,})")
    print(f"    D-Glc: {v8Tokens.get('D-Glc',0):,} -> {glcV11:,}")
    print(f"    D-Gal: {v8Tokens.get('D-Gal',0):,} -> {galV11:,}")
    print(f"    L-Rha: {v8Tokens.get('L-Rha',0):,} -> {rhaV11:,}")


# =====================================================================
# Main
# =====================================================================
def main():
    parser = argparse.ArgumentParser(
        description="V11 Full Pipeline Final Rerun")
    parser.add_argument("--skip-rescue", action="store_true",
                        help="Skip rescue phases 1 & 2")
    args = parser.parse_args()

    print("=" * 70)
    print("  V11 FINAL PIPELINE — KEGG Perfect Dictionary")
    print("  V11 终极管线 — KEGG 完美字典全库重刷")
    print("=" * 70)
    t0 = time.time()

    # Dictionary spot check
    from lib.glycan_reference_library import RAW_MONOSACCHARIDE_SMILES
    print(f"\n  Dictionary entries: {len(RAW_MONOSACCHARIDE_SMILES)}")
    print(f"  D-Tal: {RAW_MONOSACCHARIDE_SMILES.get(('D-Tal','a'),'MISSING')}")
    print(f"  D-Gal: {RAW_MONOSACCHARIDE_SMILES.get(('D-Gal','a'),'MISSING')}")

    # Load data
    print(f"\n  Loading: {INPUT_CSV}")
    df = pd.read_csv(INPUT_CSV, low_memory=False)
    print(f"  Rows: {len(df):,}, Columns: {len(df.columns)}")

    # Step 1: V10 Re-match
    df = step1_v10Rematch(df)

    # Step 2-3: Rescue (optional)
    if not args.skip_rescue:
        df = step2_rescuePhase1(df)
        df = step3_rescuePhase2(df)
    else:
        print("\n  [SKIP] Rescue phases (--skip-rescue)")

    # Step 4: Final report
    step4_finalReport(df)

    # Clean temp columns
    if "Sugar_Sequence_v8" in df.columns:
        df.drop(columns=["Sugar_Sequence_v8"], inplace=True)

    # Save
    print(f"\n  Saving: {OUTPUT_CSV}")
    df.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")
    print(f"  Saved: {len(df):,} rows")

    totalTime = time.time() - t0
    print(f"\n  Total time: {totalTime:.0f}s ({totalTime/60:.1f}min)")
    print("=" * 70)
    print("  V11 FINAL PIPELINE COMPLETE!")
    print("=" * 70)


if __name__ == "__main__":
    main()
