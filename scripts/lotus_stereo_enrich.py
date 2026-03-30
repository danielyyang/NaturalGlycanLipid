"""
LOTUS 立体化学补全: COCONUT → LOTUS 手性信息比对与富集
LOTUS Stereochemistry Enrichment: Compare & enrich COCONUT SMILES with LOTUS chirality

策略 (Strategy):
  1. 用 InChIKey Block1 (前14位, 2D 骨架匹配) 关联两库
  2. 比较 COCONUT canonical_smiles 与 LOTUS structure_smiles 的手性信息
  3. 当 LOTUS 手性更完整时, 补全到新列 LOTUS_Enriched_SMILES
  4. 输出 Enriched CSV, 所有富集列均有 [Enriched] 标注

Design intent: COCONUT's canonical_smiles often lacks stereochemistry
(wedge/dash bonds) because many structures were auto-converted from 2D.
LOTUS has manually curated or better-converted 3D structures.
"""
import os
import sys
import time

import pandas as pd
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from rdkit import Chem

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
REPORT_DIR = os.path.join(BASE_DIR, "reports")
INPUT_CSV = os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched_FINAL.csv")
LOTUS_GZ = os.path.join(BASE_DIR, "data", "230106_frozen_metadata.csv.gz")
OUTPUT_CSV = os.path.join(REPORT_DIR, "GlycoNP_LOTUS_Enriched.csv")


def countStereocenters(smiles: str) -> int:
    """统计 SMILES 中指定手性的原子数 (@ 和 @@ 标记)。
    Count number of specified stereocenters in a SMILES string."""
    if not smiles or smiles == "nan":
        return 0
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    count = 0
    for atom in mol.GetAtoms():
        if atom.HasProp("_CIPCode"):
            count += 1
    return count


def main():
    print("=" * 70)
    print("  LOTUS Stereochemistry Enrichment Pipeline")
    print("  LOTUS 立体化学补全管线")
    print("=" * 70)
    t0 = time.time()

    # ================================================================
    # Step 1: 加载 GlycoNP FINAL (Load GlycoNP FINAL)
    # ================================================================
    print(f"\n  [1/6] Loading GlycoNP FINAL: {INPUT_CSV}")
    df = pd.read_csv(INPUT_CSV, low_memory=False)
    total = len(df)
    print(f"    Rows: {total:,}")

    # 确保有 InChIKey block1 列 (Ensure InChIKey block1 column)
    if "non_isomeric_inchikey_block1" not in df.columns:
        # 先从 standard_inchi_key 截取
        if "standard_inchi_key" in df.columns:
            df["non_isomeric_inchikey_block1"] = (
                df["standard_inchi_key"].astype(str).str[:14]
            )
        else:
            print("  [ERROR] No InChIKey column available!")
            return

    # ================================================================
    # Step 2: 加载 LOTUS (Load LOTUS)
    # ================================================================
    print(f"\n  [2/6] Loading LOTUS: {LOTUS_GZ}")
    lotusCols = [
        "structure_inchikey",
        "structure_smiles",        # 带手性的 3D SMILES
        "structure_smiles_2D",     # 平面 2D SMILES
        "structure_stereocenters_total",
        "structure_stereocenters_unspecified",
    ]
    lotusDf = pd.read_csv(LOTUS_GZ, dtype=str, low_memory=False,
                          usecols=lotusCols)
    print(f"    LOTUS rows: {len(lotusDf):,}")

    # 生成 LOTUS Block1 (Generate LOTUS Block1)
    lotusDf["_block1"] = lotusDf["structure_inchikey"].astype(str).str[:14]
    lotusDf = lotusDf.dropna(subset=["_block1"])
    lotusDf = lotusDf[lotusDf["_block1"] != "nan"]

    # 去重: 同一 Block1 取 unspecified 最少的 (优先选手性完整的)
    # Dedup: for same Block1, pick the one with fewest unspecified stereocenters
    lotusDf["_unspec"] = pd.to_numeric(
        lotusDf["structure_stereocenters_unspecified"], errors="coerce"
    ).fillna(999)
    lotusDf["_total"] = pd.to_numeric(
        lotusDf["structure_stereocenters_total"], errors="coerce"
    ).fillna(0)
    lotusDf = lotusDf.sort_values("_unspec").drop_duplicates(
        subset=["_block1"], keep="first"
    )
    print(f"    Unique LOTUS Block1 entries: {len(lotusDf):,}")

    # ================================================================
    # Step 3: 匹配 (Merge on Block1)
    # ================================================================
    print(f"\n  [3/6] Matching by InChIKey Block1...")
    merged = df.merge(
        lotusDf[["_block1", "structure_smiles", "structure_smiles_2D",
                 "_total", "_unspec"]],
        left_on="non_isomeric_inchikey_block1",
        right_on="_block1",
        how="left",
    )
    hitCount = merged["_block1"].notna().sum()
    print(f"    Block1 hits: {hitCount:,} / {total:,} ({hitCount/total*100:.1f}%)")

    # ================================================================
    # Step 4: 立体化学比对 (Stereochemistry comparison)
    # ================================================================
    print(f"\n  [4/6] Comparing stereochemistry (COCONUT vs LOTUS)...")

    # 初始化 Enriched 列 (Initialize enriched columns)
    merged["[Enriched]_LOTUS_SMILES_3D"] = np.nan
    merged["[Enriched]_Stereo_Source"] = np.nan
    merged["[Enriched]_COCONUT_Stereocenters"] = np.nan
    merged["[Enriched]_LOTUS_Stereocenters_Total"] = np.nan
    merged["[Enriched]_LOTUS_Stereocenters_Unspecified"] = np.nan
    merged["[Enriched]_Stereo_Gain"] = np.nan

    enriched = 0
    noGain = 0
    noLotus = 0
    lotusWorse = 0
    processingErrors = 0

    # 只处理有 LOTUS 匹配的行
    hasMask = merged["_block1"].notna()
    hasIdxs = merged.index[hasMask]
    batchSize = 5000
    totalToProcess = len(hasIdxs)

    for batchStart in range(0, totalToProcess, batchSize):
        batchEnd = min(batchStart + batchSize, totalToProcess)
        batchIdxs = hasIdxs[batchStart:batchEnd]

        for idx in batchIdxs:
            coconutSmi = str(merged.at[idx, "canonical_smiles"]) \
                if pd.notna(merged.at[idx, "canonical_smiles"]) else ""
            lotusSmi = str(merged.at[idx, "structure_smiles"]) \
                if pd.notna(merged.at[idx, "structure_smiles"]) else ""
            lotusTotal = merged.at[idx, "_total"]
            lotusUnspec = merged.at[idx, "_unspec"]

            if not lotusSmi or lotusSmi == "nan":
                noLotus += 1
                continue

            try:
                coconutStereo = countStereocenters(coconutSmi)
                lotusStereo = int(lotusTotal) - int(lotusUnspec) \
                    if lotusTotal is not None and lotusUnspec is not None \
                    else countStereocenters(lotusSmi)

                stereoGain = lotusStereo - coconutStereo

                # 记录比较数据 (Record comparison data)
                merged.at[idx, "[Enriched]_COCONUT_Stereocenters"] = coconutStereo
                merged.at[idx, "[Enriched]_LOTUS_Stereocenters_Total"] = int(lotusTotal) if pd.notna(lotusTotal) else None
                merged.at[idx, "[Enriched]_LOTUS_Stereocenters_Unspecified"] = int(lotusUnspec) if pd.notna(lotusUnspec) and lotusUnspec != 999 else None
                merged.at[idx, "[Enriched]_Stereo_Gain"] = stereoGain

                if stereoGain > 0:
                    # LOTUS 手性更完整, 补全 (LOTUS has better chirality)
                    merged.at[idx, "[Enriched]_LOTUS_SMILES_3D"] = lotusSmi
                    merged.at[idx, "[Enriched]_Stereo_Source"] = "LOTUS_BETTER"
                    enriched += 1
                elif stereoGain == 0:
                    merged.at[idx, "[Enriched]_Stereo_Source"] = "EQUAL"
                    noGain += 1
                else:
                    merged.at[idx, "[Enriched]_Stereo_Source"] = "COCONUT_BETTER"
                    lotusWorse += 1

            except Exception:
                processingErrors += 1

        elapsed = time.time() - t0
        print(f"    Processed {batchEnd:,}/{totalToProcess:,} "
              f"({elapsed:.0f}s) | Enriched: {enriched:,}")

    # ================================================================
    # Step 5: 统计报告 (Statistics Report)
    # ================================================================
    print(f"\n{'='*70}")
    print(f"  STEREO ENRICHMENT REPORT")
    print(f"  立体化学补全报告")
    print(f"{'='*70}")
    print(f"  Total rows:             {total:,}")
    print(f"  LOTUS matched:          {hitCount:,} ({hitCount/total*100:.1f}%)")
    print(f"  ---")
    print(f"  LOTUS better (enriched):{enriched:,}")
    print(f"  Equal stereo:           {noGain:,}")
    print(f"  COCONUT better:         {lotusWorse:,}")
    print(f"  No LOTUS SMILES:        {noLotus:,}")
    print(f"  Processing errors:      {processingErrors:,}")

    # 手性提升分布 (Stereo gain distribution)
    gainCol = merged["[Enriched]_Stereo_Gain"].dropna()
    if len(gainCol) > 0:
        print(f"\n  Stereo Gain Distribution:")
        gainCounts = gainCol.value_counts().sort_index()
        for gain, count in gainCounts.items():
            if count >= 100:
                label = "LOTUS better" if gain > 0 else ("equal" if gain == 0 else "COCONUT better")
                print(f"    Gain={int(gain):+3d}: {count:>7,}  ({label})")

    # ================================================================
    # Step 6: 保存 (Save)
    # ================================================================
    # 清理临时列 (Clean temp columns)
    dropCols = [c for c in merged.columns if c.startswith("_")]
    merged.drop(columns=dropCols, inplace=True, errors="ignore")

    print(f"\n  [6/6] Saving: {OUTPUT_CSV}")
    merged.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")
    print(f"    Saved: {len(merged):,} rows, {len(merged.columns)} columns")

    # 列出所有 Enriched 列 (List all enriched columns)
    enrichedCols = [c for c in merged.columns if "[Enriched]" in c]
    print(f"\n  Enriched columns ({len(enrichedCols)}):")
    for col in enrichedCols:
        validCount = merged[col].notna().sum()
        print(f"    {col}: {validCount:,} values")

    totalTime = time.time() - t0
    print(f"\n  Total time: {totalTime:.0f}s ({totalTime/60:.1f}min)")
    print("=" * 70)
    print("  LOTUS Stereo Enrichment COMPLETE!")
    print("=" * 70)


if __name__ == "__main__":
    main()
