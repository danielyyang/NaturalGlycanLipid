"""
Deep Architecture Upgrade — Phases 2 + 3 + 5
=============================================

Phase 2: Scaffold Abstraction + Fine-Grained Classification
Phase 3: Sugar Chain Topology Analytics + Charts
Phase 5: ML-Ready Dataset Export

Usage:
  python scripts/deep_upgrade_235.py
"""
import argparse
import os
import re
import sys
from collections import Counter
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
REPORT_DIR = os.path.join(BASE_DIR, "reports")

# =====================================================================
# Phase 2: Scaffold Abstraction
# =====================================================================

# 设计意图: 将 np_classifier_superclass 映射到 10 个易于理解的"超级家族"。
# 不需要 SMARTS — COCONUT 的 NPClassifier 已经提供了高质量分类。

SUPERCLASS_MAP = {
    # === 萜类 (Terpenoids) ===
    "Triterpenoids":                "Triterpenoid saponins",
    "Sesquiterpenoids":             "Sesqui- & Sesterterpenoids",
    "Sesterterpenoids":             "Sesqui- & Sesterterpenoids",
    "Diterpenoids":                 "Diterpenoid glycosides",
    "Monoterpenoids":               "Monoterpenoid glycosides",

    # === 甾体 (Steroids) ===
    "Steroids":                     "Steroidal saponins",

    # === 黄酮 (Flavonoids) ===
    "Flavonoids":                   "Flavonoid glycosides",

    # === 聚酮 (Polyketides) ===
    "Macrolides":                   "Macrolide glycosides",
    "Polycyclic aromatic polyketides": "Aromatic polyketide glycosides",

    # === 苯丙素 (Phenylpropanoids) ===
    "Lignans":                      "Lignan glycosides",
    "Coumarins":                    "Coumarin glycosides",
    "Phenylpropanoids (C6-C3)":     "Phenylpropanoid glycosides",

    # === 糖类 (Carbohydrates) ===
    "Saccharides":                  "Oligosaccharides",

    # === 生物碱 (Alkaloids) ===
    "Indole alkaloids":             "Alkaloid glycosides",
    "Purine alkaloids":             "Alkaloid glycosides",
    "Steroidal alkaloids":          "Alkaloid glycosides",
    "Isoquinoline alkaloids":       "Alkaloid glycosides",
}

DEFAULT_SUPERCLASS = "Other glycosides"


def runPhase2(df: pd.DataFrame) -> pd.DataFrame:
    """
    Phase 2: 骨架降维与精细分类。
    Scaffold Abstraction & Fine-Grained Classification.
    """
    print("\n" + "=" * 70)
    print("  Phase 2: Scaffold Abstraction & Classification")
    print("  骨架降维与精细分类")
    print("=" * 70)

    # Detailed_NP_Class: 直接使用 np_classifier_class
    df["Detailed_NP_Class"] = df["np_classifier_class"].copy()

    # Super_Scaffold_Class: 映射 np_classifier_superclass
    def mapSuperclass(val: str) -> str:
        s = str(val).strip()
        if s in ("", "nan", "None"):
            return ""
        return SUPERCLASS_MAP.get(s, DEFAULT_SUPERCLASS)

    df["Super_Scaffold_Class"] = df["np_classifier_superclass"].apply(
        mapSuperclass)

    # 统计 (Stats)
    detailedValid = (df["Detailed_NP_Class"].notna() &
                     ~df["Detailed_NP_Class"].astype(str).isin(
                         ["", "nan", "None"])).sum()
    superValid = (df["Super_Scaffold_Class"] != "").sum()

    print(f"\n  Detailed_NP_Class coverage:      {detailedValid:,} / "
          f"{len(df):,} ({detailedValid/len(df)*100:.1f}%)")
    print(f"  Super_Scaffold_Class coverage:   {superValid:,} / "
          f"{len(df):,} ({superValid/len(df)*100:.1f}%)")

    # Top distribution
    print(f"\n  Super_Scaffold_Class distribution:")
    dist = df["Super_Scaffold_Class"].value_counts()
    for cls, count in dist.head(12).items():
        if cls:
            print(f"    {cls:35s}: {count:>6,}")

    print(f"\n  Detailed_NP_Class top 10:")
    dist2 = df["Detailed_NP_Class"].value_counts()
    for cls, count in dist2.head(10).items():
        print(f"    {cls:40s}: {count:>6,}")

    return df


# =====================================================================
# Phase 3: Sugar Chain Topology Analytics
# =====================================================================

SUGAR_TOKEN_RE = re.compile(
    r'Neu5Ac|Neu5Gc|KDO|'
    r'[DL]-[A-Za-z]+(?:\([^)]*\))?|'
    r'(?:Hex|Pen|dHex|HexA|HexN|Non|Oct|Hept)(?:[A-Z])?(?:\([^)]*\))?'
)


def parseSugarTopology(seq: str) -> Dict:
    """
    解析糖链拓扑特征。
    Parse sugar chain topology from a Sugar_Sequence string.

    规则:
    - '; ' 分隔不同连接位点的糖链 (independent sub-chains)
    - 每条子链中的糖 token 数 = 该链长度

    Returns:
      {
        "total": 总单糖数,
        "chains": [每条子链的长度],
        "max": 最长链的长度,
      }
    """
    if not seq or str(seq) in ("nan", "", "None"):
        return {"total": 0, "chains": [], "max": 0}

    # 按 ' ; ' 分割为子链
    subChains = [s.strip() for s in str(seq).split(" ; ") if s.strip()]

    chainLengths = []
    totalCount = 0

    for chain in subChains:
        tokens = SUGAR_TOKEN_RE.findall(chain)
        chainLen = len(tokens)
        chainLengths.append(chainLen)
        totalCount += chainLen

    return {
        "total": totalCount,
        "chains": chainLengths,
        "max": max(chainLengths) if chainLengths else 0,
    }


def runPhase3(df: pd.DataFrame) -> pd.DataFrame:
    """
    Phase 3: 糖链拓扑统计与图表渲染。
    Sugar Chain Topology Analytics & Chart Rendering.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns

    print("\n" + "=" * 70)
    print("  Phase 3: Sugar Chain Topology Analytics")
    print("  糖链拓扑解析与双轨统计")
    print("=" * 70)

    seqCol = "Consensus_Sugar_Sequence"
    if seqCol not in df.columns:
        seqCol = "Sugar_Sequence"

    # 解析全库 (Parse all records)
    topologies = df[seqCol].apply(parseSugarTopology)
    df["Total_Sugar_Count"] = topologies.apply(lambda t: t["total"])
    df["Sub_Chain_Lengths"] = topologies.apply(
        lambda t: "|".join(map(str, t["chains"])))
    df["Max_Chain_Length"] = topologies.apply(lambda t: t["max"])

    # 统计
    totalCounts = df["Total_Sugar_Count"].astype(int)
    maxChains = df["Max_Chain_Length"].astype(int)

    print(f"\n  Total Sugar Count distribution:")
    print(f"    Mean:   {totalCounts.mean():.2f}")
    print(f"    Median: {totalCounts.median():.1f}")
    print(f"    Max:    {totalCounts.max()}")
    print(f"    == 0:   {(totalCounts == 0).sum():,}")
    print(f"    == 1:   {(totalCounts == 1).sum():,}")
    print(f"    == 2:   {(totalCounts == 2).sum():,}")
    print(f"    == 3:   {(totalCounts == 3).sum():,}")
    print(f"    >= 4:   {(totalCounts >= 4).sum():,}")
    print(f"    >= 6:   {(totalCounts >= 6).sum():,}")
    print(f"    >= 10:  {(totalCounts >= 10).sum():,}")

    # === Chart A: 单糖链长度分布 (Sub-chain length distribution) ===
    allChainLens = []
    for chainStr in df["Sub_Chain_Lengths"]:
        if str(chainStr) in ("", "nan", "None"):
            continue
        for l in str(chainStr).split("|"):
            if l.isdigit() and int(l) > 0:
                allChainLens.append(int(l))

    chainLenDist = Counter(allChainLens)

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))

    # Chart A
    ax = axes[0]
    maxLen = min(max(chainLenDist.keys(), default=1), 15)
    xLabels = list(range(1, maxLen + 1))
    yValues = [chainLenDist.get(i, 0) for i in xLabels]

    barLabels = [f"Mono\n(1)", "Di\n(2)", "Tri\n(3)", "Tetra\n(4)",
                 "Penta\n(5)"]
    barLabels += [f"{i}" for i in range(6, maxLen + 1)]
    barLabels = barLabels[:len(xLabels)]

    colors = sns.color_palette("viridis", len(xLabels))
    bars = ax.bar(range(len(xLabels)), yValues, color=colors, edgecolor="white",
                  linewidth=0.5)
    ax.set_xticks(range(len(xLabels)))
    ax.set_xticklabels(barLabels, fontsize=9)
    ax.set_xlabel("Sugar Chain Length (per attachment site)", fontsize=11)
    ax.set_ylabel("Frequency", fontsize=11)
    ax.set_title("Chart A: Sub-Chain Length Distribution\n"
                 "(Each attachment site counted independently)",
                 fontsize=12, fontweight="bold")
    # 数值标注
    for bar, val in zip(bars, yValues):
        if val > 0:
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                    f"{val:,}", ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0, max(yValues) * 1.15 if yValues else 1)

    # Chart B: 每化合物总糖数分布
    ax2 = axes[1]
    totalCountDist = totalCounts[totalCounts > 0].value_counts().sort_index()
    maxTotal = min(totalCountDist.index.max() if len(totalCountDist) > 0 else 1,
                   20)
    xLabelsB = list(range(1, maxTotal + 1))
    yValuesB = [int(totalCountDist.get(i, 0)) for i in xLabelsB]

    colors2 = sns.color_palette("magma", len(xLabelsB))
    bars2 = ax2.bar(range(len(xLabelsB)), yValuesB, color=colors2,
                    edgecolor="white", linewidth=0.5)
    ax2.set_xticks(range(len(xLabelsB)))
    ax2.set_xticklabels([str(i) for i in xLabelsB], fontsize=9)
    ax2.set_xlabel("Total Sugar Count per Compound", fontsize=11)
    ax2.set_ylabel("Frequency", fontsize=11)
    ax2.set_title("Chart B: Total Sugar Load Distribution\n"
                  "(All sugars per compound summed)",
                  fontsize=12, fontweight="bold")
    for bar, val in zip(bars2, yValuesB):
        if val > 0:
            ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                     f"{val:,}", ha="center", va="bottom", fontsize=8)
    ax2.set_ylim(0, max(yValuesB) * 1.15 if yValuesB else 1)

    plt.tight_layout()
    chartPath = os.path.join(REPORT_DIR, "Sugar_Topology_Charts.png")
    fig.savefig(chartPath, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"\n  Charts saved: {chartPath}")

    return df


# =====================================================================
# Phase 5: ML-Ready Dataset Export
# =====================================================================

def runPhase5(df: pd.DataFrame) -> pd.DataFrame:
    """
    Phase 5: 导出 ML-Ready 数据集。
    Export GlycoNP_DeepChem_Ready.csv
    """
    print("\n" + "=" * 70)
    print("  Phase 5: ML-Ready Dataset Export")
    print("  机器学习数据集导出")
    print("=" * 70)

    # 选择核心列 (Select core columns)
    coreColumns = [
        "coconut_id",
        "canonical_smiles",
        "standard_inchi_key",
        "name",
        "organisms",
        "Sugar_Sequence",
        "Consensus_Sugar_Sequence",
        "Murcko_Scaffold",
        "Super_Scaffold_Class",
        "Detailed_NP_Class",
        "np_classifier_pathway",
        "Total_Sugar_Count",
        "Max_Chain_Length",
        "Sub_Chain_Lengths",
        "NLP_Bioactivity_Profile",
        "NLP_Targets",
        "LOTUS_kingdom",
        "LOTUS_family",
    ]

    # 添加可选的活性列 (Add optional bioactivity columns)
    optionalCols = [
        "bioactivity_summary",
        "pchembl_value",
        "chembl_target_name",
        "Is_Stereo_Upgraded",
        "Is_3D_Rescued",
    ]
    for col in optionalCols:
        if col in df.columns:
            coreColumns.append(col)

    # 仅保留存在的列 (Only keep existing columns)
    existingCols = [c for c in coreColumns if c in df.columns]

    mlDf = df[existingCols].copy()

    # 数据质量标注 (Data quality flags)
    mlDf["Has_Precise_Sugar"] = (~mlDf["Sugar_Sequence"].str.contains(
        r'\b(?:Hex|Pen|dHex|HexA|Non|Oct|Hept)\b',
        na=True, regex=True)).astype(int)

    mlDf["Has_Bioactivity"] = (
        mlDf.get("bioactivity_summary", pd.Series(dtype=str)).notna() &
        ~mlDf.get("bioactivity_summary", pd.Series(dtype=str)).astype(str
            ).isin(["", "nan", "None"])
    ).astype(int)

    # 导出全量 (Export full dataset)
    fullPath = os.path.join(REPORT_DIR, "GlycoNP_DeepChem_Ready.csv")
    mlDf.to_csv(fullPath, index=False, encoding="utf-8-sig")
    print(f"  Full dataset: {fullPath}")
    print(f"    Rows: {len(mlDf):,}")
    print(f"    Columns: {len(mlDf.columns)}")

    # 导出高活性子集 (Export high-potency subset)
    if "pchembl_value" in mlDf.columns:
        try:
            mlDf["_pchembl_float"] = pd.to_numeric(
                mlDf["pchembl_value"], errors="coerce")
            highPotency = mlDf[mlDf["_pchembl_float"] > 7.0].drop(
                columns=["_pchembl_float"])
            if len(highPotency) > 0:
                hpPath = os.path.join(REPORT_DIR,
                                      "GlycoNP_HighPotency_pChEMBL7.csv")
                highPotency.to_csv(hpPath, index=False, encoding="utf-8-sig")
                print(f"\n  High-potency subset (pChEMBL > 7):")
                print(f"    Rows: {len(highPotency):,}")
                print(f"    Saved: {hpPath}")
        except Exception:
            pass

    # 数据质量摘要 (Quality summary)
    print(f"\n  Data Quality Summary:")
    print(f"    Has precise sugar:   "
          f"{(mlDf['Has_Precise_Sugar']==1).sum():,} / {len(mlDf):,}")
    print(f"    Has bioactivity:     "
          f"{(mlDf['Has_Bioactivity']==1).sum():,} / {len(mlDf):,}")
    print(f"    Has scaffold class:  "
          f"{(mlDf['Super_Scaffold_Class']!='').sum():,} / {len(mlDf):,}")
    print(f"    Has sugar topology:  "
          f"{(mlDf['Total_Sugar_Count'].astype(int)>0).sum():,} / "
          f"{len(mlDf):,}")

    return df


# =====================================================================
# Main
# =====================================================================

def main():
    print("=" * 70)
    print("  GlycoNP Deep Architecture Upgrade — Phases 2 + 3 + 5")
    print("  深度架构升级 — 快速通道")
    print("=" * 70)

    # 自动选取最新数据 (Auto-select latest)
    candidates = [
        "GlycoNP_NLP_Enriched.csv",
        "GlycoNP_LLM_Rescued.csv",
        "GlycoNP_Fully_Enriched.csv",
    ]
    inputCsv = None
    for c in candidates:
        p = os.path.join(REPORT_DIR, c)
        if os.path.exists(p):
            inputCsv = p
            break

    if not inputCsv:
        print("  [ERROR] No input CSV found")
        return

    df = pd.read_csv(inputCsv, dtype=str, low_memory=False,
                     encoding="utf-8-sig")
    print(f"  Loaded: {len(df):,} rows from {os.path.basename(inputCsv)}")

    # ============ Phase 2 ============
    df = runPhase2(df)

    # ============ Phase 3 ============
    df = runPhase3(df)

    # ============ Save enriched ============
    enrichedPath = os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched.csv")
    df.to_csv(enrichedPath, index=False, encoding="utf-8-sig")
    print(f"\n  Enriched dataset: {enrichedPath}")

    # ============ Phase 5 ============
    runPhase5(df)

    print("\n" + "=" * 70)
    print("  Deep Architecture Upgrade Complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
