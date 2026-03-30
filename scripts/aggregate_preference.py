"""
偏好性聚合脚本 — 生成 "物种×糖链" 和 "骨架×糖链" 频次交叉表
Preference Aggregation — Generate organism×sugar and scaffold×sugar crosstabs

输出文件 / Output Files:
  1. reports/organism_sugar_matrix.csv — 物种 (Kingdom/Phylum/Family) × 糖链序列 交叉表
  2. reports/aglycone_sugar_matrix.csv — Murcko 骨架 (或大类) × 糖链序列 交叉表

这两个矩阵可直接用于绘制热力图和桑基图 (Sankey Diagram)。
These matrices can be directly used for heatmaps and Sankey diagrams.

使用方法 / Usage:
  python scripts/aggregate_preference.py
  python scripts/aggregate_preference.py --top-sugars 30 --top-organisms 50
"""
import argparse
import os
import re
import sys
import time
from typing import Optional

import pandas as pd
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


# =====================================================================
# 1. 物种分类学层级提取 (Taxonomy Level Extraction)
# =====================================================================

def extractTaxonomy(
    organismStr: str,
    level: str = "family",
) -> str:
    """
    从 organisms 字段中提取分类学层级。
    Extract taxonomy level from organisms field.

    策略: organisms 字段格式多样 (COCONUT 原始格式)。
    常见形式:
      - "Glycyrrhiza glabra" (直接学名)
      - "Glycyrrhiza glabra|Glycyrrhiza uralensis" (多物种, | 分隔)
      - "Fabaceae" (已为科级别)

    对于无法解析的, 返回 "Unknown"。
    """
    if not organismStr or str(organismStr) in ("nan", "", "None"):
        return "Unknown"

    # 取第一个物种 (Take first organism if multiple)
    org = str(organismStr).split("|")[0].strip()
    if not org or org == "nan":
        return "Unknown"

    # 如果已经看起来像科名 (以 -aceae, -idae 结尾)
    # If already looks like a family name
    if re.search(r'(aceae|idae|ales|phyta|mycota)$', org, re.IGNORECASE):
        return org

    # 返回属名 (genus) 作为近似分类
    # Return genus name as approximate taxonomy
    parts = org.strip().split()
    if parts:
        return parts[0]  # 属名 (genus name)

    return "Unknown"


def extractKingdom(organismStr: str) -> str:
    """
    根据物种名推断界级别 (粗分类)。
    Infer kingdom level from organism name (coarse classification).

    使用简单的后缀/关键词规则:
      - *phyta, *mycota → 对应界
      - 已知属名映射
    """
    if not organismStr or str(organismStr) in ("nan", "", "None"):
        return "Unknown"

    org = str(organismStr).lower().strip()

    # 关键词推断 (Keyword inference)
    if any(k in org for k in ["mycota", "fungus", "fungi", "aspergillus",
                               "penicillium", "fusarium", "trichoderma"]):
        return "Fungi"
    if any(k in org for k in ["bacteria", "streptomyces", "bacillus",
                               "pseudomonas", "escherichia", "actinomyces"]):
        return "Bacteria"
    if any(k in org for k in ["sponge", "coral", "ascidian", "nudibranch",
                               "starfish", "sea urchin", "holothuria"]):
        return "Animalia(Marine)"

    # 默认植物 (Default: Plantae — 大部分天然产物来自植物)
    return "Plantae"


# =====================================================================
# 2. 糖链序列标准化 (Sugar Sequence Normalization)
# =====================================================================

def normalizeSugarSequence(
    seqStr: str,
    maxTokens: int = 4,
) -> str:
    """
    标准化糖链序列以减少交叉表维度。
    Normalize sugar sequence to reduce crosstab dimensionality.

    策略 / Strategy:
      1. 截断到 maxTokens 个糖单元 (避免超长序列占比太低)
      2. 移除空白和多余分隔符
      3. 统一分隔符为 " → "
    """
    if not seqStr or str(seqStr) in ("nan", "", "None"):
        return "Unknown"

    seq = str(seqStr).strip()

    # 可能的分隔符 (Possible delimiters)
    for sep in [" → ", "→", " -> ", "->", " - ", ","]:
        if sep in seq:
            tokens = [t.strip() for t in seq.split(sep) if t.strip()]
            if len(tokens) > maxTokens:
                tokens = tokens[:maxTokens]
                tokens.append("...")
            return " → ".join(tokens)

    return seq


# =====================================================================
# 3. 交叉频次表生成 (Crosstab Generation)
# =====================================================================

def buildOrganismSugarCrosstab(
    df: pd.DataFrame,
    topNOrganisms: int = 50,
    topNSugars: int = 30,
) -> pd.DataFrame:
    """
    生成 "物种 × 糖链" 频次交叉表。
    Generate organism × sugar sequence frequency crosstab.

    行: 物种 (属/科级别), 列: 糖链序列, 值: 出现次数
    Rows: organism (genus/family), Columns: sugar sequence, Values: count
    """
    print("  Building organism × sugar crosstab...")

    # 准备数据 (Prepare data)
    workDf = df[["organisms", "Final_Sugar_Sequence"]].copy()
    workDf = workDf.dropna(subset=["Final_Sugar_Sequence"])
    workDf = workDf[workDf["Final_Sugar_Sequence"] != ""]

    # 提取分类学 (Extract taxonomy)
    workDf["organism_group"] = workDf["organisms"].apply(extractTaxonomy)
    workDf["sugar_norm"] = workDf["Final_Sugar_Sequence"].apply(
        normalizeSugarSequence)

    # 过滤 Unknown (Filter Unknown)
    workDf = workDf[(workDf["organism_group"] != "Unknown") &
                     (workDf["sugar_norm"] != "Unknown")]

    if workDf.empty:
        print("    [WARN] No valid data for crosstab")
        return pd.DataFrame()

    # Top N 过滤 (Top N filtering)
    topOrganisms = workDf["organism_group"].value_counts().head(topNOrganisms).index
    topSugars = workDf["sugar_norm"].value_counts().head(topNSugars).index

    filteredDf = workDf[
        workDf["organism_group"].isin(topOrganisms) &
        workDf["sugar_norm"].isin(topSugars)
    ]

    # 交叉表 (Crosstab)
    crosstab = pd.crosstab(
        filteredDf["organism_group"],
        filteredDf["sugar_norm"],
    )

    print(f"    Shape: {crosstab.shape[0]} organisms × {crosstab.shape[1]} sugars")
    print(f"    Non-zero cells: {(crosstab > 0).sum().sum():,} / "
          f"{crosstab.shape[0] * crosstab.shape[1]:,}")

    return crosstab


def buildScaffoldSugarCrosstab(
    df: pd.DataFrame,
    topNScaffolds: int = 50,
    topNSugars: int = 30,
    scaffoldCol: str = "Murcko_Scaffold",
) -> pd.DataFrame:
    """
    生成 "骨架 × 糖链" 频次交叉表。
    Generate scaffold × sugar sequence frequency crosstab.

    行: Murcko 骨架 (或 Superclass), 列: 糖链序列, 值: 出现次数
    Rows: Murcko scaffold (or Superclass), Columns: sugar sequence, Values: count
    """
    print("  Building scaffold × sugar crosstab...")

    # 如果 Murcko_Scaffold 缺失或全空, 回退到 Superclass
    # If Murcko_Scaffold missing or all empty, fallback to Superclass
    useCol = scaffoldCol
    if useCol not in df.columns:
        useCol = "Superclass"
        print(f"    [INFO] '{scaffoldCol}' not found, using '{useCol}'")

    validCount = df[useCol].notna().sum()
    if validCount < 100:
        useCol = "Superclass"
        print(f"    [INFO] Too few valid '{scaffoldCol}', using '{useCol}'")

    # 准备数据 (Prepare data)
    workDf = df[[useCol, "Final_Sugar_Sequence"]].copy()
    workDf = workDf.dropna(subset=["Final_Sugar_Sequence", useCol])
    workDf = workDf[(workDf["Final_Sugar_Sequence"] != "") &
                     (workDf[useCol] != "")]

    # 清理 Superclass (Clean Superclass if applicable)
    if useCol == "Superclass":
        workDf[useCol] = workDf[useCol].apply(
            lambda x: re.sub(r'\(Tanimoto=.*?\)', '', str(x)).strip()
        )

    # 标准化糖链序列 (Normalize sugar sequences)
    workDf["sugar_norm"] = workDf["Final_Sugar_Sequence"].apply(
        normalizeSugarSequence)
    workDf = workDf[workDf["sugar_norm"] != "Unknown"]

    if workDf.empty:
        print("    [WARN] No valid data for crosstab")
        return pd.DataFrame()

    # Top N 过滤 (Top N filtering)
    topScaffolds = workDf[useCol].value_counts().head(topNScaffolds).index
    topSugars = workDf["sugar_norm"].value_counts().head(topNSugars).index

    filteredDf = workDf[
        workDf[useCol].isin(topScaffolds) &
        workDf["sugar_norm"].isin(topSugars)
    ]

    # 交叉表 (Crosstab)
    crosstab = pd.crosstab(
        filteredDf[useCol],
        filteredDf["sugar_norm"],
    )

    scaffoldLabel = "scaffolds" if useCol == scaffoldCol else "superclasses"
    print(f"    Shape: {crosstab.shape[0]} {scaffoldLabel} × "
          f"{crosstab.shape[1]} sugars")
    print(f"    Non-zero cells: {(crosstab > 0).sum().sum():,} / "
          f"{crosstab.shape[0] * crosstab.shape[1]:,}")

    return crosstab


# =====================================================================
# 4. Final_Sugar_Sequence 列构建
# Build Final_Sugar_Sequence column
# =====================================================================

def buildFinalSugarSequence(df: pd.DataFrame) -> pd.DataFrame:
    """
    优先使用 Rescued_Sugar_Sequence, 回退到 Sugar_Sequence。
    Prefer Rescued_Sugar_Sequence, fallback to Sugar_Sequence.
    """
    def _pick(row):
        rescued = str(row.get("Rescued_Sugar_Sequence", ""))
        original = str(row.get("Sugar_Sequence", ""))
        if rescued and rescued not in ("", "nan", "None"):
            return rescued
        if original and original not in ("", "nan", "None"):
            return original
        return ""

    df["Final_Sugar_Sequence"] = df.apply(_pick, axis=1)
    validCount = (df["Final_Sugar_Sequence"] != "").sum()
    print(f"  Final_Sugar_Sequence: {validCount:,} / {len(df):,} valid")
    return df


# =====================================================================
# Main
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Preference Aggregation Crosstab Generator")
    parser.add_argument("--input", type=str, default=None,
                        help="Enriched GlycoNP CSV")
    parser.add_argument("--top-sugars", type=int, default=30)
    parser.add_argument("--top-organisms", type=int, default=50)
    parser.add_argument("--top-scaffolds", type=int, default=50)
    args = parser.parse_args()

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")

    # 尝试多个可能的输入路径 (Try multiple input paths)
    # 优先级: Fully_Enriched > ChEMBL_Enriched > Pipeline_Full_Cleaned
    inputPath = args.input
    if not inputPath:
        candidates = [
            os.path.join(reportDir, "GlycoNP_Fully_Enriched.csv"),
            os.path.join(reportDir, "GlycoNP_ChEMBL_Enriched.csv"),
            os.path.join(reportDir, "GlycoNP_Pipeline_Full_Cleaned.csv"),
        ]
        for path in candidates:
            if os.path.exists(path):
                inputPath = path
                break

    if not inputPath or not os.path.exists(inputPath):
        print(f"  [ERROR] No input CSV found. Run pipeline first.")
        return

    print("=" * 70)
    print("  GlycoNP Preference Aggregation — Crosstab Generator")
    print("  糖缀合物偏好性聚合 — 交叉频次表生成器")
    print("=" * 70)
    print(f"  Input: {inputPath}")

    t0 = time.time()
    df = pd.read_csv(inputPath, low_memory=False, dtype=str,
                     encoding="utf-8-sig")
    print(f"  Loaded: {len(df):,} rows")

    # 构建 Final_Sugar_Sequence (Build Final_Sugar_Sequence)
    df = buildFinalSugarSequence(df)

    # ---- Crosstab 1: Organism × Sugar ----
    print(f"\n{'='*70}")
    print(f"  Crosstab 1: Organism × Sugar Sequence")
    print(f"{'='*70}")

    orgSugarMatrix = buildOrganismSugarCrosstab(
        df, topNOrganisms=args.top_organisms, topNSugars=args.top_sugars,
    )

    orgMatrixPath = os.path.join(reportDir, "organism_sugar_matrix.csv")
    if not orgSugarMatrix.empty:
        orgSugarMatrix.to_csv(orgMatrixPath, encoding="utf-8-sig")
        print(f"  Saved: {orgMatrixPath}")
    else:
        print(f"  [SKIP] Empty organism matrix")

    # ---- Crosstab 2: Scaffold × Sugar ----
    print(f"\n{'='*70}")
    print(f"  Crosstab 2: Scaffold/Superclass × Sugar Sequence")
    print(f"{'='*70}")

    scaffoldSugarMatrix = buildScaffoldSugarCrosstab(
        df, topNScaffolds=args.top_scaffolds, topNSugars=args.top_sugars,
    )

    scaffoldMatrixPath = os.path.join(reportDir, "aglycone_sugar_matrix.csv")
    if not scaffoldSugarMatrix.empty:
        scaffoldSugarMatrix.to_csv(scaffoldMatrixPath, encoding="utf-8-sig")
        print(f"  Saved: {scaffoldMatrixPath}")
    else:
        print(f"  [SKIP] Empty scaffold matrix")

    # ---- Summary ----
    elapsed = time.time() - t0
    print(f"\n{'='*70}")
    print(f"  Aggregation Complete!")
    print(f"  Time: {elapsed:.0f}s")
    if not orgSugarMatrix.empty:
        print(f"  1. {orgMatrixPath} "
              f"({orgSugarMatrix.shape[0]}×{orgSugarMatrix.shape[1]})")
    if not scaffoldSugarMatrix.empty:
        print(f"  2. {scaffoldMatrixPath} "
              f"({scaffoldSugarMatrix.shape[0]}×{scaffoldSugarMatrix.shape[1]})")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
