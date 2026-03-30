"""
GlycoNP 数据体检与进度报告 (Data Quality Audit & Roadmap Report)

四维度审计 (Four-Dimension Audit):
  1. 核心数据缺失率 (Missing Values)
  2. 糖链解析精度 — Hex 占比 (Sugar Precision)
  3. 异常碎片排查 — 苷元重原子数 (Aglycon Sanity)
  4. 功能开发进度盘点 (Roadmap Check)

使用方法 (Usage):
  python scripts/audit_data_quality.py [--input PATH] [--imputed]
"""
import argparse
import os
import re
import sys
from collections import Counter
from typing import Dict, List, Tuple

import pandas as pd
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


def main():
    parser = argparse.ArgumentParser(
        description="GlycoNP Data Quality Audit")
    parser.add_argument("--input", type=str, default=None)
    parser.add_argument("--imputed", action="store_true",
                        help="Use imputed version")
    args = parser.parse_args()

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")

    if args.input:
        inputPath = args.input
    elif args.imputed:
        inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full_Imputed.csv")
    else:
        inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full.csv")

    df = pd.read_csv(inputPath, low_memory=False, dtype=str,
                     encoding="utf-8-sig")
    total = len(df)

    mdLines: List[str] = []
    mdLines.append("# GlycoNP Data Quality Audit Report")
    mdLines.append("# GlycoNP 数据质量体检报告\n")
    mdLines.append(f"> Source: `{os.path.basename(inputPath)}`")
    mdLines.append(f"> Total compounds: **{total:,}**\n")

    print("=" * 70)
    print("  GlycoNP Data Quality Audit")
    print("=" * 70)
    print(f"  Input: {inputPath}")
    print(f"  Total: {total:,}\n")

    # =================================================================
    # 1. 核心数据缺失率 (Missing Values Report)
    # =================================================================
    print("  [1/4] Missing Values Report")
    mdLines.append("## 1. Missing Values Report / 核心数据缺失率\n")

    missingCols = [
        ("organisms", "物种来源"),
        ("organism_taxonomy_05family", "科级分类 (Family)"),
        ("np_classifier_superclass", "NP Classifier 大类"),
        ("Superclass", "Superclass (Pipeline)"),
        ("dois", "文献链接 (DOI)"),
        ("Sugar_Sequence", "糖链序列"),
        ("Glycan_Modifications", "糖链修饰"),
        ("Aglycon_SMILES", "苷元 SMILES"),
        ("Murcko_Scaffold", "Murcko 骨架"),
        ("canonical_smiles", "SMILES"),
        ("alogp", "LogP"),
        ("np_likeness", "NP-Likeness"),
    ]

    mdLines.append("| Column | NaN Count | % | Status |")
    mdLines.append("|:-------|----------:|--:|:------:|")

    for col, desc in missingCols:
        if col in df.columns:
            # 检查真空值: NaN, 空字符串, "nan" 字符串
            mask = df[col].isna() | (df[col].astype(str).str.strip().isin(
                ["", "nan", "None", "NULL"]))
            nMissing = mask.sum()
            pct = nMissing / total * 100
            status = "✅" if pct < 5 else "⚠️" if pct < 30 else "🔴"
            print(f"    {col:40s} {nMissing:>6} ({pct:5.1f}%) {status}")
            mdLines.append(
                f"| `{col}` ({desc}) | {nMissing:,} | {pct:.1f}% | {status} |")
        else:
            print(f"    {col:40s} COLUMN MISSING ❌")
            mdLines.append(f"| `{col}` ({desc}) | — | — | ❌ Missing |")

    # Superclass "Unclassified" 特别统计
    if "Superclass" in df.columns:
        unclassified = (df["Superclass"].astype(str).str.strip().str.lower()
                        .isin(["unclassified", ""]))
        nUncl = unclassified.sum()
        pctUncl = nUncl / total * 100
        mdLines.append(
            f"\n> **Superclass = 'Unclassified'**: {nUncl:,} ({pctUncl:.1f}%)")
        print(f"\n    Superclass='Unclassified': {nUncl:,} ({pctUncl:.1f}%)")

    # =================================================================
    # 2. 糖链解析精度 — Hex 比率 (Sugar Precision)
    # =================================================================
    print("\n  [2/4] Sugar Precision Report")
    mdLines.append("\n## 2. Sugar Precision / 糖链解析精度\n")

    if "Sugar_Sequence" in df.columns:
        seqSeries = df["Sugar_Sequence"].dropna().astype(str)
        seqSeries = seqSeries[~seqSeries.isin(["", "nan", "None"])]
        totalSugarRows = len(seqSeries)

        # 提取所有糖单元 (Extract all monosaccharide tokens)
        allTokens = []
        for seq in seqSeries:
            # 去除连接信息, 提取糖名 (Strip linkage info)
            tokens = re.findall(
                r'[DL]-[A-Z][a-z]+[A-Z]?[a-z]*|Hex|dHex|Pen|HexA|HexN|'
                r'HexNAc|Unk|Invalid|Error',
                seq)
            allTokens.extend(tokens)

        totalTokens = len(allTokens)

        # 分类: 精确 vs 退避 (Classify: precise vs fallback)
        FALLBACK_PATTERNS = {"Hex", "dHex", "Pen", "HexA", "HexN", "HexNAc",
                             "Unk", "Invalid", "Error"}

        preciseTokens = [t for t in allTokens if t not in FALLBACK_PATTERNS]
        fallbackTokens = [t for t in allTokens if t in FALLBACK_PATTERNS]

        nPrecise = len(preciseTokens)
        nFallback = len(fallbackTokens)
        preciseRatio = nPrecise / totalTokens * 100 if totalTokens else 0
        fallbackRatio = nFallback / totalTokens * 100 if totalTokens else 0

        print(f"    Rows with Sugar_Sequence: {totalSugarRows:,}")
        print(f"    Total monosaccharide tokens: {totalTokens:,}")
        print(f"    Precise (D-Glc, L-Rha, etc.): {nPrecise:,} ({preciseRatio:.1f}%)")
        print(f"    Fallback (Hex, dHex, Pen...): {nFallback:,} ({fallbackRatio:.1f}%)")

        mdLines.append(f"| Metric | Count | % |")
        mdLines.append(f"|:-------|------:|--:|")
        mdLines.append(f"| Rows with sugar | {totalSugarRows:,} | — |")
        mdLines.append(f"| Total tokens | {totalTokens:,} | 100% |")
        mdLines.append(
            f"| **Precise** (D-Glc, L-Rha, L-Col...) | {nPrecise:,} "
            f"| **{preciseRatio:.1f}%** |")
        mdLines.append(
            f"| **Fallback** (Hex, dHex, Pen...) | {nFallback:,} "
            f"| **{fallbackRatio:.1f}%** |")

        # Top 精确糖 (Top precise sugars)
        preciseCounter = Counter(preciseTokens)
        mdLines.append(f"\n### Top 15 Precise Sugars\n")
        mdLines.append("| Sugar | Count | % of precise |")
        mdLines.append("|:------|------:|-------------:|")
        print(f"\n    Top precise sugars:")
        for sugar, cnt in preciseCounter.most_common(15):
            pct = cnt / nPrecise * 100 if nPrecise else 0
            print(f"      {sugar:20s} {cnt:>6} ({pct:.1f}%)")
            mdLines.append(f"| `{sugar}` | {cnt:,} | {pct:.1f}% |")

        # Top 退避糖 (Top fallback)
        fallbackCounter = Counter(fallbackTokens)
        mdLines.append(f"\n### Fallback Breakdown\n")
        mdLines.append("| Fallback | Count | % of fallback |")
        mdLines.append("|:---------|------:|--------------:|")
        print(f"\n    Fallback breakdown:")
        for sugar, cnt in fallbackCounter.most_common():
            pct = cnt / nFallback * 100 if nFallback else 0
            print(f"      {sugar:20s} {cnt:>6} ({pct:.1f}%)")
            mdLines.append(f"| `{sugar}` | {cnt:,} | {pct:.1f}% |")

    # =================================================================
    # 3. 苷元重原子数排查 (Aglycon Sanity Check)
    # =================================================================
    print("\n  [3/4] Aglycon Sanity Check")
    mdLines.append("\n## 3. Aglycon Sanity Check / 异常碎片排查\n")

    if "Aglycon_SMILES" in df.columns:
        from rdkit import Chem

        heavyAtomCounts = []
        errorCount = 0

        for smiles in df["Aglycon_SMILES"]:
            s = str(smiles)
            if s in ("nan", "", "None", "NULL"):
                heavyAtomCounts.append(-1)
                continue
            try:
                mol = Chem.MolFromSmiles(s)
                if mol:
                    heavyAtomCounts.append(mol.GetNumHeavyAtoms())
                else:
                    heavyAtomCounts.append(-1)
                    errorCount += 1
            except Exception:
                heavyAtomCounts.append(-1)
                errorCount += 1

        df["_ha_count"] = heavyAtomCounts
        valid = df[df["_ha_count"] >= 0]

        # 重原子数分布 (Heavy atom count distribution)
        tinyAglycons = valid[valid["_ha_count"] <= 5]
        smallAglycons = valid[(valid["_ha_count"] > 5) & (valid["_ha_count"] <= 10)]
        normalAglycons = valid[valid["_ha_count"] > 10]

        nTiny = len(tinyAglycons)
        nSmall = len(smallAglycons)
        nNormal = len(normalAglycons)

        print(f"    Valid aglycon SMILES: {len(valid):,}")
        print(f"    Parse errors: {errorCount}")
        print(f"    Heavy atoms ≤ 5 (suspicious): {nTiny:,} ({nTiny/len(valid)*100:.1f}%)")
        print(f"    Heavy atoms 6-10 (small): {nSmall:,}")
        print(f"    Heavy atoms > 10 (normal): {nNormal:,}")

        mdLines.append("| Category | Count | % |")
        mdLines.append("|:---------|------:|--:|")
        mdLines.append(f"| Valid aglycon SMILES | {len(valid):,} | — |")
        mdLines.append(
            f"| **≤ 5 heavy atoms** (suspicious) | **{nTiny:,}** | "
            f"**{nTiny/len(valid)*100:.1f}%** |")
        mdLines.append(
            f"| 6–10 heavy atoms (small) | {nSmall:,} | "
            f"{nSmall/len(valid)*100:.1f}% |")
        mdLines.append(
            f"| > 10 heavy atoms (normal) | {nNormal:,} | "
            f"{nNormal/len(valid)*100:.1f}% |")

        # Top tiny aglycons
        if nTiny > 0:
            mdLines.append(f"\n### Top 10 Suspicious Tiny Aglycons (HA ≤ 5)\n")
            mdLines.append("| SMILES | HA | Count | Likely Issue |")
            mdLines.append("|:-------|---:|------:|:-------------|")
            tinySmilesCounts = tinyAglycons["Aglycon_SMILES"].value_counts()
            print(f"\n    Top suspicious tiny aglycons:")
            for smi, cnt in tinySmilesCounts.head(10).items():
                try:
                    mol = Chem.MolFromSmiles(str(smi))
                    ha = mol.GetNumHeavyAtoms() if mol else 0
                except Exception:
                    ha = 0
                issue = ("Methyl/Ethyl fragment" if ha <= 2
                         else "Small substituent" if ha <= 4
                         else "Borderline")
                print(f"      {str(smi)[:40]:40s} HA={ha} Count={cnt}")
                mdLines.append(f"| `{str(smi)[:50]}` | {ha} | {cnt:,} | {issue} |")

        # Histogram buckets
        mdLines.append(f"\n### Heavy Atom Count Distribution\n")
        mdLines.append("| Range | Count | % |")
        mdLines.append("|:------|------:|--:|")
        bins = [(0, 5), (6, 10), (11, 20), (21, 30), (31, 50), (51, 100), (101, 9999)]
        for low, high in bins:
            n = len(valid[(valid["_ha_count"] >= low) & (valid["_ha_count"] <= high)])
            label = f"{low}–{high}" if high < 9999 else f"{low}+"
            mdLines.append(f"| {label} | {n:,} | {n/len(valid)*100:.1f}% |")

        df.drop(columns=["_ha_count"], inplace=True)

    # =================================================================
    # 4. 功能开发进度盘点 (Roadmap Check)
    # =================================================================
    print("\n  [4/4] Roadmap Check")
    mdLines.append("\n## 4. Development Roadmap / 功能开发进度\n")

    completed = [
        ("Phase 1", "InChIKey 去重 + 含糖筛选", "✅"),
        ("Phase 2", "糖苷键切割 + 糖/苷元分离", "✅"),
        ("Phase 3", "核苷酸/肽键二级扫描", "✅ (peptide removed due to false positive)"),
        ("Phase 4", "LOTUS 分类学填充 + 名称退避", "✅"),
        ("Phase 5", "Morgan FP + Murcko Scaffold + 修饰基团", "✅"),
        ("Phase 6", "NP Classifier + Superclass", "✅"),
        ("Phase 7", "四色可视化 + Expert Review Panel", "✅"),
        ("UMAP", "化学空间聚类 (10K sample)", "✅"),
        ("Network", "科级网络图 (8 families)", "✅"),
        ("Sankey", "糖流桑基图 (Top 10 classes)", "✅"),
        ("Heatmap", "交叉热力图 (golden pairs)", "✅"),
        ("ChEMBL", "真实活性整合 (7,137 hits, 812 targets)", "✅"),
        ("Bio-Value", "糖链生物活性评分 (678 combos)", "✅"),
        ("Mod Landscape", "同糖不同修饰分析", "✅"),
        ("Multi-Sheet", "按 Superclass 分 Sheet Excel", "✅"),
        ("Imputation", "三层数据填补 (Tanimoto rescue)", "✅"),
    ]

    pending = [
        ("ΔLogP", "糖链剥离前后 LogP 差值计算", "📋 Not started"),
        ("Glycan Similarity", "糖链 Morgan 相似度矩阵", "📋 Not started"),
        ("Aglycon Similarity", "苷元 Tanimoto 相似度", "📋 Not started"),
        ("UMAP Markers", "按相似类型区分点形状", "📋 Not started"),
        ("OpenAlex 文献", "DOI → 摘要/关键词提取", "📋 Candidates prepared"),
        ("Organism API", "PubChem/GBIF 物种二级查询", "📋 Not started"),
        ("Phase 6 Image Export", "按分类层级存储图片", "📋 Not started"),
        ("SMARTS Depth Scan", "遗留测试文件 SMARTS 提取", "📋 Not started"),
    ]

    mdLines.append("### Completed Modules / 已完成\n")
    mdLines.append("| Phase | Module | Status |")
    mdLines.append("|:------|:-------|:------:|")
    for phase, desc, status in completed:
        mdLines.append(f"| {phase} | {desc} | {status} |")
        print(f"    {status} {phase}: {desc}")

    mdLines.append("\n### Pending / 待办\n")
    mdLines.append("| Module | Description | Status |")
    mdLines.append("|:-------|:------------|:------:|")
    for module, desc, status in pending:
        mdLines.append(f"| {module} | {desc} | {status} |")
        print(f"    {status} {module}: {desc}")

    # =================================================================
    # 5. Debug 优先级建议 (Recommendations)
    # =================================================================
    mdLines.append("\n## 5. Debug Priority Recommendations\n")
    mdLines.append("> Based on audit results:\n")
    mdLines.append(
        "1. **Sugar Precision**: If fallback (Hex/dHex/Pen) ratio is high, "
        "expand the monosaccharide reference library in `lib/monosaccharide_identifier.py`.\n")
    mdLines.append(
        "2. **Tiny Aglycons**: Compounds with HA ≤ 5 aglycons should be "
        "reviewed — they may be free sugar derivatives, not true glycosides.\n")
    mdLines.append(
        "3. **DOI Coverage**: 58% missing DOIs limit literature-based "
        "activity mining. Consider OpenAlex batch API or CrossRef.\n")
    mdLines.append(
        "4. **Organism Coverage**: 49% missing organisms limit "
        "taxonomy-based analysis. Name-based dictionary + PubChem API recommended.\n")

    # 保存报告 (Save report)
    mdPath = os.path.join(reportDir, "Data_Quality_Audit.md")
    with open(mdPath, "w", encoding="utf-8") as f:
        f.write("\n".join(mdLines))

    print(f"\n  Report saved: {mdPath}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
