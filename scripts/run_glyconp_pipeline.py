"""
GlycoNP Pipeline — 全量管线 (Phase 1-7)
GlycoNP Pipeline — Full Pipeline (Phase 1-7)

Phase 1: 糖苷过滤 (使用已有的 Coconut_Sugar_Check.csv)
Phase 2: 糖苷/苷元切分 (glycosidic_cleavage)
Phase 3: 核苷酸/肽键二级扫描 (phase3_secondary_scan)
Phase 4: 分类学填补 (local_taxonomy — 已预处理)
Phase 5: 骨架/指纹/糖序列提取 (phase5_features)
Phase 6: 智能化学分类 (phase6_classifier)
Phase 7: 可视化报告 (phase7_visualizer — 可选)

使用方法 (Usage):
  python run_glyconp_pipeline.py [--limit N] [--output PATH] [--report]
"""
import argparse
import os
import sys
import time

import pandas as pd
from tqdm import tqdm

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from rdkit import Chem

from lib.glycan_topology import find_mapped_sugar_units
from lib.bond_cleavage_engine import cleaveWithConservation
from lib.secondary_fragment_scanner import scanSecondaryFragments
from lib.feature_extractor import characterizeGlycan, extractMurckoScaffold
from lib.chemical_classifier import classifyAglycon, buildReferenceFingerprints
from lib.modification_scanner import scanAndFormat as scanGlycanMods


# =====================================================================
# 核心行处理函数 (Core per-row processing function)
# =====================================================================

def processRow(smiles: str, refFps, refLabels) -> dict:
    """
    对单行 SMILES 执行 Phase 2→6 全流程。
    Execute Phase 2→6 on a single SMILES.

    Args:
        smiles: canonical_smiles
        refFps: Tanimoto 参考指纹列表
        refLabels: 参考标签列表

    Returns:
        dict with all Phase output columns
    """
    result = {
        "Glycan_SMILES": "", "Aglycon_SMILES": "",
        "Sugar_Sequence": "", "Glycan_Modifications": "",
        "Has_Nucleotide": False, "Nucleotide_Detail": "",
        "Has_Peptide": False, "Peptide_Detail": "",
        "Murcko_Scaffold": "", "Superclass": "",
        "Classification_Method": "", "Glycolipid_Flag": "",
    }

    if not smiles or smiles in ("NULL", "nan", ""):
        return result

    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return result

        # Phase 2: 切分 (Cleavage)
        units = find_mapped_sugar_units(mol)
        glycan, aglycon, meta = cleaveWithConservation(mol, units)
        result["Glycan_SMILES"] = glycan
        result["Aglycon_SMILES"] = aglycon

        # Phase 3: 核苷酸/肽键扫描 (Secondary scan)
        secResult = scanSecondaryFragments(smiles)
        result["Has_Nucleotide"] = secResult["Has_Nucleotide"]
        result["Nucleotide_Detail"] = secResult["Nucleotide_Detail"]
        result["Has_Peptide"] = secResult["Has_Peptide"]
        result["Peptide_Detail"] = secResult["Peptide_Detail"]

        # Phase 5: 糖序列 + 骨架 (Sugar sequence + scaffold)
        glycanFeats = characterizeGlycan(glycan)
        result["Sugar_Sequence"] = glycanFeats.get("sugar_sequence", "")
        result["Murcko_Scaffold"] = extractMurckoScaffold(aglycon)

        # 糖链修饰扫描 (Glycan modification scan — ONLY on Glycan_SMILES)
        result["Glycan_Modifications"] = scanGlycanMods(glycan)

        # Phase 6: 智能分类 (Classification)
        classResult = classifyAglycon(smiles, refFps, refLabels)
        superclass = classResult.get("classification", "Unclassified")
        result["Classification_Method"] = classResult.get("classification_method", "")
        result["Glycolipid_Flag"] = classResult.get("glycolipid_flag", "")

        # 强制覆盖: 核苷酸糖/糖肽 (Override: nucleotide sugar / glycopeptide)
        if result["Has_Nucleotide"]:
            superclass = "Nucleotide Sugar"
        elif result["Has_Peptide"]:
            superclass = "Glycopeptide"

        result["Superclass"] = superclass

    except Exception as e:
        result["Superclass"] = f"Error: {str(e)[:50]}"

    return result


# =====================================================================
# 主管线 (Main Pipeline)
# =====================================================================

def runPipeline(
    inputPath: str,
    outputPath: str,
    limit: int = 0,
    generateReport: bool = False,
):
    """
    执行 GlycoNP 全量管线。
    Execute the full GlycoNP pipeline.

    Args:
        inputPath: 输入 CSV (Coconut_Sugar_Check.csv)
        outputPath: 输出 CSV
        limit: 最大行数 (0=全量)
        generateReport: 是否生成 HTML 可视化报告
    """
    print("=" * 70)
    print("  GlycoNP Pipeline — Full Execution")
    print("=" * 70)
    t0 = time.time()

    # Load data
    print(f"\n[Phase 1] Loading data: {inputPath}")
    df = pd.read_csv(inputPath, low_memory=False, dtype=str, encoding="utf-8-sig")
    print(f"  Total rows: {len(df)}")

    if limit > 0:
        df = df.head(limit)
        print(f"  Limited to: {limit} rows")

    smilesCol = "canonical_smiles"
    classCol = "np_classifier_superclass"

    # Build Tanimoto reference from existing classifications
    print(f"\n[Phase 6 Prep] Building Tanimoto reference...")
    existingMask = (
        df[classCol].notna()
        & (~df[classCol].isin(["", "nan", "NULL"]))
        & df[smilesCol].notna()
    )
    refPool = df[existingMask]
    refSample = refPool.sample(n=min(500, len(refPool)), random_state=42) if len(refPool) > 0 else refPool
    refFps, refLabels = buildReferenceFingerprints(refSample, smilesCol, classCol)
    print(f"  Reference FPs: {len(refFps)}")

    # Process all rows
    print(f"\n[Phase 2→6] Processing {len(df)} compounds...")
    tProc = time.time()

    allResults = []
    for _, row in tqdm(df.iterrows(), total=len(df), desc="Pipeline"):
        smiles = str(row.get(smilesCol, ""))
        r = processRow(smiles, refFps, refLabels)
        allResults.append(r)

    # Append new columns
    resultDf = pd.DataFrame(allResults)
    for col in resultDf.columns:
        df[col] = resultDf[col].values

    procTime = time.time() - tProc
    print(f"  Processing time: {procTime:.1f}s ({procTime/len(df)*1000:.1f}ms/compound)")

    # Save output
    print(f"\n[Output] Saving to: {outputPath}")
    df.to_csv(outputPath, index=False, encoding="utf-8-sig")
    print(f"  Saved {len(df)} rows, {len(df.columns)} columns")

    # Statistics
    print(f"\n{'='*70}")
    print("  Pipeline Statistics")
    print(f"{'='*70}")

    nucCount = df["Has_Nucleotide"].sum() if "Has_Nucleotide" in df else 0
    pepCount = df["Has_Peptide"].sum() if "Has_Peptide" in df else 0
    modCount = sum(1 for v in df["Glycan_Modifications"] if v and v != "")
    seqCount = sum(1 for v in df["Sugar_Sequence"] if v and v != "" and v != "nan")

    classDist = df["Superclass"].value_counts().head(10)
    print(f"\n  Nucleotide Sugars: {nucCount}")
    print(f"  Glycopeptides: {pepCount}")
    print(f"  With Glycan Modifications: {modCount}")
    print(f"  Sugar Sequences generated: {seqCount}")
    print(f"\n  [Superclass Distribution (top 10)]")
    for cls, count in classDist.items():
        print(f"    {str(cls):45s} {count:5d}")

    # Optional HTML report
    if generateReport:
        reportPath = outputPath.replace(".csv", "_Report.html")
        print(f"\n[Phase 7] Generating HTML report: {reportPath}")
        from lib.molecular_visualizer import generateHtmlReport
        generateHtmlReport(
            df, reportPath,
            smilesCol=smilesCol,
            inchikeyCol="standard_inchi_key",
            sugarSeqCol="Sugar_Sequence",
            scaffoldCol="Murcko_Scaffold",
            classCol="Superclass",
            glycolipidCol="Glycolipid_Flag",
            maxRows=100,
        )

    totalTime = time.time() - t0
    print(f"\n  Total time: {totalTime:.1f}s")
    print(f"{'='*70}")

    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GlycoNP Pipeline")
    parser.add_argument("--limit", type=int, default=0, help="Max rows (0=all)")
    parser.add_argument("--output", type=str, default=None, help="Output CSV path")
    parser.add_argument("--report", action="store_true", help="Generate HTML report")
    args = parser.parse_args()

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    inputPath = os.path.join(baseDir, "reports", "Coconut_Sugar_Check.csv")

    if args.output:
        outputPath = args.output
    else:
        suffix = f"_{args.limit}" if args.limit > 0 else "_Full"
        outputPath = os.path.join(baseDir, "reports", f"GlycoNP_Pipeline{suffix}.csv")

    runPipeline(inputPath, outputPath, limit=args.limit, generateReport=args.report)
