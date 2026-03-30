"""
GlycoNP 手性完整度审计 (Chirality Completeness Audit)
=====================================================

量化 94K 含糖 SMILES 中立体化学信息的缺失程度。
Quantify stereochemistry deficiency across 94K sugar-containing SMILES.

审计维度 (Audit Dimensions):
  1. 分子级: 每个含糖分子的手性覆盖状态 (FULL_3D / PARTIAL_3D / FLAT_2D)
  2. 糖环级: 每个独立糖环的手性碳 (C2-C5) 定义率
  3. 关联分析: 手性缺失 vs 当前 Sugar_Sequence 中泛指标签 (Hex/Pen) 的相关性

使用方法 (Usage):
  python scripts/audit_stereo_completeness.py [--input PATH] [--limit N]
"""
import argparse
import os
import re
import sys
import time
from collections import Counter
from typing import Dict, List, Tuple, Optional

import pandas as pd
from rdkit import Chem

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from lib.glycan_topology import find_mapped_sugar_units


# =====================================================================
# 1. 核心审计引擎 (Core Audit Engine)
# =====================================================================

def auditSugarRingChirality(
    mol: Chem.Mol,
    ringAtoms: List[int],
) -> Dict:
    """审计单个糖环的手性完整度。
    Audit chirality completeness for a single sugar ring.

    分析环内碳原子 (排除环氧和异头碳) 的 ChiralTag:
    - FULL_3D: 所有非异头碳环碳有 @/@@
    - PARTIAL_3D: 部分有手性标记
    - FLAT_2D: 所有环碳都无手性标记

    Args:
        mol: RDKit Mol 对象
        ringAtoms: 糖环原子索引列表

    Returns:
        审计结果字典
    """
    ringSet = set(ringAtoms)

    # 找到环氧 (Ring oxygen)
    ringOxygenIdx = None
    for idx in ringAtoms:
        if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8:
            ringOxygenIdx = idx
            break

    # 收集环碳 (Ring carbons, excluding ring O)
    ringCarbons = [idx for idx in ringAtoms
                   if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]

    if not ringCarbons:
        return {"status": "NO_CARBON", "defined": 0, "total": 0}

    # 异头碳检测: 连环O + 环外杂原子的碳
    # Anomeric carbon: ring C adjacent to ring O AND exocyclic heteroatom
    anomericCarbonIdx = None
    if ringOxygenIdx is not None:
        for nbr in mol.GetAtomWithIdx(ringOxygenIdx).GetNeighbors():
            if nbr.GetIdx() in ringSet and nbr.GetAtomicNum() == 6:
                hasExoHetero = any(
                    n.GetAtomicNum() in (7, 8, 16) and n.GetIdx() not in ringSet
                    for n in nbr.GetNeighbors()
                )
                if hasExoHetero:
                    anomericCarbonIdx = nbr.GetIdx()
                    break

    # 只检查非异头碳的环碳 (C2-C5) 的手性
    # Only check chirality on non-anomeric ring carbons (C2-C5)
    checkableCarbons = [idx for idx in ringCarbons if idx != anomericCarbonIdx]

    definedCount = 0
    for idx in checkableCarbons:
        atom = mol.GetAtomWithIdx(idx)
        tag = atom.GetChiralTag()
        if tag != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            definedCount += 1

    totalCheckable = len(checkableCarbons)

    if totalCheckable == 0:
        status = "NO_CHECKABLE"
    elif definedCount == totalCheckable:
        status = "FULL_3D"
    elif definedCount > 0:
        status = "PARTIAL_3D"
    else:
        status = "FLAT_2D"

    return {
        "status": status,
        "defined": definedCount,
        "total": totalCheckable,
        "ratio": definedCount / totalCheckable if totalCheckable > 0 else 0.0,
    }


def auditMolecule(
    smiles: str,
) -> Dict:
    """审计单个分子的所有糖环手性完整度。
    Audit all sugar rings in a single molecule for chirality completeness.

    Args:
        smiles: canonical SMILES

    Returns:
        分子级审计结果
    """
    result = {
        "status": "ERROR",
        "sugar_count": 0,
        "ring_results": [],
        "worst_status": "FULL_3D",
    }

    if not smiles or smiles in ("nan", "None", "", "NULL"):
        return result

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return result

        units = find_mapped_sugar_units(mol)
        if not units:
            result["status"] = "NO_SUGAR"
            return result

        result["sugar_count"] = len(units)
        statusPriority = {"FLAT_2D": 0, "NO_CHECKABLE": 1,
                          "PARTIAL_3D": 2, "FULL_3D": 3, "NO_CARBON": 4}
        worstPriority = 99

        for unit in units:
            ringAtoms = unit.get("ring_atoms", [])
            ringResult = auditSugarRingChirality(mol, ringAtoms)
            result["ring_results"].append(ringResult)

            priority = statusPriority.get(ringResult["status"], 5)
            if priority < worstPriority:
                worstPriority = priority
                result["worst_status"] = ringResult["status"]

        result["status"] = result["worst_status"]
        return result

    except Exception:
        return result


# =====================================================================
# 2. 关联分析: 手性缺失 vs Hex 退避 (Correlation Analysis)
# =====================================================================

GENERIC_PAT = re.compile(r'\b(Hex|Pen|dHex|HexA|HexN|HexNAc|Non|Oct|Hept)\b')


def hasGenericSugarLabel(sequence: str) -> bool:
    """检查糖序列是否包含泛指标签。
    Check if sugar sequence contains generic fallback labels.
    """
    if not sequence or sequence in ("nan", "None", ""):
        return False
    return bool(GENERIC_PAT.search(sequence))


# =====================================================================
# 3. 主入口 (Main Entry)
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        description="GlycoNP Chirality Completeness Audit")
    parser.add_argument("--input", type=str, default=None,
                        help="Input CSV path")
    parser.add_argument("--limit", type=int, default=0,
                        help="Limit rows for quick test (0=all)")
    parser.add_argument("--output", type=str, default=None,
                        help="Output markdown report path")
    args = parser.parse_args()

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")

    if args.input:
        inputPath = args.input
    else:
        # 自动检测最新报告文件 (Auto-detect latest report file)
        candidates = [
            os.path.join(reportDir, "GlycoNP_Deep_Enriched_v12.csv"),
            os.path.join(reportDir, "GlycoNP_Pipeline_Full.csv"),
            os.path.join(reportDir, "GlycoNP_Deep_Enriched_FINAL.csv"),
        ]
        inputPath = next((p for p in candidates if os.path.exists(p)), None)
        if inputPath is None:
            print("ERROR: No report CSV found. Use --input to specify.")
            sys.exit(1)

    print("=" * 70)
    print("  GlycoNP Chirality Completeness Audit")
    print("  手性完整度审计")
    print("=" * 70)

    # 读取数据 (Load data)
    print(f"  Loading: {inputPath}")
    df = pd.read_csv(inputPath, low_memory=False, dtype=str,
                     encoding="utf-8-sig")
    total = len(df)
    print(f"  Total rows: {total:,}")

    if args.limit > 0:
        df = df.head(args.limit)
        total = len(df)
        print(f"  Limited to: {total:,}")

    # 筛选含糖行 (Filter sugar-containing rows)
    if "contains_sugar" in df.columns:
        hasSugar = df["contains_sugar"].astype(str).str.lower().isin(
            ["true", "1", "yes"])
        sugarDf = df[hasSugar].copy()
    else:
        sugarDf = df.copy()

    sugarCount = len(sugarDf)
    print(f"  Sugar-containing rows: {sugarCount:,}\n")

    # === 审计循环 (Audit loop) ===
    from tqdm import tqdm
    t0 = time.time()

    statusCounter = Counter()
    ringStatusCounter = Counter()
    totalRings = 0
    totalDefinedCenters = 0
    totalCheckableCenters = 0

    # 关联分析计数器 (Correlation counters)
    flatWithGeneric = 0
    flatWithoutGeneric = 0
    fullWithGeneric = 0
    fullWithoutGeneric = 0

    for idx in tqdm(sugarDf.index, desc="  Auditing", ncols=80):
        smilesVal = str(sugarDf.at[idx, "canonical_smiles"]) if pd.notna(
            sugarDf.at[idx, "canonical_smiles"]) else ""

        auditResult = auditMolecule(smilesVal)
        statusCounter[auditResult["status"]] += 1

        # 糖环级统计 (Ring-level stats)
        for ringRes in auditResult.get("ring_results", []):
            ringStatusCounter[ringRes["status"]] += 1
            totalRings += 1
            totalDefinedCenters += ringRes.get("defined", 0)
            totalCheckableCenters += ringRes.get("total", 0)

        # 关联分析 (Correlation)
        seqVal = str(sugarDf.at[idx, "Sugar_Sequence"]) if (
            "Sugar_Sequence" in sugarDf.columns and
            pd.notna(sugarDf.at[idx, "Sugar_Sequence"])
        ) else ""

        isGeneric = hasGenericSugarLabel(seqVal)
        isFlat = auditResult["status"] == "FLAT_2D"

        if isFlat and isGeneric:
            flatWithGeneric += 1
        elif isFlat and not isGeneric:
            flatWithoutGeneric += 1
        elif not isFlat and isGeneric:
            fullWithGeneric += 1
        elif not isFlat and not isGeneric:
            fullWithoutGeneric += 1

    elapsed = time.time() - t0

    # === 输出报告 (Output report) ===
    print("\n" + "=" * 70)
    print("  RESULTS / 审计结果")
    print("=" * 70)

    # 分子级 (Molecule level)
    print("\n  [1] Molecule-Level Chirality Status / 分子级手性状态")
    for status in ["FULL_3D", "PARTIAL_3D", "FLAT_2D", "NO_SUGAR", "ERROR"]:
        count = statusCounter.get(status, 0)
        pct = count / sugarCount * 100 if sugarCount else 0
        bar = "█" * int(pct / 2)
        print(f"    {status:15s} {count:>8,} ({pct:5.1f}%) {bar}")

    # 糖环级 (Ring level)
    print(f"\n  [2] Ring-Level Stats / 糖环级统计")
    print(f"    Total sugar rings:    {totalRings:>8,}")
    print(f"    Defined chiral C:     {totalDefinedCenters:>8,}")
    print(f"    Checkable chiral C:   {totalCheckableCenters:>8,}")
    if totalCheckableCenters > 0:
        overallRatio = totalDefinedCenters / totalCheckableCenters * 100
        print(f"    Overall chirality:    {overallRatio:>7.1f}%")

    for status in ["FULL_3D", "PARTIAL_3D", "FLAT_2D", "NO_CHECKABLE"]:
        count = ringStatusCounter.get(status, 0)
        pct = count / totalRings * 100 if totalRings else 0
        print(f"    Ring {status:15s} {count:>8,} ({pct:5.1f}%)")

    # 关联分析 (Correlation)
    print(f"\n  [3] Correlation: Chirality vs Generic Labels / 手性-泛指标签关联")
    print(f"    FLAT_2D + Generic:     {flatWithGeneric:>8,}")
    print(f"    FLAT_2D + Precise:     {flatWithoutGeneric:>8,}")
    print(f"    Has-Chiral + Generic:  {fullWithGeneric:>8,}")
    print(f"    Has-Chiral + Precise:  {fullWithoutGeneric:>8,}")

    if flatWithGeneric + flatWithoutGeneric > 0:
        flatGenericRate = flatWithGeneric / (flatWithGeneric + flatWithoutGeneric) * 100
        print(f"    → FLAT_2D → Generic rate: {flatGenericRate:.1f}%")
    if fullWithGeneric + fullWithoutGeneric > 0:
        fullGenericRate = fullWithGeneric / (fullWithGeneric + fullWithoutGeneric) * 100
        print(f"    → Has-Chiral → Generic rate: {fullGenericRate:.1f}%")

    # 决策建议 (Decision recommendation)
    flatRate = statusCounter.get("FLAT_2D", 0) / sugarCount * 100 if sugarCount else 0
    print(f"\n  [4] Decision Point")
    if flatRate > 30:
        print(f"    [!] FLAT_2D rate = {flatRate:.1f}% (>30%)")
        print(f"    -> Recommend activating ChEMBL InChIKey index first (Plan C)")
    elif flatRate > 10:
        print(f"    [!] FLAT_2D rate = {flatRate:.1f}% (10-30%)")
        print(f"    -> Recommend parallel: Virtual Deprotection (Step 2) + ChEMBL")
    else:
        print(f"    [OK] FLAT_2D rate = {flatRate:.1f}% (<10%)")
        print(f"    -> Chirality loss not primary bottleneck, proceed to Step 2")

    print(f"\n  Time: {elapsed:.1f}s")
    print("=" * 70)

    # 保存 Markdown 报告 (Save markdown report)
    outputPath = args.output or os.path.join(
        reportDir, "Stereo_Completeness_Audit.md")

    mdLines = [
        "# GlycoNP Chirality Completeness Audit",
        "# 手性完整度审计报告\n",
        f"> Source: `{os.path.basename(inputPath)}`",
        f"> Total sugar-containing rows: **{sugarCount:,}**",
        f"> Total sugar rings: **{totalRings:,}**",
        f"> Audit time: {elapsed:.1f}s\n",
        "## 1. Molecule-Level Status / 分子级手性状态\n",
        "| Status | Count | % |",
        "|:-------|------:|--:|",
    ]
    for status in ["FULL_3D", "PARTIAL_3D", "FLAT_2D", "NO_SUGAR", "ERROR"]:
        count = statusCounter.get(status, 0)
        pct = count / sugarCount * 100 if sugarCount else 0
        mdLines.append(f"| **{status}** | {count:,} | {pct:.1f}% |")

    mdLines.extend([
        "\n## 2. Ring-Level Stats / 糖环级统计\n",
        "| Metric | Value |",
        "|:-------|------:|",
        f"| Total sugar rings | {totalRings:,} |",
        f"| Defined chiral centers | {totalDefinedCenters:,} |",
        f"| Checkable chiral centers | {totalCheckableCenters:,} |",
    ])
    if totalCheckableCenters > 0:
        overallRatio = totalDefinedCenters / totalCheckableCenters * 100
        mdLines.append(f"| **Overall chirality ratio** | **{overallRatio:.1f}%** |")

    mdLines.extend([
        "\n## 3. Correlation / 手性-泛指标签关联\n",
        "| Condition | Count |",
        "|:----------|------:|",
        f"| FLAT_2D + Generic label | {flatWithGeneric:,} |",
        f"| FLAT_2D + Precise label | {flatWithoutGeneric:,} |",
        f"| Has-Chiral + Generic label | {fullWithGeneric:,} |",
        f"| Has-Chiral + Precise label | {fullWithoutGeneric:,} |",
    ])

    mdLines.extend([
        f"\n## 4. Decision / 决策\n",
        f"> **FLAT_2D rate = {flatRate:.1f}%**\n",
    ])

    with open(outputPath, "w", encoding="utf-8") as f:
        f.write("\n".join(mdLines))

    print(f"  Report saved: {outputPath}")


if __name__ == "__main__":
    main()
