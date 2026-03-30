"""
verify_csv_swap.py
==================
验证 V13 CSV 列交换后 Glycan_SMILES 确实含糖环, Aglycon_SMILES 确实是苷元。
Verify post-swap CSV: Glycan_SMILES should contain sugar rings, Aglycon_SMILES should be aglycone.

检查方法 (Verification method):
1. 糖环特征: 含氧杂环 (五/六元) + 多个 -OH
2. 苷元特征: 通常含更多碳骨架, 更少 -OH 密度
3. 抽样 100 行做统计

[TEST DATA ONLY]
"""
import os
import sys
import re

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

REPORT_DIR = os.path.join(os.path.dirname(__file__), "..", "reports")


def hasSugarRing(mol):
    """检测分子是否含有糖环特征 (五/六元含氧饱和环)。
    Detect if molecule has sugar ring features (5/6-membered O-containing saturated ring).
    """
    if mol is None:
        return False
    ringInfo = mol.GetRingInfo()
    for ring in ringInfo.AtomRings():
        if len(ring) not in (5, 6):
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        symbols = [a.GetSymbol() for a in atoms]
        # 糖环: 恰好 1 个 O + 其余为 C (Sugar ring: exactly 1 O + rest C)
        if symbols.count("O") == 1 and symbols.count("C") == len(ring) - 1:
            return True
    return False


def analyzeSmiles(smi):
    """分析 SMILES 的化学特征。
    Analyze chemical features of a SMILES string.
    """
    if pd.isna(smi) or not smi or smi == "nan":
        return {"valid": False}
    mol = Chem.MolFromSmiles(str(smi))
    if mol is None:
        return {"valid": False}

    heavyAtomCount = mol.GetNumHeavyAtoms()
    oCount = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "O")
    cCount = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "C")
    nCount = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "N")
    hasSugar = hasSugarRing(mol)
    oRatio = oCount / heavyAtomCount if heavyAtomCount > 0 else 0

    return {
        "valid": True,
        "heavy": heavyAtomCount,
        "C": cCount,
        "O": oCount,
        "N": nCount,
        "O_ratio": oRatio,
        "has_sugar_ring": hasSugar,
    }


def verifyFile(csvPath, name):
    """验证单个 CSV 文件。"""
    print(f"\n{'='*70}")
    print(f"验证 {name}: {csvPath}")
    print(f"{'='*70}")

    if not os.path.exists(csvPath):
        print(f"  [SKIP] 文件不存在")
        return

    df = pd.read_csv(csvPath, low_memory=False, nrows=500)
    print(f"  读取前 {len(df)} 行")

    if "Glycan_SMILES" not in df.columns or "Aglycon_SMILES" not in df.columns:
        print(f"  [SKIP] 缺少 Glycan_SMILES 或 Aglycon_SMILES 列")
        return

    # 筛选两列都非空的行 (Filter rows where both columns are non-empty)
    mask = df["Glycan_SMILES"].notna() & df["Aglycon_SMILES"].notna()
    mask = mask & (df["Glycan_SMILES"] != "") & (df["Aglycon_SMILES"] != "")
    mask = mask & (df["Glycan_SMILES"] != "nan") & (df["Aglycon_SMILES"] != "nan")
    validDf = df[mask].head(200)
    print(f"  有效双列行数: {len(validDf)}")

    glycanSugarCount = 0
    aglyconSugarCount = 0
    glycanHighO = 0
    aglyconHighO = 0
    sampleResults = []

    for idx, row in validDf.iterrows():
        glyInfo = analyzeSmiles(row["Glycan_SMILES"])
        aglInfo = analyzeSmiles(row["Aglycon_SMILES"])

        if not glyInfo["valid"] or not aglInfo["valid"]:
            continue

        if glyInfo["has_sugar_ring"]:
            glycanSugarCount += 1
        if aglInfo["has_sugar_ring"]:
            aglyconSugarCount += 1
        if glyInfo["O_ratio"] > 0.3:
            glycanHighO += 1
        if aglInfo["O_ratio"] > 0.3:
            aglyconHighO += 1

        sampleResults.append({
            "idx": idx,
            "glycan_sugar": glyInfo["has_sugar_ring"],
            "glycan_O_ratio": glyInfo["O_ratio"],
            "glycan_heavy": glyInfo["heavy"],
            "aglycon_sugar": aglInfo["has_sugar_ring"],
            "aglycon_O_ratio": aglInfo["O_ratio"],
            "aglycon_heavy": aglInfo["heavy"],
        })

    total = len(sampleResults)
    if total == 0:
        print("  [WARN] 无有效样本")
        return

    print(f"\n  --- 统计 ({total} 有效样本) ---")
    print(f"  Glycan_SMILES 含糖环: {glycanSugarCount}/{total} ({glycanSugarCount/total*100:.1f}%)")
    print(f"  Aglycon_SMILES 含糖环: {aglyconSugarCount}/{total} ({aglyconSugarCount/total*100:.1f}%)")
    print(f"  Glycan_SMILES 高氧比(>0.3): {glycanHighO}/{total} ({glycanHighO/total*100:.1f}%)")
    print(f"  Aglycon_SMILES 高氧比(>0.3): {aglyconHighO}/{total} ({aglyconHighO/total*100:.1f}%)")

    # 判定 (Verdict)
    if glycanSugarCount > aglyconSugarCount:
        print(f"\n  ✅ [CORRECT] Glycan 列含糖环比例 ({glycanSugarCount}) > Aglycon 列 ({aglyconSugarCount})")
        print(f"     列交换正确！")
    elif glycanSugarCount < aglyconSugarCount:
        print(f"\n  ❌ [WRONG] Aglycon 列含糖环比例 ({aglyconSugarCount}) > Glycan 列 ({glycanSugarCount})")
        print(f"     列可能仍然是反的！")
    else:
        print(f"\n  ⚠️ [INCONCLUSIVE] 两列糖环比例相同 ({glycanSugarCount} vs {aglyconSugarCount})")

    # 显示前 5 个样本详情 (Show first 5 sample details)
    print(f"\n  --- 前 5 个样本详情 ---")
    for i, sr in enumerate(sampleResults[:5]):
        print(f"  Row {sr['idx']}: "
              f"Glycan(sugar={sr['glycan_sugar']}, O%={sr['glycan_O_ratio']:.2f}, atoms={sr['glycan_heavy']}) | "
              f"Aglycon(sugar={sr['aglycon_sugar']}, O%={sr['aglycon_O_ratio']:.2f}, atoms={sr['aglycon_heavy']})")


# ─── 主程序 ───
verifyFile(
    os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched_v13_Pruned.csv"),
    "V13 Pruned"
)
verifyFile(
    os.path.join(REPORT_DIR, "GlycoNP_Saponin_DB_v13.csv"),
    "Saponin DB"
)
