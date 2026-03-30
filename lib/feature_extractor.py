"""
Phase 5: 骨架提取与糖序列 — 统一特征化入口
Phase 5: Scaffold Extraction & Sugar Sequence — Unified Feature Entry Point

苷元特征化 (Aglycone Features):
  - Morgan Fingerprint (radius=2, 2048-bit)
  - Bemis-Murcko Scaffold (含脂肪链退避标记)

糖基特征化 (Glycan Features):
  - 伪 IUPAC 糖序列（三层退避匹配）
    Tier 1: 精确立体化学 SMILES 库匹配 (α-D-Glc, β-D-Gal, etc.)
    Tier 2: CIP 构型图遍历救援 (D-Gal(Rescued), D-Glc(Rescued))
    Tier 3: 骨架级别退避 (Hex, Pen, dHex)
"""
import os
import sys
from typing import Optional, Tuple, Dict, List

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__))))

from monosaccharide_identifier import analyze_glycan


# =====================================================================
# 1. 苷元特征化 (Aglycone Feature Extraction)
# =====================================================================

def computeMorganFingerprint(
    smiles: str,
    radius: int = 2,
    nBits: int = 2048,
) -> Optional[str]:
    """
    从 SMILES 计算 Morgan 指纹 (ECFP-like) 并返回 bit 字符串。
    Compute Morgan fingerprint from SMILES and return as bit string.

    Args:
        smiles: Aglycan 或任意分子的 SMILES
        radius: 指纹半径 (default=2, 等同于 ECFP4)
        nBits: 位数 (default=2048)

    Returns:
        2048 位 bit 字符串, 或 None (如果 SMILES 无效)
    """
    if not smiles or smiles in ("NULL", "nan", "", "*"):
        return None
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return None
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
        return fp.ToBitString()
    except Exception:
        return None


def extractMurckoScaffold(smiles: str) -> str:
    """
    提取 Bemis-Murcko 骨架 SMILES。
    Extract Bemis-Murcko scaffold SMILES.

    退避策略 (Fallback Strategy):
    - 如果骨架为空（直链脂肪酸/简单链状分子），返回 "Aliphatic Chain"
    - 如果 SMILES 无效，返回 "NULL"

    Args:
        smiles: Aglycan SMILES

    Returns:
        Murcko Scaffold SMILES, "Aliphatic Chain", or "NULL"
    """
    if not smiles or smiles in ("NULL", "nan", "", "*"):
        return "NULL"
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return "NULL"

        # 提取 Murcko 骨架 (Extract Murcko scaffold)
        scaffold = GetScaffoldForMol(mol)
        scaffoldSmiles = Chem.MolToSmiles(scaffold)

        # 退避: 空骨架 → 脂肪链标记 (Fallback: empty scaffold → aliphatic chain)
        if not scaffoldSmiles or scaffoldSmiles == "":
            # 检查是否为脂肪链: 无环且碳数 ≥ 4
            # Check if aliphatic chain: no rings and carbon count ≥ 4
            ri = mol.GetRingInfo()
            carbonCount = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
            if ri.NumRings() == 0 and carbonCount >= 4:
                return "Aliphatic Chain"
            else:
                return "NULL"

        return scaffoldSmiles

    except Exception:
        return "NULL"


def getTopologyScaffoldSmiles(molOrSmiles) -> str:
    """
    全碳简化拓扑骨架提取 (Topology-only Murcko Scaffold).

    四步转换 (Four-step transformation):
      1. 提取 Bemis-Murcko 骨架 — 去除所有侧链
         Extract Murcko scaffold — remove all side chains
      2. 杂原子 (N, O, S, P 等) → C — 全碳化
         Replace ALL non-carbon heavy atoms with carbon
      3. 非芳香环内双键 → 单键 — 去除局部构象干扰
         Saturate non-aromatic ring double bonds (e.g. db5/db10 in solanidine)
      4. 去除手性/键立体 — 纯拓扑连接
         Strip chirality and bond stereo for topology-only representation

    设计意图 (Design Rationale):
      为下游《数据可视化与图谱库》中的热力图/聚类分析提供纯碳拓扑指纹。
      消除杂原子和局部构象差异, 仅保留环-桥-链的连接拓扑。
      Provides a pure-carbon topology fingerprint for downstream heatmap/clustering
      analysis, removing heteroatom and conformational noise to retain only the
      ring-bridge-chain connectivity topology.

    Args:
        molOrSmiles: RDKit Mol 对象或 SMILES 字符串

    Returns:
        全碳拓扑骨架 SMILES (non-chiral, canonical), "NULL", 或 "Aliphatic Chain"
    """
    # 输入标准化: 同时接受 Mol 和 SMILES (Input normalization)
    if isinstance(molOrSmiles, str):
        if not molOrSmiles or molOrSmiles in ("NULL", "nan", "", "*"):
            return "NULL"
        mol = Chem.MolFromSmiles(molOrSmiles)
    else:
        mol = molOrSmiles

    if mol is None:
        return "NULL"

    try:
        from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol

        # Step 1: 提取 Murcko 骨架 — 去除所有侧链
        # Extract Murcko scaffold — strips all side chains, keeps ring system + linkers
        scaffold = GetScaffoldForMol(mol)
        if scaffold is None or scaffold.GetNumAtoms() == 0:
            ri = mol.GetRingInfo()
            carbonCount = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
            if ri.NumRings() == 0 and carbonCount >= 4:
                return "Aliphatic Chain"
            return "NULL"

        # Step 2: 杂原子全替换 — 所有非碳重原子 → C(6)
        # Replace ALL non-carbon heavy atoms (N=7, O=8, S=16, P=15, etc.) with C
        rw = Chem.RWMol(scaffold)
        for atom in rw.GetAtoms():
            atomicNum = atom.GetAtomicNum()
            if atomicNum != 6 and atomicNum != 1:  # 非碳且非氢
                atom.SetAtomicNum(6)
                atom.SetFormalCharge(0)
                atom.SetNumExplicitHs(0)
                atom.SetNoImplicit(False)
            # Step 4a: 去除原子手性 (Strip atom chirality)
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

        # Step 3: 非芳香环内双键/三键 → 单键
        # Saturate non-aromatic double/triple bonds to single bonds
        for bond in rw.GetBonds():
            # Step 4b: 去除键立体 (Strip bond stereo)
            bond.SetStereo(Chem.BondStereo.STEREONONE)
            if not bond.GetIsAromatic():
                if bond.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE):
                    bond.SetBondType(Chem.BondType.SINGLE)

        try:
            Chem.SanitizeMol(rw)
        except Exception:
            pass  # 尽力而为 (best-effort sanitization)

        # Step 4c: 输出非手性规范 SMILES (Non-chiral canonical SMILES)
        topoSmiles = Chem.MolToSmiles(rw, isomericSmiles=False)
        return topoSmiles if topoSmiles else "NULL"

    except Exception:
        return "NULL"


def calculateHierarchicalScore(
    molOrSmiles1,
    molOrSmiles2,
) -> Dict[str, float]:
    """
    分层骨架相似性评分 (Hierarchical Scaffold Similarity Score).

    两层评分策略 (Two-tier scoring):
      Part 1 — 精确骨架匹配 (Exact Scaffold Match):
        全碳拓扑骨架 SMILES 完全一致 → 100, 否则 → 0
      Part 2 — 相似度评分 (Likeness Score, 仅在 Part 1 = 100 时计算):
        原始分子的 Morgan 指纹 Tanimoto 相似度 × 100

    设计意图 (Design Rationale):
      先用全碳拓扑骨架做"硬筛选"(同骨架 vs 异骨架),
      再用分子指纹做"软排序"(同骨架内的修饰差异排序)。
      First hard-filter by topology scaffold identity, then soft-rank by
      Morgan fingerprint similarity for molecules sharing the same scaffold.

    Args:
        molOrSmiles1: 分子1 (RDKit Mol 或 SMILES)
        molOrSmiles2: 分子2 (RDKit Mol 或 SMILES)

    Returns:
        {"scaffold_match_score": 0|100, "likeness_score": 0.0~100.0}
    """
    result = {"scaffold_match_score": 0, "likeness_score": 0.0}

    # 输入标准化 (Input normalization)
    def _toMol(x):
        if isinstance(x, str):
            return Chem.MolFromSmiles(x) if x and x not in ("NULL", "nan") else None
        return x

    mol1 = _toMol(molOrSmiles1)
    mol2 = _toMol(molOrSmiles2)
    if mol1 is None or mol2 is None:
        return result

    # Part 1: 精确骨架匹配 (Exact Scaffold Match)
    topo1 = getTopologyScaffoldSmiles(mol1)
    topo2 = getTopologyScaffoldSmiles(mol2)

    if topo1 == "NULL" or topo2 == "NULL":
        return result

    if topo1 == topo2:
        result["scaffold_match_score"] = 100
    else:
        # 骨架不同 → 两个分数都是 0 (Different scaffold → both scores 0)
        return result

    # Part 2: 相似度评分 (Likeness Score via Morgan Tanimoto)
    # 仅当骨架完全匹配时才进入此步骤
    # Only computed when scaffolds match exactly
    try:
        from rdkit import DataStructs
        fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius=2, nBits=2048)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius=2, nBits=2048)
        tanimoto = DataStructs.TanimotoSimilarity(fp1, fp2)
        result["likeness_score"] = round(tanimoto * 100.0, 2)
    except Exception:
        result["likeness_score"] = 0.0

    return result



def characterizeAglycon(smiles: str) -> Dict[str, str]:
    """
    苷元一站式特征化: Morgan FP + Murcko Scaffold + 拓扑骨架 + 基本描述符。
    One-stop aglycone characterization: Morgan FP + Murcko + Topology Scaffold + descriptors.

    Args:
        smiles: Aglycon SMILES

    Returns:
        字典包含: morgan_fp, murcko_scaffold, topology_scaffold, molecular_weight, ring_count
    """
    result = {
        "morgan_fp": "",
        "murcko_scaffold": "NULL",
        "topology_scaffold": "NULL",
        "aglycon_mw": "",
        "aglycon_ring_count": "",
    }

    if not smiles or smiles in ("NULL", "nan", "", "*"):
        return result

    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return result

        # Morgan 指纹 (Morgan fingerprint)
        fp = computeMorganFingerprint(smiles)
        if fp:
            result["morgan_fp"] = fp

        # Murcko 骨架 (Murcko scaffold)
        result["murcko_scaffold"] = extractMurckoScaffold(smiles)

        # 全碳拓扑骨架 (Topology-only scaffold)
        result["topology_scaffold"] = getTopologyScaffoldSmiles(mol)

        # 分子量 (Molecular weight)
        result["aglycon_mw"] = f"{Descriptors.ExactMolWt(mol):.2f}"

        # 环数 (Ring count)
        result["aglycon_ring_count"] = str(mol.GetRingInfo().NumRings())

    except Exception:
        pass

    return result


# =====================================================================
# 2. 糖基特征化 (Glycan Feature Extraction)
# =====================================================================

def characterizeGlycan(smiles: str) -> Dict[str, str]:
    """
    糖基一站式特征化: 三层退避的伪 IUPAC 序列。
    One-stop glycan characterization: tiered fallback pseudo-IUPAC sequence.

    匹配层级 (Matching Tiers):
      Tier 1: 精确 SMILES 库 → α-D-Glc, β-D-Gal 等
      Tier 2: CIP 图遍历救援 → D-Gal(Rescued) 等
      Tier 3: 骨架退避 → Hex, Pen, dHex

    Args:
        smiles: Glycan SMILES (可含多个糖环)

    Returns:
        字典包含:
          sugar_sequence: 人类可读的伪 IUPAC 序列 (e.g. "D-Glc-(b1→?)-L-Rha")
          sugar_functional_group: 详细修饰标注 (e.g. "D-Glc_1() ; L-Rha_2(*deoxy)")
    """
    result = {
        "sugar_sequence": "",
        "sugar_functional_group": "",
    }

    if not smiles or smiles in ("NULL", "nan", "", "*"):
        return result

    try:
        seq, mods = analyze_glycan(smiles)
        result["sugar_sequence"] = seq if seq else ""
        result["sugar_functional_group"] = mods if mods else ""
    except Exception:
        result["sugar_sequence"] = "Error"

    return result


# =====================================================================
# 3. 统一 Phase 5 入口 (Unified Phase 5 Entry Point)
# =====================================================================

def processPhase5Row(
    glycanSmiles: str,
    aglycanSmiles: str,
) -> Dict[str, str]:
    """
    对一行数据执行完整 Phase 5 特征化。
    Execute full Phase 5 characterization for one row of data.

    Args:
        glycanSmiles: Phase 2 拆分出的 Glycan SMILES
        aglycanSmiles: Phase 2 拆分出的 Aglycan SMILES

    Returns:
        字典包含所有 Phase 5 产出列
    """
    glycanFeatures = characterizeGlycan(glycanSmiles)
    aglycanFeatures = characterizeAglycon(aglycanSmiles)

    return {**glycanFeatures, **aglycanFeatures}


# =====================================================================
# 4. 批量处理 (Batch Processing)
# =====================================================================

def processPhase5Batch(
    df: pd.DataFrame,
    glycanCol: str = "Glycan_SMILES",
    aglycanCol: str = "Aglycan_SMILES",
) -> pd.DataFrame:
    """
    对 DataFrame 批量执行 Phase 5，将特征列直接追加到 DataFrame。
    Batch process Phase 5 for a DataFrame, appending feature columns.

    Args:
        df: 输入 DataFrame (必须包含 Glycan_SMILES 和 Aglycan_SMILES 列)
        glycanCol: 糖基 SMILES 列名
        aglycanCol: 苷元 SMILES 列名

    Returns:
        添加了 Phase 5 特征列的 DataFrame
    """
    from tqdm import tqdm
    tqdm.pandas(desc="Phase 5 Feature Extraction")

    results = df.progress_apply(
        lambda row: processPhase5Row(
            str(row.get(glycanCol, "")),
            str(row.get(aglycanCol, "")),
        ),
        axis=1,
        result_type="expand",
    )

    return pd.concat([df, results], axis=1)


if __name__ == "__main__":
    # 快速自测 (Quick self-test)
    print("=" * 70)
    print("Phase 5 Quick Self-Test")
    print("=" * 70)

    # Rutin (芦丁): quercetin-3-O-rutinoside
    # 含一个葡萄糖 (Glc) + 一个鼠李糖 (Rha), 通过 1→6 连接
    rutinSmiles = "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@H]2[C@@H](O)[C@H](O)[C@@H](O)[C@@H](C)O2)C(CO)O1"
    aglycanSmiles = "O=c1cc(-c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"

    print(f"\nInput Glycan SMILES: {rutinSmiles}")
    print(f"Input Aglycan SMILES: {aglycanSmiles}")

    result = processPhase5Row(rutinSmiles, aglycanSmiles)

    print(f"\n--- Glycan Features ---")
    print(f"  Sugar Sequence:        {result['sugar_sequence']}")
    print(f"  Sugar Functional Group: {result['sugar_functional_group']}")

    print(f"\n--- Aglycon Features ---")
    print(f"  Murcko Scaffold:       {result['murcko_scaffold']}")
    print(f"  Aglycon MW:            {result['aglycon_mw']}")
    print(f"  Aglycon Ring Count:    {result['aglycon_ring_count']}")
    print(f"  Morgan FP (first 50):  {result['morgan_fp'][:50]}...")
