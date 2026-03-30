"""
Phase 3: 二级片段扫描器 — 核苷酸糖 & 糖肽识别
Phase 3: Secondary Fragment Scanner — Nucleotide Sugars & Glycopeptides

核苷酸糖 (Nucleotide Sugars):
  识别含有嘌呤/嘧啶碱基 + 核糖 + 磷酸基团的糖缀合物。
  典型例子: UDP-Glc, GDP-Man, CMP-Neu5Ac

糖肽/氨基酸糖苷 (Glycopeptides / Amino Acid Glycosides):
  识别含有肽键或 α-氨基酸骨架的糖缀合物。
  关键: 区分真正的肽键与糖上的 N-乙酰基 (NAc)。

设计意图 (Design Rationale):
  - SMARTS 模式经过精细调优，避免误匹配
  - NAc 排除使用 "邻碳环境限制" 而非简单黑名单
  - 核苷酸检测要求 碱基 + 磷酸 共存 (两者缺一则不算)
"""
import os
import sys
from typing import Dict, List, Optional, Tuple
from collections import OrderedDict

from rdkit import Chem

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__))))


# =====================================================================
# 1. SMARTS 定义 (Pattern Definitions)
# =====================================================================

# --- 核苷酸碱基 SMARTS (Nucleobase SMARTS) ---
# 设计: 每种碱基用独立的高特异性 SMARTS
NUCLEOBASE_SMARTS: Dict[str, str] = {
    # 嘌呤骨架 (Purine) — 腺嘌呤, 鸟嘌呤, 次黄嘌呤
    # 使用小写芳香 SMILES, RDKit 在 N-糖苷键场景中将嘌呤视为芳香
    "Purine":           "n1cnc2c1ncnc2",

    # 嘧啶二酮 (Pyrimidine dione) — 尿嘧啶, 胸腺嘧啶
    # 使用 ~ (any bond) 同时匹配芳香态和 Kekulé 态
    "Uracil_Thymine":   "[#8]=[#6]1~[#7]~[#6](=[#8])~[#7]~[#6]~[#6]1",

    # 胞嘧啶 (Cytosine) — 4-氨基嘧啶-2-酮
    # 同样用 ~ 兼容芳香/非芳香
    "Cytosine":         "[#7]~[#6]1~[#6]~[#6]~[#7]~[#6](=[#8])~[#7]~1",
}

# --- 磷酸基团 SMARTS (Phosphate) ---
PHOSPHATE_SMARTS = "P(=O)([O,OH,$([O-])])([O,OH,$([O-])])"

# --- 肽键 SMARTS (Peptide Bond) ---
# 关键防误判设计:
#   真正的肽键: N-C(=O)-Cα-Cβ, 其中 Cα 至少连接一个 H 和一个侧链碳/杂原子
#   NAc 伪肽键: N-C(=O)-CH3, 其中 C(=O) 后面接的是纯甲基
#
# 策略: 要求酰胺碳连接的 α-碳 (Cα) 还必须连接至少一个非 H 的原子
# (NAc 的甲基碳 CH3 只有 3个H, 不满足 Cα 的侧链条件)
PEPTIDE_BOND_SMARTS: Dict[str, str] = {
    # 标准肽键: N-C(=O)-Cα, 且 Cα 不是 CH3 (排除 NAc)
    # [CX4;!CH3] 要求 α-碳至少有一个非氢取代基 (即不是甲基)
    "Peptide_Bond":     "[NX3;!$(NC=S)]-[CX3](=[OX1])-[CX4;!CH3]",

    # α-氨基酸骨架: NH2-Cα(H)(R)-COOH 或其衍生物
    # 这匹配游离氨基酸或氨基酸糖苷 (如 Ser-糖苷)
    "Alpha_AminoAcid":  "[NX3;H2,H1]-[CX4;H1](-[*])-[CX3](=[OX1])-[OX2,OX1,$([O-])]",
}

# --- 额外特征 (Aromatic amino acid residues) ---
AROMATIC_RESIDUE_SMARTS = "c1ccccc1-[CH2]-[CX4](-[NX3])-[CX3](=[OX1])"


# =====================================================================
# 2. 预编译 (Pre-compile All Patterns)
# =====================================================================

_COMPILED: Dict[str, Optional[Chem.Mol]] = {}


def _ensureCompiled() -> Dict[str, Optional[Chem.Mol]]:
    """懒加载编译所有 SMARTS (Lazy compile all patterns)."""
    global _COMPILED
    if _COMPILED:
        return _COMPILED

    allSmarts = {}
    for name, sma in NUCLEOBASE_SMARTS.items():
        allSmarts[f"NB_{name}"] = sma
    allSmarts["Phosphate"] = PHOSPHATE_SMARTS
    for name, sma in PEPTIDE_BOND_SMARTS.items():
        allSmarts[f"PEP_{name}"] = sma
    allSmarts["Aromatic_Residue"] = AROMATIC_RESIDUE_SMARTS

    for name, sma in allSmarts.items():
        mol = Chem.MolFromSmarts(sma)
        if mol is None:
            print(f"  [WARNING] Phase 3: Failed to compile SMARTS '{name}': {sma}")
        _COMPILED[name] = mol

    return _COMPILED


# =====================================================================
# 3. 核苷酸糖检测 (Nucleotide Sugar Detection)
# =====================================================================

def detectNucleotideSugar(smiles: str) -> Tuple[bool, str]:
    """
    检测分子是否含有核苷酸糖特征:
    必须同时满足: (1) 碱基存在 + (2) 磷酸基团存在。
    Detect if molecule contains nucleotide sugar features.
    Requires BOTH: nucleobase AND phosphate group.

    Args:
        smiles: 全分子 SMILES

    Returns:
        (is_nucleotide_sugar, detail_string)
        detail: e.g. "Uracil_Thymine+Phosphate" or ""
    """
    if not smiles or smiles in ("NULL", "nan", "", "*"):
        return False, ""

    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return False, ""

    patterns = _ensureCompiled()

    # Step 1: 检测碱基 (Detect nucleobase)
    foundBases = []
    for key in ["NB_Purine", "NB_Uracil_Thymine", "NB_Cytosine"]:
        pat = patterns.get(key)
        if pat and mol.HasSubstructMatch(pat):
            baseName = key.replace("NB_", "")
            foundBases.append(baseName)

    # Step 2: 检测磷酸 (Detect phosphate)
    hasPhos = False
    phPat = patterns.get("Phosphate")
    if phPat and mol.HasSubstructMatch(phPat):
        hasPhos = True

    # 核苷酸糖判定: 碱基 + 磷酸必须共存
    # Nucleotide sugar: BOTH base AND phosphate required
    if foundBases and hasPhos:
        detail = "+".join(foundBases) + "+Phosphate"
        return True, detail

    return False, ""


# =====================================================================
# 4. 糖肽/氨基酸检测 (Glycopeptide / Amino Acid Detection)
# =====================================================================

def detectPeptideOrAminoAcid(smiles: str) -> Tuple[bool, str]:
    """
    检测分子是否含有肽键或氨基酸骨架。
    内置 NAc 防误判: 排除 N-C(=O)-CH3 模式。
    Detect peptide bonds or amino acid skeletons.
    Built-in NAc anti-false-positive: excludes N-C(=O)-CH3.

    Args:
        smiles: 全分子 SMILES (或 Aglycan SMILES)

    Returns:
        (has_peptide_feature, detail_string)
    """
    if not smiles or smiles in ("NULL", "nan", "", "*"):
        return False, ""

    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return False, ""

    patterns = _ensureCompiled()
    foundFeatures = []

    # 肽键检测 (Peptide bond)
    pepPat = patterns.get("PEP_Peptide_Bond")
    if pepPat:
        matches = mol.GetSubstructMatches(pepPat)
        if matches:
            foundFeatures.append(f"Peptide_Bond(x{len(matches)})")

    # α-氨基酸骨架 (Alpha amino acid skeleton)
    aaPat = patterns.get("PEP_Alpha_AminoAcid")
    if aaPat:
        matches = mol.GetSubstructMatches(aaPat)
        if matches:
            foundFeatures.append(f"Alpha_AminoAcid(x{len(matches)})")

    # 芳香族残基 (Aromatic residue — Phe/Tyr/Trp)
    arPat = patterns.get("Aromatic_Residue")
    if arPat and mol.HasSubstructMatch(arPat):
        foundFeatures.append("Aromatic_Residue")

    if foundFeatures:
        return True, ", ".join(foundFeatures)
    return False, ""


# =====================================================================
# 5. 统一扫描接口 (Unified Scan Interface)
# =====================================================================

def scanSecondaryFragments(smiles: str) -> Dict[str, object]:
    """
    统一二级片段扫描: 核苷酸 + 肽键。
    Unified secondary fragment scan: nucleotides + peptides.

    Args:
        smiles: 全分子 SMILES

    Returns:
        {
            "Has_Nucleotide": bool,
            "Nucleotide_Detail": str,
            "Has_Peptide": bool,
            "Peptide_Detail": str,
        }
    """
    isNuc, nucDetail = detectNucleotideSugar(smiles)
    isPep, pepDetail = detectPeptideOrAminoAcid(smiles)

    return {
        "Has_Nucleotide": isNuc,
        "Nucleotide_Detail": nucDetail,
        "Has_Peptide": isPep,
        "Peptide_Detail": pepDetail,
    }


# =====================================================================
# 6. 批量处理 (Batch Processing)
# =====================================================================

def batchScanSecondary(
    df: "pd.DataFrame",
    smilesCol: str = "canonical_smiles",
) -> "pd.DataFrame":
    """
    对 DataFrame 批量扫描二级片段，追加 4 列。
    Batch scan secondary fragments, appending 4 new columns.

    新增列 (New Columns):
      - Has_Nucleotide (bool)
      - Nucleotide_Detail (str)
      - Has_Peptide (bool)
      - Peptide_Detail (str)
    """
    from tqdm import tqdm

    nucFlags, nucDetails, pepFlags, pepDetails = [], [], [], []
    for smiles in tqdm(df[smilesCol].fillna(""), desc="Phase 3: Secondary Scan"):
        r = scanSecondaryFragments(str(smiles))
        nucFlags.append(r["Has_Nucleotide"])
        nucDetails.append(r["Nucleotide_Detail"])
        pepFlags.append(r["Has_Peptide"])
        pepDetails.append(r["Peptide_Detail"])

    df["Has_Nucleotide"] = nucFlags
    df["Nucleotide_Detail"] = nucDetails
    df["Has_Peptide"] = pepFlags
    df["Peptide_Detail"] = pepDetails

    return df
