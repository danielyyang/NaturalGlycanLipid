"""
Phase 6: 智能化学分类引擎
Phase 6: Intelligent Chemical Classification Engine

三层分类策略 (Three-tier Classification Strategy):
  Tier 1: SMARTS 骨架匹配 — 甾体/黄酮/香豆素/生物碱/萜类 等
  Tier 2: 糖脂捕获 — 长链脂肪烃(>C10)/神经酰胺/甘油骨架
  Tier 3: Tanimoto 最近邻 — Morgan FP 相似度 ≥ 0.85 继承分类

设计原则 (Design Principles):
  - Tier 1 命中即停止，不进入后续层级
  - Tier 2 独立于 Tier 1 运行（糖脂可能同时含有环状骨架）
  - Tier 3 作为最后的救援机制
"""
import os
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__))))


# =====================================================================
# 1. SMARTS 骨架字典 (Core Skeleton SMARTS Dictionary)
# =====================================================================

# 使用渐进细化策略: 先匹配最特异的子结构，再匹配宽泛骨架
# Progressive refinement: match most specific substructures first, then broader skeletons
CLASSIFICATION_SMARTS = {
    # ----------------------------------------------------------------
    # 甾体与三萜 (Steroids & Triterpenoids)
    # ----------------------------------------------------------------
    # 甾体核: 环戊烷并全氢化菲 (Cyclopentanoperhydrophenanthrene)
    # 四环稠合 A/B/C/D rings — 广义匹配 (Generic steroid nucleus)
    "Steroids": [
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~[#6]~1",
    ],
    # 三萜: 五环三萜骨架 (Pentacyclic triterpenoid, e.g. oleanane/ursane)
    "Triterpenoids": [
        # 齐墩果烷型 (Oleanane-type): 5 fused rings, ≥25 carbons
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~5~[#6]~[#6]~[#6]~[#6]~[#6]~5~[#6]~4~[#6]~3~[#6]~2~[#6]~1",
    ],

    # ----------------------------------------------------------------
    # 黄酮类 (Flavonoids & Isoflavonoids)
    # ----------------------------------------------------------------
    # 黄酮核心: 2-苯基色原酮 (2-Phenylchromanone / chromen-4-one)
    "Flavonoids": [
        # flavone / flavanone 核心 C6-C3-C6
        "c1ccc(-c2cc(=O)c3ccccc3o2)cc1",    # flavone
        "c1ccc(-c2oc3ccccc3c(=O)c2)cc1",     # isoflavone
        "c1ccc(C2CC(=O)c3ccccc3O2)cc1",      # flavanone
        "c1ccc(-c2cc(=O)c3cc(O)cc(O)c3o2)cc1",  # flavone with OH
        "O=c1cc(-c2ccccc2)oc2ccccc12",       # generic chromen-4-one
        "O=C1CC(-c2ccccc2)Oc2ccccc21",       # flavanone relaxed
    ],

    # ----------------------------------------------------------------
    # 香豆素 (Coumarins)
    # ----------------------------------------------------------------
    # 苯并吡喃-2-酮 (2H-chromen-2-one)
    "Coumarins": [
        "O=c1ccc2ccccc2o1",   # coumarin core
        "O=c1ccc2cc(O)ccc2o1",  # hydroxy-coumarin
        "O=c1oc2ccccc2cc1",   # iso-coumarin
    ],

    # ----------------------------------------------------------------
    # 生物碱 (Alkaloids)
    # ----------------------------------------------------------------
    # 含氮杂环 (Nitrogen heterocycles in ring systems)
    "Alkaloids": [
        # 吲哚生物碱 (Indole alkaloids)
        "c1ccc2[nH]ccc2c1",
        # 异喹啉 (Isoquinoline)
        "c1ccc2ncccc2c1",
        # 喹啉 (Quinoline)
        "c1ccc2cnccc2c1",
        # 吡啶 (Pyridine) — 宽泛匹配
        "c1ccncc1",
        # 吡咯烷 (Pyrrolidine ring with N)
        "C1CCNC1",
        # 哌啶 (Piperidine)
        "C1CCNCC1",
    ],

    # ----------------------------------------------------------------
    # 蒽醌 (Anthraquinones)
    # ----------------------------------------------------------------
    "Anthraquinones": [
        "O=C1c2ccccc2C(=O)c2ccccc21",
    ],

    # ----------------------------------------------------------------
    # 木脂素 (Lignans)
    # ----------------------------------------------------------------
    "Lignans": [
        # 二苄基丁内酯 (Dibenzylbutyrolactone) — 常见木脂素
        "c1ccc(CC2COC(=O)C2Cc2ccccc2)cc1",
    ],

    # ----------------------------------------------------------------
    # 苯丙素 (Phenylpropanoids) — 宽泛
    # ----------------------------------------------------------------
    "Phenylpropanoids": [
        "c1ccc(/C=C/C=O)cc1",    # cinnamaldehyde
        "c1ccc(/C=C/C(=O)O)cc1", # cinnamic acid
    ],
}

# 预编译 SMARTS (Pre-compile SMARTS patterns)
_COMPILED_SMARTS: Dict[str, List[Chem.Mol]] = {}


def _getCompiledSmarts() -> Dict[str, List[Chem.Mol]]:
    """懒加载预编译 SMARTS 模式 (Lazy-load pre-compiled SMARTS patterns)."""
    global _COMPILED_SMARTS
    if not _COMPILED_SMARTS:
        for className, patterns in CLASSIFICATION_SMARTS.items():
            compiled = []
            for p in patterns:
                mol = Chem.MolFromSmarts(p)
                if mol is not None:
                    compiled.append(mol)
            _COMPILED_SMARTS[className] = compiled
    return _COMPILED_SMARTS


# =====================================================================
# 2. Tier 1: SMARTS 骨架匹配 (SMARTS Skeleton Matching)
# =====================================================================

def classifyBySMARTS(mol: Chem.Mol) -> Optional[str]:
    """
    使用 SMARTS 子结构匹配对苷元进行分类。
    Classify aglycone using SMARTS substructure matching.

    匹配顺序按特异性递减排列 (Matching order: most specific first).
    第一个命中即返回 (Returns on first match).

    Args:
        mol: RDKit Mol 对象

    Returns:
        分类名称或 None (如果未命中)
    """
    smartsDict = _getCompiledSmarts()
    for className, patterns in smartsDict.items():
        for pattern in patterns:
            if mol.HasSubstructMatch(pattern):
                return className
    return None


# =====================================================================
# 3. Tier 2: 糖脂捕获 (Glycolipid Catcher)
# =====================================================================

# 糖脂 SMARTS 骨架 (Glycolipid SMARTS patterns)
GLYCOLIPID_SMARTS = {
    # 神经酰胺 (Ceramide): 长链碱基 + 酰胺键 + 长链脂肪酸
    # Sphingoid base + amide bond + fatty acid chain
    "Sphingolipid": "CCCCCCCC(O)C(NC(=O)CCCC)CO",

    # 甘油酯 (Glycerolipid ester): 甘油骨架 + 至少两个酯键 + 长链
    # Glycerol backbone with at least two ester bonds
    "Glycerolipid_ester": "OCC(OC(=O)CCCCCCCC)COC(=O)CCCCCCCC",

    # 鞘氨醇 (Sphingosine backbone) — 含长链不饱和
    "Sphingoid_base": "CCCCCCCC/C=C/C(O)C(N)CO",
}

# 预编译
_COMPILED_GLYCOLIPID: Dict[str, Chem.Mol] = {}


def _getCompiledGlycolipid() -> Dict[str, Optional[Chem.Mol]]:
    """懒加载糖脂 SMARTS (Lazy-load glycolipid SMARTS)."""
    global _COMPILED_GLYCOLIPID
    if not _COMPILED_GLYCOLIPID:
        for name, smarts in GLYCOLIPID_SMARTS.items():
            _COMPILED_GLYCOLIPID[name] = Chem.MolFromSmarts(smarts)
    return _COMPILED_GLYCOLIPID


def classifyGlycolipid(mol: Chem.Mol) -> Optional[str]:
    """
    糖脂特异性捕获: 检测长链脂肪烃 (>C10) 或神经酰胺/甘油骨架。
    Glycolipid-specific catcher: detect long aliphatic chains (>C10) or
    ceramide/glycerol backbones.

    判定逻辑 (Decision Logic):
    1. 检查是否含有已知糖脂骨架 (Ceramide, Glycerol)
    2. 否则，检查最长连续脂肪碳链是否 > 10

    Args:
        mol: RDKit Mol 对象

    Returns:
        糖脂亚类名或 None
    """
    # 1. 检查已知骨架 (Check known glycolipid skeletons)
    glycolipidSmarts = _getCompiledGlycolipid()
    for subtype, pattern in glycolipidSmarts.items():
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return f"Glycolipid ({subtype})"

    # 2. 长链脂肪烃检测 (Long aliphatic chain detection)
    # 构建碳原子邻接图，找最长纯碳链 (Build C-atom graph, find longest pure-C chain)
    longestChain = _findLongestAliphaticChain(mol)
    if longestChain > 10:
        return f"Glycolipid (Aliphatic C{longestChain})"

    return None


def _findLongestAliphaticChain(mol: Chem.Mol) -> int:
    """
    在分子中查找最长的连续脂肪碳链 (非芳香、非环内)。
    Find the longest contiguous aliphatic (non-aromatic, non-ring) carbon chain.

    使用 DFS 遍历所有 sp3/sp2 碳原子之间的链式连接。
    Uses DFS to traverse chain connections between sp3/sp2 carbon atoms.

    Args:
        mol: RDKit Mol 对象

    Returns:
        最长脂肪碳链长度
    """
    ri = mol.GetRingInfo()
    ringAtoms = set()
    for ring in ri.AtomRings():
        ringAtoms.update(ring)

    # 收集非环、非芳香的碳原子 (Collect non-ring, non-aromatic carbon atoms)
    aliphaticCarbons = set()
    for atom in mol.GetAtoms():
        if (atom.GetAtomicNum() == 6
                and atom.GetIdx() not in ringAtoms
                and not atom.GetIsAromatic()):
            aliphaticCarbons.add(atom.GetIdx())

    if not aliphaticCarbons:
        return 0

    # 构建邻接表 (Build adjacency list for aliphatic carbons only)
    adjacency: Dict[int, List[int]] = {idx: [] for idx in aliphaticCarbons}
    for idx in aliphaticCarbons:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            nbrIdx = nbr.GetIdx()
            if nbrIdx in aliphaticCarbons:
                adjacency[idx].append(nbrIdx)

    # DFS 找最长路径 (DFS for longest path)
    maxLength = 0

    def dfs(node: int, visited: set) -> int:
        best = 0
        for nbrIdx in adjacency[node]:
            if nbrIdx not in visited:
                visited.add(nbrIdx)
                length = 1 + dfs(nbrIdx, visited)
                best = max(best, length)
                visited.discard(nbrIdx)
        return best

    for startNode in aliphaticCarbons:
        length = 1 + dfs(startNode, {startNode})
        maxLength = max(maxLength, length)

    return maxLength


# =====================================================================
# 4. Tier 3: Tanimoto 最近邻救援 (Tanimoto Nearest-Neighbor Rescue)
# =====================================================================

def buildReferenceFingerprints(
    df: pd.DataFrame,
    smilesCol: str = "Aglycan_SMILES",
    classCol: str = "np_classifier_superclass",
) -> Tuple[List[DataStructs.ExplicitBitVect], List[str]]:
    """
    从已有分类数据构建参考指纹库。
    Build reference fingerprint library from existing classified data.

    Args:
        df: 含有 SMILES 和分类列的 DataFrame
        smilesCol: SMILES 列名
        classCol: 分类列名

    Returns:
        (fingerprint_list, class_label_list)
    """
    refFps = []
    refLabels = []

    for _, row in df.iterrows():
        smiles = str(row.get(smilesCol, ""))
        label = str(row.get(classCol, ""))

        if not smiles or smiles in ("NULL", "nan", ""):
            continue
        if not label or label in ("NULL", "nan", ""):
            continue

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        try:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
            refFps.append(fp)
            refLabels.append(label)
        except Exception:
            continue

    return refFps, refLabels


def classifyByTanimoto(
    mol: Chem.Mol,
    referenceFps: List,
    referenceLabels: List[str],
    threshold: float = 0.85,
) -> Optional[str]:
    """
    通过 Tanimoto 最近邻搜索进行分类继承。
    Classify by Tanimoto nearest-neighbor search.

    Args:
        mol: 待分类分子
        referenceFps: 参考指纹列表
        referenceLabels: 参考分类标签列表
        threshold: 相似度阈值 (default=0.85)

    Returns:
        最近邻分类名或 None
    """
    if not referenceFps:
        return None

    try:
        queryFp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    except Exception:
        return None

    similarities = DataStructs.BulkTanimotoSimilarity(queryFp, referenceFps)

    maxSim = max(similarities)
    if maxSim >= threshold:
        bestIdx = similarities.index(maxSim)
        return f"{referenceLabels[bestIdx]} (Tanimoto={maxSim:.3f})"

    return None


# =====================================================================
# 5. 统一分类入口 (Unified Classification Entry Point)
# =====================================================================

def classifyAglycon(
    smiles: str,
    referenceFps: Optional[List] = None,
    referenceLabels: Optional[List[str]] = None,
    tanimotoThreshold: float = 0.85,
) -> Dict[str, str]:
    """
    三层分类: SMARTS → 糖脂捕获 → Tanimoto 救援。
    Three-tier classification: SMARTS → Glycolipid → Tanimoto rescue.

    Args:
        smiles: 苷元 SMILES
        referenceFps: Tanimoto 参考指纹库 (可选)
        referenceLabels: Tanimoto 参考标签库 (可选)
        tanimotoThreshold: Tanimoto 阈值

    Returns:
        字典包含: classification, classification_method, glycolipid_flag
    """
    result = {
        "classification": "Unclassified",
        "classification_method": "none",
        "glycolipid_flag": "",
    }

    if not smiles or smiles in ("NULL", "nan", "", "*"):
        return result

    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return result

    # Tier 2: 糖脂捕获 — 独立运行，不依赖 Tier 1 结果
    # Tier 2: Glycolipid catcher — runs independently
    glycolipidClass = classifyGlycolipid(mol)
    if glycolipidClass:
        result["glycolipid_flag"] = glycolipidClass

    # Tier 1: SMARTS 骨架匹配 (SMARTS skeleton matching)
    smartsClass = classifyBySMARTS(mol)
    if smartsClass:
        result["classification"] = smartsClass
        result["classification_method"] = "SMARTS"
        return result

    # 如果 Tier 1 未命中但 Tier 2 命中，直接用糖脂分类
    # If Tier 1 missed but Tier 2 hit, use glycolipid classification
    if glycolipidClass:
        result["classification"] = glycolipidClass
        result["classification_method"] = "Glycolipid_Rule"
        return result

    # Tier 3: Tanimoto 最近邻救援 (Tanimoto nearest-neighbor rescue)
    if referenceFps and referenceLabels:
        tanimotoClass = classifyByTanimoto(
            mol, referenceFps, referenceLabels, tanimotoThreshold
        )
        if tanimotoClass:
            result["classification"] = tanimotoClass
            result["classification_method"] = "Tanimoto"
            return result

    return result
