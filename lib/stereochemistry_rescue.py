"""
立体化学确证流水线 (Stereo-Rescue Pipeline)
============================================

多维立体化学确证模块，用于解决 2D SMILES 丢失手性导致的糖链识别模糊问题。
Multi-dimensional stereo-rescue pipeline for resolving ambiguous sugar
identifications (Hex, Pen, Non) caused by lost chirality in 2D SMILES.

模块 1: ChEMBL InChIKey 结构换头术 — 用高质量 SMILES 替换丢失手性的结构
Module 1: ChEMBL InChIKey Structure Replacement — replace with chirality-rich SMILES

模块 2: LOTUS InChIKey 交叉验证 — 利用 LOTUS (Wiki) 数据库进行交叉匹配
Module 2: LOTUS InChIKey Cross-validation — cross-match via LOTUS (Wiki) database

模块 3: IUPAC/通用名糖词提取 — 从名称字段正则提取糖的确切类型
Module 3: IUPAC/Name Sugar Word Extraction — regex extract sugar types from name fields
"""
import os
import re
import sqlite3
import gzip
from typing import Optional, Dict, List, Tuple
from functools import lru_cache

import pandas as pd
from rdkit import Chem
from rdkit.Chem.inchi import MolFromInchi, MolToInchi, InchiToInchiKey


# =====================================================================
# 全局配置 (Global Configuration)
# =====================================================================

# ChEMBL SQLite 路径 (ChEMBL SQLite database path)
CHEMBL_DB_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data", "chembl_36", "chembl_36_sqlite", "chembl_36.db"
)

# LOTUS 压缩 CSV 路径 (LOTUS compressed CSV path)
LOTUS_CSV_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "data", "230106_frozen.csv.gz"
)

# 泛指糖标签 — 需要被 rescue 的标签 (Generic sugar labels needing rescue)
GENERIC_SUGAR_LABELS = frozenset({"Hex", "Pen", "Non", "Oct", "Hept", "dHex", "HexN"})

# IUPAC 糖根词映射白名单 (Sugar root-word → standard name mapping)
# 设计原则: 只映射确定性最高的根词, 避免假阳性
# Design: only map the most unambiguous root words
SUGAR_ROOT_WORDS: Dict[str, str] = {
    "glucopyran":    "D-Glc",
    "glucofuran":    "D-Glc",
    "gluco":         "D-Glc",
    "galactopyran":  "D-Gal",
    "galactofuran":  "D-Gal",
    "galacto":       "D-Gal",
    "mannopyran":    "D-Man",
    "mannofuran":    "D-Man",
    "manno":         "D-Man",
    "rhamnopyran":   "L-Rha",
    "rhamno":        "L-Rha",
    "arabinopyran":  "L-Ara",
    "arabinofuran":  "L-Ara",
    "arabino":       "L-Ara",
    "xylopyran":     "D-Xyl",
    "xylofuran":     "D-Xyl",
    "xylo":          "D-Xyl",
    "fructopyran":   "D-Fru",
    "fructofuran":   "D-Fru",
    "fructo":        "D-Fru",
    "ribopyran":     "D-Rib",
    "ribofuran":     "D-Rib",
    "ribo":          "D-Rib",
    "fucopyran":     "L-Fuc",
    "fucofuran":     "L-Fuc",
    "fuco":          "L-Fuc",
    "glucurono":     "D-GlcA",
    "glucuronic":    "D-GlcA",
    "galacturono":   "D-GalA",
    "galacturonic":  "D-GalA",
    "quinovopyran":  "D-Qui",
    "quinovo":       "D-Qui",
}

# 预编译正则: 按长度降序排列，确保最长匹配优先
# Pre-compiled regex: sorted by length descending for longest-match-first
_SUGAR_ROOT_PATTERN = re.compile(
    "|".join(sorted(SUGAR_ROOT_WORDS.keys(), key=len, reverse=True)),
    re.IGNORECASE
)


# =====================================================================
# 模块 1: ChEMBL InChIKey 结构换头术 (ChEMBL InChIKey Structure Replacement)
# =====================================================================

def _getChemblConnection() -> Optional[sqlite3.Connection]:
    """
    获取 ChEMBL SQLite 连接 (Get ChEMBL SQLite connection).
    如果数据库不存在或打开失败, 返回 None。
    Returns None if database doesn't exist or can't be opened.
    """
    if not os.path.exists(CHEMBL_DB_PATH):
        return None
    try:
        conn = sqlite3.connect(CHEMBL_DB_PATH, timeout=10)
        conn.execute("PRAGMA query_only = ON")   # 只读模式 (read-only mode)
        conn.execute("PRAGMA cache_size = -50000")  # 50MB 缓存 (cache)
        return conn
    except sqlite3.Error:
        return None


def rescueViaChembl(smiles: str) -> Optional[str]:
    """
    通过 InChIKey 第一区块在 ChEMBL 中查找带完整立体化学的 SMILES。
    Look up a chirality-rich SMILES in ChEMBL using the first InChIKey block.

    Args:
        smiles: 当前数据库中的 canonical SMILES (可能缺失手性)

    Returns:
        带完整立体化学的 SMILES 字符串, 或 None (如果未找到或无改善)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    try:
        inchi = MolToInchi(mol)
        if inchi is None:
            return None
        inchiKey = InchiToInchiKey(inchi)
        if inchiKey is None:
            return None
    except Exception:
        return None

    # 提取第一区块 (前14字符, 代表 2D 骨架, 不含立体化学)
    # Extract first block (14 chars, represents 2D skeleton without stereo)
    firstBlock = inchiKey.split("-")[0]

    conn = _getChemblConnection()
    if conn is None:
        return None

    try:
        cursor = conn.execute(
            "SELECT canonical_smiles FROM compound_structures "
            "WHERE standard_inchi_key LIKE ? || '%'",
            (firstBlock,)
        )
        candidates = [row[0] for row in cursor.fetchall() if row[0]]
    except sqlite3.Error:
        return None
    finally:
        conn.close()

    if not candidates:
        return None

    # 选择立体中心数量最多的候选分子
    # Select the candidate with the most defined stereocenters
    bestSmiles = None
    bestStereoCount = -1
    for candSmi in candidates:
        candMol = Chem.MolFromSmiles(candSmi)
        if candMol is None:
            continue
        # 计算已定义的立体中心数量
        # Count defined stereocenters
        Chem.AssignStereochemistry(candMol, cleanIt=True, force=True)
        stereoCount = sum(
            1 for atom in candMol.GetAtoms()
            if atom.HasProp('_CIPCode')
        )
        if stereoCount > bestStereoCount:
            bestStereoCount = stereoCount
            bestSmiles = candSmi

    # 只有当 ChEMBL 的版本确实有更多立体中心时才替换
    # Only replace if ChEMBL version genuinely has more stereocenters
    if bestSmiles and bestStereoCount > 0:
        origMol = Chem.MolFromSmiles(smiles)
        if origMol:
            Chem.AssignStereochemistry(origMol, cleanIt=True, force=True)
            origStereo = sum(
                1 for atom in origMol.GetAtoms()
                if atom.HasProp('_CIPCode')
            )
            if bestStereoCount > origStereo:
                return bestSmiles

    return None


# =====================================================================
# 模块 2: LOTUS InChIKey 交叉验证 (LOTUS InChIKey Cross-validation)
# =====================================================================

# 全局缓存: LOTUS InChIKey → organism 映射
# Global cache: LOTUS InChIKey → organism mapping
_LOTUS_INCHIKEY_SET: Optional[set] = None


def _loadLotusInchiKeys() -> set:
    """
    加载 LOTUS 数据库的所有 InChIKey (Load all InChIKeys from LOTUS database).
    LOTUS 包含经过人工验证的天然产物, 其 InChIKey 是可信的。
    LOTUS contains manually validated natural products; its InChIKeys are reliable.
    """
    global _LOTUS_INCHIKEY_SET
    if _LOTUS_INCHIKEY_SET is not None:
        return _LOTUS_INCHIKEY_SET

    if not os.path.exists(LOTUS_CSV_PATH):
        _LOTUS_INCHIKEY_SET = set()
        return _LOTUS_INCHIKEY_SET

    try:
        lotusDF = pd.read_csv(
            LOTUS_CSV_PATH,
            compression='gzip',
            usecols=['structure_inchikey'],
            dtype={'structure_inchikey': str}
        )
        # 提取第一区块 (前14字符) 用于模糊匹配
        # Extract first block (14 chars) for fuzzy matching
        _LOTUS_INCHIKEY_SET = set(
            key.split("-")[0]
            for key in lotusDF['structure_inchikey'].dropna().unique()
            if isinstance(key, str) and "-" in key
        )
        print(f"  [Stereo-Rescue] Loaded {len(_LOTUS_INCHIKEY_SET):,} LOTUS InChIKey first-blocks")
    except Exception as e:
        print(f"  [Stereo-Rescue] Failed to load LOTUS: {e}")
        _LOTUS_INCHIKEY_SET = set()

    return _LOTUS_INCHIKEY_SET


def isInLotus(smiles: str) -> bool:
    """
    检查该分子是否存在于 LOTUS 数据库中 (基于 InChIKey 第一区块匹配)。
    Check if this molecule exists in the LOTUS database (via InChIKey first-block match).

    用于筛选: 如果在 LOTUS 中找到, 则该分子是已知天然产物, 值得进一步 rescue。
    Filtering use: if found in LOTUS, the molecule is a known NP, worth further rescue.
    """
    lotusKeys = _loadLotusInchiKeys()
    if not lotusKeys:
        return False

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False

    try:
        inchi = MolToInchi(mol)
        if inchi is None:
            return False
        inchiKey = InchiToInchiKey(inchi)
        if inchiKey is None:
            return False
        firstBlock = inchiKey.split("-")[0]
        return firstBlock in lotusKeys
    except Exception:
        return False


# =====================================================================
# 模块 3: IUPAC/通用名糖词提取 (IUPAC/Name Sugar Word Extraction)
# =====================================================================

def rescueViaName(
    name: Optional[str],
    iupacName: Optional[str]
) -> Dict[str, str]:
    """
    从化合物名称和 IUPAC 名中提取糖类根词, 推断泛指标签的确切身份。
    Extract sugar root words from compound name and IUPAC name to infer
    the exact identity of generic sugar labels.

    Args:
        name: 化合物通用名 (compound common name)
        iupacName: IUPAC 系统名 (IUPAC systematic name)

    Returns:
        字典: {糖根词 → 标准糖名}, e.g. {"gluco": "D-Glc", "rhamno": "L-Rha"}
        Dict: {sugar_root_word → standard_sugar_name}
    """
    combinedText = " ".join(filter(None, [
        str(name) if name and str(name) != "nan" else None,
        str(iupacName) if iupacName and str(iupacName) != "nan" else None,
    ]))

    if not combinedText.strip():
        return {}

    matches = _SUGAR_ROOT_PATTERN.findall(combinedText)
    if not matches:
        return {}

    # 去重并保留映射关系 (Deduplicate and keep mapping)
    result: Dict[str, str] = {}
    for match in matches:
        matchLower = match.lower()
        for root, sugarName in SUGAR_ROOT_WORDS.items():
            if matchLower.startswith(root.lower()):
                result[root] = sugarName
                break

    return result


def replaceGenericLabelsInSequence(
    sequence: str,
    inferredSugars: Dict[str, str]
) -> str:
    """
    在糖序列字符串中替换泛指标签为推断出的具体糖名。
    Replace generic labels in sugar sequence with inferred specific sugar names.

    策略: 按名称在文本中出现的顺序, 与序列中泛指标签的位置进行对齐。
    Strategy: align inferred sugars by positional order in name text with
    generic labels in the sequence.

    修复: 先按 " ; " 分割多糖链, 再按糖苷键分割; 使用位置映射而非频率去重
    Fix: split by " ; " first, then by linkage; use positional mapping, not
    frequency-based dedup.

    Args:
        sequence: 原始糖序列, e.g. "Non ; Non" 或 "Hex-(a1-4)-Hex-(b1-3)-D-Glc"
        inferredSugars: 推断结果, e.g. {"arabinopyran": "L-Ara", "glucopyran": "D-Glc"}

    Returns:
        替换后的序列字符串
    """
    if not inferredSugars or not sequence:
        return sequence

    genericLabels = {"Hex", "Pen", "Non", "Oct", "Hept", "dHex", "HexN",
                     "HexA", "HexN", "HexNAc"}

    # 按名称出现顺序排列的推断列表 (保留位置对应)
    # Positional order: dict insertion order == order found in name text
    inferredList = list(inferredSugars.values())
    if not inferredList:
        return sequence

    # 去掉完全重复但保留顺序 (deduplicate preserving order)
    seen = set()
    uniqueInferred = []
    for s in inferredList:
        if s not in seen:
            seen.add(s)
            uniqueInferred.append(s)

    # 先按 " ; " 分割多糖链 (Split by " ; " first for multi-chain sequences)
    chains = sequence.split(" ; ")

    replaced = False
    sugarIdx = 0
    newChains = []

    for chain in chains:
        # 在每条链内, 按糖苷键分割 (Split within each chain by linkage pattern)
        tokens = re.split(r'(-\([^)]+\)-)', chain)
        newTokens = []
        for token in tokens:
            stripped = token.strip()
            if stripped in genericLabels and sugarIdx < len(uniqueInferred):
                newTokens.append(uniqueInferred[sugarIdx])
                sugarIdx += 1
                replaced = True
            elif stripped in genericLabels and uniqueInferred:
                # 超出推断数量时, 循环使用 (cycle if more generics than inferred)
                newTokens.append(uniqueInferred[sugarIdx % len(uniqueInferred)])
                sugarIdx += 1
                replaced = True
            else:
                newTokens.append(token)
        newChains.append("".join(newTokens))

    if replaced:
        return " ; ".join(newChains)
    return sequence


# =====================================================================
# 主入口: 综合 Rescue 流水线 (Main Entry: Combined Rescue Pipeline)
# =====================================================================

def rescueSugarSequence(
    smiles: str,
    currentSequence: str,
    name: Optional[str] = None,
    iupacName: Optional[str] = None
) -> Tuple[str, str]:
    """
    对含泛指标签的糖序列执行多维 rescue。
    Execute multi-dimensional rescue on sugar sequences containing generic labels.

    执行顺序 (Execution Order):
    1. ChEMBL InChIKey 换头 → 完全替换 SMILES 并重新分析
    2. IUPAC/Name 糖词提取 → 在序列层面替换泛指标签

    Args:
        smiles: 当前分子的 canonical SMILES
        currentSequence: 当前的 Sugar_Sequence 字符串
        name: 化合物通用名
        iupacName: IUPAC 系统名

    Returns:
        (rescuedSequence, rescueMethod) 元组
        rescuedSequence: rescue 后的序列 (如果无变化则与 currentSequence 相同)
        rescueMethod: 使用的 rescue 方法描述 ("ChEMBL", "Name-Inferred", "")
    """
    # 快速检查: 是否包含需要 rescue 的泛指标签
    # Quick check: does the sequence contain any generic labels?
    if not currentSequence or not any(
        label in currentSequence for label in GENERIC_SUGAR_LABELS
    ):
        return currentSequence, ""

    # --- 策略 1: IUPAC/Name 糖词提取 (快速, 无 IO) ---
    # --- Strategy 1: IUPAC/Name sugar word extraction (fast, no IO) ---
    inferredSugars = rescueViaName(name, iupacName)
    if inferredSugars:
        rescuedSeq = replaceGenericLabelsInSequence(currentSequence, inferredSugars)
        if rescuedSeq != currentSequence:
            return rescuedSeq, "Name-Inferred"

    # --- 策略 2: ChEMBL InChIKey 换头术 (需要 SQLite IO) ---
    # --- Strategy 2: ChEMBL InChIKey replacement (requires SQLite IO) ---
    # 注意: ChEMBL SQLite 30GB 无索引, LIKE 查询 22s/行, 暂时禁用
    # NOTE: ChEMBL SQLite 30GB has no index on standard_inchi_key.
    #       LIKE query takes 22s/row → disabled until index is built.
    #       To enable: CREATE INDEX idx_inchikey_prefix ON compound_structures
    #                  (substr(standard_inchi_key, 1, 14));
    # betterSmiles = rescueViaChembl(smiles)
    # if betterSmiles:
    #     return currentSequence, f"ChEMBL-Candidate:{betterSmiles}"

    return currentSequence, ""
