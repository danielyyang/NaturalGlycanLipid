"""
糖链修饰基团 SMARTS 扫描器
Glycan Modification SMARTS Scanner

扫描糖链片段上的常见化学修饰基团并生成格式化标注。
Scan glycan fragments for common chemical modifications and produce
formatted annotation strings.

已整合的修饰基团 (Integrated Modification Groups):
  - O-Acetyl (O-Ac):    酯键乙酰化
  - N-Acetyl (NAc):     酰胺键乙酰化 (GlcNAc 等)
  - O-Methyl (O-Me):    甲基化
  - Sulfate (S):        硫酸酯化
  - Phosphate (P):      磷酸酯化
  - Carboxyl (COOH):    羧基 (醛糖酸)
  - Deoxy:              脱氧 (鼠李糖等)
  - Amino (NH2):        游离氨基
"""
import os
import sys
from typing import Dict, List, Optional, Tuple
from collections import OrderedDict

from rdkit import Chem

# =====================================================================
# 1. SMARTS 修饰基团字典 (Modification SMARTS Dictionary)
# =====================================================================
# 修饰 SMARTS 字典现已统一定义在 glycan_reference_library.py 中。
# Modification SMARTS are now centralized in glycan_reference_library.py.
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__))))
try:
    from glycan_reference_library import MODIFICATION_SMARTS
except ImportError:
    from lib.glycan_reference_library import MODIFICATION_SMARTS

# 预编译 (Pre-compile)
_COMPILED_MOD_SMARTS: Dict[str, Chem.Mol] = {}


def _getCompiledModSmarts() -> Dict[str, Optional[Chem.Mol]]:
    """懒加载预编译修饰 SMARTS (Lazy-load compiled modification SMARTS)."""
    global _COMPILED_MOD_SMARTS
    if not _COMPILED_MOD_SMARTS:
        for name, smarts in MODIFICATION_SMARTS.items():
            mol = Chem.MolFromSmarts(smarts)
            if mol is None:
                print(f"[WARNING] Failed to compile SMARTS for '{name}': {smarts}")
            _COMPILED_MOD_SMARTS[name] = mol
    return _COMPILED_MOD_SMARTS


# =====================================================================
# 2. 单分子扫描 (Single-Molecule Scan)
# =====================================================================

def scanGlycanModifications(smiles: str) -> Dict[str, int]:
    """
    扫描糖链 SMILES 中的所有修饰基团并计数。
    Scan glycan SMILES for all modification groups and count matches.

    Args:
        smiles: 糖链 SMILES (Phase 2 切分产物)

    Returns:
        OrderedDict {修饰名: 匹配次数} (只含 count>0 的)
    """
    result: Dict[str, int] = OrderedDict()

    if not smiles or smiles in ("NULL", "nan", "", "*"):
        return result

    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return result

    smartsDict = _getCompiledModSmarts()

    for modName, pattern in smartsDict.items():
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            result[modName] = len(matches)

    # ─── 芳香酯去重 (Aromatic Ester Deduplication) ───────────────────
    # 设计意图: O-Bz 的 SMARTS 是 O-Gall/O-pCou/O-Caf/O-Fer 的子结构,
    # 同一个酰基氧原子会同时命中 O-Bz 和更特异的模式。需减去重叠计数。
    # Design: O-Bz SMARTS is a subgraph of O-Gall/O-pCou/O-Caf/O-Fer.
    # Same acyl oxygen hits both. Subtract overlap from O-Bz count.
    SPECIFIC_AROMATIC_ESTERS = ("O-Gall", "O-pCou", "O-Caf", "O-Fer")
    if "O-Bz" in result:
        overlap = sum(result.get(k, 0) for k in SPECIFIC_AROMATIC_ESTERS)
        result["O-Bz"] = max(0, result["O-Bz"] - overlap)
        if result["O-Bz"] == 0:
            del result["O-Bz"]

    return result


def scanWithAtomIndices(smiles: str) -> Tuple[Dict[str, int], set]:
    """
    扫描修饰基团并返回修饰原子的索引集合 (用于可视化高亮)。
    Scan modifications and return atom indices of modification groups.

    Args:
        smiles: 糖链 SMILES

    Returns:
        (modCounts, modificationAtomIndices)
        modCounts: {修饰名: 匹配数}
        modificationAtomIndices: 所有修饰基团涉及的原子索引集合
    """
    modCounts: Dict[str, int] = OrderedDict()
    modAtomIndices: set = set()

    if not smiles or smiles in ("NULL", "nan", "", "*"):
        return modCounts, modAtomIndices

    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return modCounts, modAtomIndices

    smartsDict = _getCompiledModSmarts()

    for modName, pattern in smartsDict.items():
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            modCounts[modName] = len(matches)
            for match in matches:
                modAtomIndices.update(match)

    return modCounts, modAtomIndices


def formatModifications(modCounts: Dict[str, int]) -> str:
    """
    将修饰计数格式化为可读字符串。
    Format modification counts into readable string.

    格式 (Format): "O-Ac: 2, NAc: 1"
    如果无修饰则返回空字符串。

    Args:
        modCounts: {修饰名: 计数}

    Returns:
        格式化字符串, 或 ""
    """
    if not modCounts:
        return ""
    parts = [f"{name}: {count}" for name, count in modCounts.items()]
    return ", ".join(parts)


def scanAndFormat(smiles: str) -> str:
    """
    一步完成: 扫描 + 格式化。
    One-step: scan + format.

    Args:
        smiles: 糖链 SMILES

    Returns:
        "O-Ac: 2, NAc: 1" 或 ""
    """
    return formatModifications(scanGlycanModifications(smiles))


# =====================================================================
# 3. 批量处理 (Batch Processing)
# =====================================================================

def batchScanModifications(
    df: "pd.DataFrame",
    glycanCol: str = "Glycan_SMILES",
    outputCol: str = "Glycan_Modifications",
) -> "pd.DataFrame":
    """
    对 DataFrame 批量扫描糖链修饰，追加新列。
    Batch scan glycan modifications for a DataFrame, appending a new column.

    Args:
        df: 输入 DataFrame
        glycanCol: 糖链 SMILES 列名
        outputCol: 输出修饰列名

    Returns:
        添加了修饰列的 DataFrame
    """
    from tqdm import tqdm

    results = []
    for smiles in tqdm(df[glycanCol].fillna(""), desc="Scanning Glycan Modifications"):
        results.append(scanAndFormat(str(smiles)))

    df[outputCol] = results
    return df
