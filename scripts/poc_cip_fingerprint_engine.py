"""
动态 R/S 指纹单糖鉴定引擎 PoC (Dynamic R/S Fingerprint Identification Engine PoC)
==================================================================================

核心洞察 (Core Insight):
  SubstructMatch(useChirality=True) 在虚拟脱修饰后**根本不可行**, 因为 CIP 优先级
  会随着取代基变化而翻转 (例如 -N-SO3H → -NH2 导致 C2 的 CIP 从 S 变成 R)。

  SubstructMatch fails after demodification because CIP priorities at chiral
  centers change when substituents are removed (e.g., -N-SO3H → -NH2 flips C2).

解决方案 (Solution):
  拓扑无关的环行走算法 (Topology-Independent Ring-Walking):
  1. 找到糖环 O
  2. 沿环行走: O → C1 → C2 → C3 → C4 → C5
  3. 在**脱修饰前的原始分子**上提取 C2-C5 的 CIP
     (C1 = 异头碳, CIP 依赖异头构型, 不用于身份匹配)
     (C2 对于氨基糖CIP会变, 但 ChiralTag 不变)
  4. 用 ChiralTag(CW/CCW) + 取代基拓扑对比参考库

  关键区别: 不是比 CIP R/S label, 而是比 ChiralTag (原子本身的四面体方向)
  + 环上邻居的原子类型。ChiralTag 不受远端取代基变化影响。

Usage:
  python scripts/poc_cip_fingerprint_engine.py
"""
import sys
import os
from typing import Dict, List, Optional, Tuple, Set

from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from lib.glycan_reference_library import RAW_MONOSACCHARIDE_SMILES


# =====================================================================
# Ring Walk: 糖环标准路径提取 (Sugar Ring Canonical Path Extraction)
# =====================================================================

def findSugarRingPath(mol: Chem.Mol, ringAtoms: Optional[List[int]] = None) -> Optional[List[int]]:
    """沿糖环行走: O → C1 → C2 → C3 → C4 → C5, 返回有序索引列表。
    Walk sugar ring: O → C1 → C2 → C3 → C4 → C5, return ordered index list.

    C1 定义: 与环 O 相邻的碳, 且有环外氧 (异头碳 / anomeric carbon)。
    如果两个碳都有环外 O, 选择总连接度更高的。
    """
    if ringAtoms is None:
        # 自动查找第一个 6-元环 (Auto-find first 6-membered ring)
        rings = mol.GetRingInfo().AtomRings()
        sixRings = [r for r in rings if len(r) == 6]
        if not sixRings:
            return None
        ringAtoms = list(sixRings[0])

    ringSet = set(ringAtoms)

    # 找环 O (Find ring oxygen)
    ringOIdx = None
    for idx in ringAtoms:
        if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8:
            ringOIdx = idx
            break
    if ringOIdx is None:
        return None

    # 找 O 的两个环碳邻居 (Find O's two ring-carbon neighbors)
    oAtom = mol.GetAtomWithIdx(ringOIdx)
    ringCarbonNbrs = [n.GetIdx() for n in oAtom.GetNeighbors()
                      if n.GetIdx() in ringSet and n.GetAtomicNum() == 6]
    if len(ringCarbonNbrs) != 2:
        return None

    # C1 = 有环外 O 的碳 (anomeric carbon has exocyclic O)
    # C1 heuristic: the carbon with exocyclic O (or more exocyclic O's)
    def countExoO(cIdx: int) -> int:
        atom = mol.GetAtomWithIdx(cIdx)
        return sum(1 for n in atom.GetNeighbors()
                   if n.GetAtomicNum() == 8 and n.GetIdx() not in ringSet)

    exoO0 = countExoO(ringCarbonNbrs[0])
    exoO1 = countExoO(ringCarbonNbrs[1])

    if exoO0 > exoO1:
        c1Idx = ringCarbonNbrs[0]
    elif exoO1 > exoO0:
        c1Idx = ringCarbonNbrs[1]
    else:
        # 平局: 选总连接度更高的 (Tie-break: higher total degree)
        c1Idx = ringCarbonNbrs[0]

    # 行走: C1 → C2 → C3 → C4 → C5 (Walk ring)
    path = [c1Idx]
    prev = ringOIdx
    curr = c1Idx
    for _ in range(4):
        atom = mol.GetAtomWithIdx(curr)
        nextC = None
        for n in atom.GetNeighbors():
            nIdx = n.GetIdx()
            if nIdx in ringSet and nIdx != prev and n.GetAtomicNum() == 6:
                nextC = nIdx
                break
        if nextC is None:
            # 遇到第二个 O (走到了 C5-O) — 停止
            break
        path.append(nextC)
        prev = curr
        curr = nextC

    return path


# =====================================================================
# Chirality Fingerprint: 基于 ChiralTag 而非 CIP
# =====================================================================

def extractChiralFingerprint(
    mol: Chem.Mol,
    ringPath: List[int],
) -> Dict:
    """提取环碳 C1-C5 的立体化学指纹。
    Extract stereochemistry fingerprint at ring carbons C1-C5.

    返回 (Returns):
      Dict with:
        'cip': [C1_CIP, C2_CIP, C3_CIP, C4_CIP, C5_CIP]
        'tag': [C1_tag, C2_tag, C3_tag, C4_tag, C5_tag]
        'exo': [(element, ...), ...] — 环外取代基元素列表
    """
    ringSet = set(ringPath)
    # 加入环 O (infer from ring)
    for idx in ringPath:
        atom = mol.GetAtomWithIdx(idx)
        for n in atom.GetNeighbors():
            if n.GetAtomicNum() == 8 and n.GetIdx() not in ringSet:
                continue
            if n.GetAtomicNum() == 8:
                ringSet.add(n.GetIdx())

    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    cipList = []
    tagList = []
    exoList = []

    for idx in ringPath:
        atom = mol.GetAtomWithIdx(idx)

        # CIP code
        cip = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else "?"
        cipList.append(cip)

        # ChiralTag (CW/CCW) — 不受远端取代基影响
        # ChiralTag (CW/CCW) — NOT affected by distal substituent changes
        tag = str(atom.GetChiralTag())
        tagList.append(tag)

        # 环外取代基 (Exocyclic substituent element types)
        exo = []
        for n in atom.GetNeighbors():
            if n.GetIdx() not in ringSet:
                exo.append(n.GetAtomicNum())
        exoList.append(tuple(sorted(exo)))

    return {
        "cip": tuple(cipList),
        "tag": tuple(tagList),
        "exo": tuple(exoList),
    }


# =====================================================================
# Building Reference Fingerprint Database
# =====================================================================

def buildReferenceFingerprints() -> Dict[Tuple[str, str], Dict]:
    """从 RAW_MONOSACCHARIDE_SMILES 构建参考指纹库。
    Build reference fingerprint database from the sugar SMILES library.
    """
    refDb = {}
    for (name, anomer), smi in RAW_MONOSACCHARIDE_SMILES.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        path = findSugarRingPath(mol)
        if path is None or len(path) < 4:
            continue
        fp = extractChiralFingerprint(mol, path)
        refDb[(name, anomer)] = fp
    return refDb


def matchByCipFingerprint(
    targetFp: Dict,
    refDb: Dict[Tuple[str, str], Dict],
    skipC1: bool = True,
    skipC2Amino: bool = True,
) -> List[Tuple[str, str, int, int]]:
    """用 C2-C5 的 CIP 指纹匹配参考库。
    Match target sugar's CIP fingerprint against reference database.

    Design Intent:
      C1 (anomeric) is skipped by default — it changes with alpha/beta.
      C2 for amino sugars: CIP flips with substituent changes,
        but we can compare if both are '?' or if the ChiralTag matches.

    Returns list of (name, anomer, matchScore, totalCheckable) sorted by score desc.
    """
    targetCip = targetFp["cip"]
    targetTag = targetFp["tag"]
    targetExo = targetFp["exo"]

    results = []
    for (refName, refAnomer), refFp in refDb.items():
        refCip = refFp["cip"]
        refTag = refFp["tag"]

        if len(targetCip) != len(refCip):
            continue

        matchCount = 0
        checkable = 0

        for i in range(len(targetCip)):
            if skipC1 and i == 0:
                continue  # C1 异头碳 — 跳过

            tCip = targetCip[i]
            rCip = refCip[i]

            # 如果任一方 CIP='?', 不可比较 — 跳过
            if tCip == "?" or rCip == "?":
                continue

            checkable += 1
            if tCip == rCip:
                matchCount += 1

        if checkable == 0:
            continue

        results.append((refName, refAnomer, matchCount, checkable))

    # 排序: 先按匹配数降序, 再按 checkable 降序
    results.sort(key=lambda x: (x[2], x[3]), reverse=True)
    return results


# =====================================================================
# Test Suite
# =====================================================================

def main():
    print("=" * 70)
    print("  Dynamic R/S Fingerprint Engine PoC")
    print("  动态 R/S 指纹鉴定引擎 PoC")
    print("=" * 70)

    # Build reference database
    refDb = buildReferenceFingerprints()
    print(f"\n  Reference fingerprints built: {len(refDb)} entries")

    # Display a few key references
    for name in ["D-Glc", "D-GlcN", "D-GlcNAc", "D-AllN", "D-Gal", "D-GalN",
                  "L-Rha", "L-Fuc", "D-Man", "L-IdoA", "D-GlcA"]:
        for anom in ["a", "b"]:
            key = (name, anom)
            if key in refDb:
                fp = refDb[key]
                print(f"  {name:10s}({anom}): CIP={fp['cip']}  Exo={fp['exo']}")

    # === Test cases ===
    testCases = [
        (
            "D-GlcNAc 6-Sulfate → demod → bare D-GlcN",
            # 先脱修饰再提指纹
            "CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)O)[C@@H](O)[C@@H]1O",
            ["O-Sulfate", "N-Acetyl"],  # 全部剥离
            "D-GlcN",
        ),
        (
            "Fondaparinux片段 (N-Sulfated GlcN)",
            "O[C@H]1[C@@H](O)[C@@H](NS(=O)(=O)O)[C@H](O)O[C@@H]1COS(=O)(=O)O",
            ["O-Sulfate", "N-Sulfate"],
            "D-GlcN",
        ),
        (
            "Bare D-Glc (no mods)",
            "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
            [],
            "D-Glc",
        ),
        (
            "Bare L-Rha (no mods)",
            "C[C@@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O",
            [],
            "L-Rha",
        ),
    ]

    # Import demod from PoC
    from poc_virtual_demodify import virtualDemodify

    print("\n" + "-" * 70)
    for testName, smi, modsToDo, expectedName in testCases:
        print(f"\n  === {testName} ===")
        print(f"  Input: {smi}")

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print("  PARSE FAIL")
            continue

        # === STEP 1: 在原始分子上提取 CIP 指纹 ===
        # 关键: CIP 在原始分子上计算, 不在脱修饰后计算!
        path1 = findSugarRingPath(mol)
        if path1:
            fpOriginal = extractChiralFingerprint(mol, path1)
            print(f"  Original CIP: {fpOriginal['cip']}")

        # === STEP 2: 脱修饰 ===
        if modsToDo:
            cleanMol, detectedMods = virtualDemodify(mol, modNames=modsToDo)
            print(f"  Demod mods: {detectedMods}")
            print(f"  Demod SMILES: {Chem.MolToSmiles(cleanMol)}")

            # 在脱修饰后的分子上也提取 CIP (用于对比)
            path2 = findSugarRingPath(cleanMol)
            if path2:
                fpDemod = extractChiralFingerprint(cleanMol, path2)
                print(f"  Demod CIP:    {fpDemod['cip']}")
        else:
            fpDemod = fpOriginal if path1 else None

        # === STEP 3: 指纹匹配 ===
        # 用原始 CIP 匹配 (避免脱修饰后的 CIP 翻转)
        if path1 and fpOriginal:
            matches = matchByCipFingerprint(fpOriginal, refDb)
            topMatches = matches[:5]
            print(f"  Top matches (original CIP, C2-C5 only):")
            for mName, mAnomer, mScore, mCheck in topMatches:
                isPerfect = "***" if mScore == mCheck else "   "
                isExpected = " <<<" if mName == expectedName else ""
                print(f"    {isPerfect} {mName:10s}({mAnomer}): {mScore}/{mCheck}{isExpected}")

        # 用脱修饰 CIP 也试一次
        if modsToDo and path2 and fpDemod:
            matches2 = matchByCipFingerprint(fpDemod, refDb)
            topMatches2 = matches2[:5]
            print(f"  Top matches (demod CIP, C2-C5 only):")
            for mName, mAnomer, mScore, mCheck in topMatches2:
                isPerfect = "***" if mScore == mCheck else "   "
                isExpected = " <<<" if mName == expectedName else ""
                print(f"    {isPerfect} {mName:10s}({mAnomer}): {mScore}/{mCheck}{isExpected}")

    print(f"\n{'='*70}")
    print("  PoC Complete")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
