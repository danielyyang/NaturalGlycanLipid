"""
精确糖苷键切分引擎 — 基于异头碳检测
Precision Glycosidic Bond Cleavage Engine — Based on Anomeric Carbon Detection

核心设计原则 (Core Design Principles):
1. 仅切断糖苷键 (Glycosidic Bonds): C1(anomeric)-O/N/S-X
2. 绝不切断糖环内键 (C-C/C-O within ring) 或 C5-C6 exocyclic bond
3. 碳原子守恒断言: C_glycan + C_aglycan == C_total
4. 使用 Chem.FragmentOnBonds + Isotope-labeled dummy atoms 追踪连接关系

术语 (Terminology):
  - 异头碳 (Anomeric Carbon, C1): 糖环中连接环氧和另一个杂原子的碳
  - 糖苷键 (Glycosidic Bond): C1-O-R 中的 C1-O 键或 O-R 键
  - 桥原子 (Bridge Atom): 连接两个单元的 O/N/S 原子
"""
import os
import sys
from typing import List, Dict, Tuple, Optional, Set

from rdkit import Chem

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__))))


# =====================================================================
# 1. 异头碳检测 (Anomeric Carbon Detection)
# =====================================================================

def findAnomericCarbons(mol: Chem.Mol, ringAtoms: List[int]) -> List[int]:
    """
    在糖环中定位异头碳 (C1)。
    Locate anomeric carbon (C1) in a sugar ring.

    异头碳的判定标准 (Anomeric carbon criteria):
    - 属于糖环中的碳原子
    - 同时连接环内氧和至少一个环外杂原子 (O/N/S)
    - 因此有 ≥2 个杂原子邻居

    Args:
        mol: RDKit Mol 对象
        ringAtoms: 糖环原子索引列表

    Returns:
        异头碳索引列表 (通常只有 1 个)
    """
    ringSet = set(ringAtoms)
    anomericCarbons = []

    for idx in ringAtoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:  # 只看碳原子
            continue

        # 统计杂原子邻居数 (Count heteroatom neighbors)
        heteroNeighborCount = 0
        hasRingOxygen = False

        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() in (8, 7, 16):  # O, N, S
                heteroNeighborCount += 1
                if nbr.GetIdx() in ringSet:
                    hasRingOxygen = True

        # 异头碳: 在环中，连接环氧，且有 ≥2 个杂原子邻居
        # Anomeric carbon: in ring, connected to ring oxygen, ≥2 heteroatom neighbors
        if hasRingOxygen and heteroNeighborCount >= 2:
            anomericCarbons.append(idx)

    return anomericCarbons


# =====================================================================
# 2. 糖苷键定位 (Glycosidic Bond Localization)
# =====================================================================

def _canReachSugarRing(
    mol: Chem.Mol,
    startIdx: int,
    allSugarRingAtoms: Set[int],
    ringToUnit: Dict[int, int],
    sourceUnitIdx: int,
    maxSteps: int = 2,
) -> Tuple[bool, int]:
    """
    BFS 检查: 从 startIdx 出发，沿碳链 (C-C 键) 最多走 maxSteps 步，
    能否到达属于另一个糖环的原子。

    BFS check: from startIdx, follow C-C bonds up to maxSteps,
    can we reach an atom belonging to another sugar ring?

    这解决了 Rha-C1-O-C6(Glc) 的问题: C6 → C5(ring) 只需 1 步。
    This solves the Rha-C1-O-C6(Glc) problem: C6 → C5(ring) = 1 step.

    Args:
        mol: RDKit Mol
        startIdx: 起始碳原子索引
        allSugarRingAtoms: 所有糖环原子集合
        ringToUnit: atom_idx → unit_index 映射
        sourceUnitIdx: 出发侧的糖单元索引 (用来排除自身环)
        maxSteps: 最大遍历步数 (default=2, 覆盖 C6→C5 和更远的连接)

    Returns:
        (is_sugar_reachable, target_unit_idx)
    """
    # 快速检查: 起始碳本身就在另一个糖环中
    if startIdx in allSugarRingAtoms:
        targetUnit = ringToUnit[startIdx]
        if targetUnit != sourceUnitIdx:
            return True, targetUnit

    # BFS: 沿碳链遍历 (BFS along carbon chain)
    visited = {startIdx}
    frontier = [(startIdx, 0)]

    while frontier:
        currentIdx, depth = frontier.pop(0)
        if depth >= maxSteps:
            continue

        atom = mol.GetAtomWithIdx(currentIdx)
        for nbr in atom.GetNeighbors():
            nbrIdx = nbr.GetIdx()
            if nbrIdx in visited:
                continue
            # 只沿碳原子走 (Only walk along carbon atoms)
            if nbr.GetAtomicNum() != 6:
                continue

            visited.add(nbrIdx)

            # 检查是否到达另一个糖环 (Check if reached another sugar ring)
            if nbrIdx in allSugarRingAtoms:
                targetUnit = ringToUnit[nbrIdx]
                if targetUnit != sourceUnitIdx:
                    return True, targetUnit

            frontier.append((nbrIdx, depth + 1))

    return False, -1


def findGlycosidicBonds(
    mol: Chem.Mol,
    sugarUnits: List[Dict],
) -> List[Dict]:
    """
    精确定位所有糖苷键 (连接糖环 C1 与外部的 C-O/C-N/C-S 键)。
    Precisely locate all glycosidic bonds connecting sugar ring C1 to external moieties.

    Sugar-to-Sugar 判定 (Sugar-to-Sugar detection):
    键的另一端碳原子如果满足以下任一条件，则为 sugar_to_sugar:
    1. 直接属于另一个糖环 (e.g. C1-O-C4)
    2. 通过 ≤2 步 C-C 键可到达另一个糖环 (e.g. C1-O-C6-C5, C5 在环中)

    Args:
        mol: RDKit Mol 对象
        sugarUnits: find_mapped_sugar_units() 返回的糖单元列表

    Returns:
        糖苷键信息列表
    """
    # 收集所有糖环原子 (Collect all sugar ring atoms)
    allSugarRingAtoms: Set[int] = set()
    ringToUnit: Dict[int, int] = {}  # atom_idx → unit_index

    for unitIdx, unit in enumerate(sugarUnits):
        for atomIdx in unit["ring_atoms"]:
            allSugarRingAtoms.add(atomIdx)
            ringToUnit[atomIdx] = unitIdx

    glycosidicBonds = []

    for unitIdx, unit in enumerate(sugarUnits):
        anomericCarbons = findAnomericCarbons(mol, unit["ring_atoms"])

        for c1Idx in anomericCarbons:
            c1Atom = mol.GetAtomWithIdx(c1Idx)

            for nbr in c1Atom.GetNeighbors():
                nbrIdx = nbr.GetIdx()

                # 跳过环内原子 (Skip ring-internal atoms)
                if nbrIdx in set(unit["ring_atoms"]):
                    continue

                # 桥原子必须是 O/N/S (Bridge atom must be O/N/S)
                if nbr.GetAtomicNum() not in (8, 7, 16):
                    continue

                bridgeAtomIdx = nbrIdx
                bridgeAtom = nbr

                # 找到桥原子另一端连接的碳 (Find carbon on the other side of bridge)
                for farNbr in bridgeAtom.GetNeighbors():
                    if farNbr.GetIdx() == c1Idx:
                        continue
                    if farNbr.GetAtomicNum() != 6:
                        continue

                    farCarbonIdx = farNbr.GetIdx()

                    # 关键判定: farCarbon 能否通过短碳链到达另一个糖环?
                    # Key check: can farCarbon reach another sugar ring via short C-C chain?
                    isSugarReachable, targetUnit = _canReachSugarRing(
                        mol, farCarbonIdx, allSugarRingAtoms, ringToUnit,
                        sourceUnitIdx=unitIdx, maxSteps=2,
                    )

                    if isSugarReachable:
                        bondType = "sugar_to_sugar"
                    else:
                        targetUnit = -1
                        bondType = "sugar_to_aglycon"

                    # 找到要切断的键:
                    # 切 bridgeAtom-farCarbon 键 (保留桥原子在糖侧)
                    # Cut bridge-farCarbon bond (keep bridge atom with sugar)
                    bond = mol.GetBondBetweenAtoms(bridgeAtomIdx, farCarbonIdx)
                    if bond is None:
                        continue

                    glycosidicBonds.append({
                        "bond_idx": bond.GetIdx(),
                        "anomeric_carbon_idx": c1Idx,
                        "bridge_atom_idx": bridgeAtomIdx,
                        "far_carbon_idx": farCarbonIdx,
                        "source_unit": unitIdx,
                        "target_unit": targetUnit,
                        "bond_type": bondType,
                        "bridge_element": bridgeAtom.GetSymbol(),
                    })

    # 去重: 相同的 bond_idx 只保留一次 (Deduplicate by bond_idx)
    seen = set()
    uniqueBonds = []
    for b in glycosidicBonds:
        if b["bond_idx"] not in seen:
            seen.add(b["bond_idx"])
            uniqueBonds.append(b)

    return uniqueBonds


# =====================================================================
# 3. 精确切分 (Precision Cleavage)
# =====================================================================

def cleaveWithConservation(
    mol: Chem.Mol,
    sugarUnits: List[Dict],
) -> Tuple[str, str, Dict]:
    """
    精确切分糖苷键，保证碳原子守恒。
    Precision glycosidic bond cleavage with carbon atom conservation guarantee.

    Args:
        mol: RDKit Mol (原始完整分子)
        sugarUnits: find_mapped_sugar_units() 返回的糖单元列表

    Returns:
        (glycan_smiles, aglycon_smiles, metadata_dict)
        metadata_dict 包含: bonds_cut, carbon_check, linkage_info
    """
    metadata = {
        "bonds_cut": 0,
        "carbon_original": 0,
        "carbon_glycan": 0,
        "carbon_aglycon": 0,
        "carbon_conserved": False,
        "linkage_info": [],
    }

    if not sugarUnits:
        return "NULL", "NULL", metadata

    # Step 1: 碳原子总数 (Count original carbons)
    totalCarbon = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
    metadata["carbon_original"] = totalCarbon

    # Step 2: 定位糖苷键 (Locate glycosidic bonds)
    glycosidicBonds = findGlycosidicBonds(mol, sugarUnits)

    # 仅切断 sugar_to_aglycon 键 (Only cut sugar-to-aglycon bonds)
    # sugar_to_sugar 键保留 (Keep inter-sugar bonds intact)
    bondsToAglycon = [b for b in glycosidicBonds if b["bond_type"] == "sugar_to_aglycon"]

    if not bondsToAglycon:
        # 无苷元连接 → 整个分子就是糖链
        metadata["carbon_glycan"] = totalCarbon
        metadata["carbon_conserved"] = True
        return Chem.MolToSmiles(mol), "NULL", metadata

    # Step 3: 执行切分 (Execute cleavage)
    bondIndices = [b["bond_idx"] for b in bondsToAglycon]
    # 使用 isotope label 标记连接位点
    # Use isotope labels to mark connection sites
    # 糖侧标记 100+unit_idx, 苷元侧标记 200+unit_idx
    dummyLabels = []
    for b in bondsToAglycon:
        sugarLabel = 100 + b["source_unit"]
        aglyconLabel = 200 + b["source_unit"]
        # FragmentOnBonds: (begin_label, end_label) 对应原本的 beginAtom, endAtom
        bond = mol.GetBondWithIdx(b["bond_idx"])
        if bond.GetBeginAtomIdx() == b["bridge_atom_idx"]:
            # bridge 在 begin 侧 → begin 是糖侧
            dummyLabels.append((sugarLabel, aglyconLabel))
        else:
            dummyLabels.append((aglyconLabel, sugarLabel))

    fragMol = Chem.FragmentOnBonds(mol, bondIndices, dummyLabels=dummyLabels)

    # Step 4: 分离碎片 (Separate fragments)
    fragIndices = Chem.GetMolFrags(fragMol)

    # 收集所有真糖环原子 (Collect true sugar ring atoms)
    trueRingAtoms = set()
    for u in sugarUnits:
        trueRingAtoms.update(u["ring_atoms"])

    glycanFrags = []
    aglyconFrags = []

    for indices in fragIndices:
        isSugar = any(idx < mol.GetNumAtoms() and idx in trueRingAtoms for idx in indices)
        fragSmiles = Chem.MolFragmentToSmiles(fragMol, atomsToUse=list(indices), isomericSmiles=True)
        if isSugar:
            glycanFrags.append(fragSmiles)
        else:
            aglyconFrags.append(fragSmiles)

    glycanSmiles = ".".join(glycanFrags) if glycanFrags else "NULL"
    aglycanSmiles = ".".join(aglyconFrags) if aglyconFrags else "NULL"

    # Step 5: 碳原子守恒检查 (Carbon conservation check)
    glycanCarbon = 0
    aglycanCarbon = 0

    if glycanSmiles != "NULL":
        gMol = Chem.MolFromSmiles(glycanSmiles)
        if gMol:
            glycanCarbon = sum(1 for a in gMol.GetAtoms() if a.GetAtomicNum() == 6)

    if aglycanSmiles != "NULL":
        aMol = Chem.MolFromSmiles(aglycanSmiles)
        if aMol:
            aglycanCarbon = sum(1 for a in aMol.GetAtoms() if a.GetAtomicNum() == 6)

    metadata["carbon_glycan"] = glycanCarbon
    metadata["carbon_aglycon"] = aglycanCarbon
    metadata["carbon_conserved"] = (glycanCarbon + aglycanCarbon == totalCarbon)
    metadata["bonds_cut"] = len(bondIndices)

    # 记录连接信息 (Record linkage info)
    for b in bondsToAglycon:
        metadata["linkage_info"].append({
            "anomeric_C": b["anomeric_carbon_idx"],
            "bridge": f"{b['bridge_element']}({b['bridge_atom_idx']})",
            "far_C": b["far_carbon_idx"],
            "type": b["bond_type"],
        })

    if not metadata["carbon_conserved"]:
        print(f"[WARNING] Carbon conservation FAILED: "
              f"{totalCarbon} != {glycanCarbon} + {aglycanCarbon} "
              f"(lost {totalCarbon - glycanCarbon - aglycanCarbon})")

    return glycanSmiles, aglycanSmiles, metadata


if __name__ == "__main__":
    import sys
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
    from lib.glycan_topology import find_mapped_sugar_units

    RUTIN = "C[C@@H]1OC(OC[C@H]2OC(Oc3c(-c4ccc(O)c(O)c4)oc4cc(O)cc(O)c4c3=O)C(O)C(O)C2O)C(O)C(O)C1O"

    print("=" * 70)
    print("Precision Glycosidic Cleavage — Rutin Test")
    print("=" * 70)

    mol = Chem.MolFromSmiles(RUTIN)
    units = find_mapped_sugar_units(mol)
    glycan, aglycon, meta = cleaveWithConservation(mol, units)

    print(f"Glycan:  {glycan}")
    print(f"Aglycon: {aglycon}")
    print(f"Bonds cut: {meta['bonds_cut']}")
    print(f"Carbon: {meta['carbon_original']} = {meta['carbon_glycan']} + {meta['carbon_aglycon']}")
    print(f"Conserved: {meta['carbon_conserved']}")
    print(f"Linkage info: {meta['linkage_info']}")
