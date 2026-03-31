"""
[EN] RDKit Virtual Hydrolase & Monosaccharide Identity Engine
     Acts as the core analytical brain of the GlycoNP-Pipeline.
     Features:
      - Virtual structural cleavage (`RwMol`) to uncouple sugar units from aglycones without SMILES string corruption.
      - A comprehensive 3-Tier SMARTS stereochemical matching engine relying on over 120 reference templates.
      - Aggressive chirality checking falling back to rigorous generic labels (Hex, Pen) when explicitly 2D. 

[CN] RDKit 虚拟水解酶与单糖深度识别引擎
     整个 GlycoNP 管线的核心“鉴定大脑”。具备如下功能：
      - 利用 `RWMol` 构建虚拟分步水解操作，剥离糖环与其他结构并保持拓扑守恒。
      - 3 级 SMARTS 立体匹配模型：在超过 120+ 黄金模板上实施精确手性比对。
      - 当输入的 SMILES 确实剥夺了手性纵深时，严格拒绝利用骨架猜想，降级返回至科学准确的 Hex/Pen 集合。
"""
import sys
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import FindMolChiralCenters
from rdkit.Chem import rdFMCS
from rdkit.Chem import AllChem
from tqdm import tqdm
import copy
import networkx as nx
import numpy as np

# Dynamically add lib path to support imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

try:
    import glycan_topology
except ImportError:
    from lib import glycan_topology

# ---------- 1. Reference Library Definitions ----------
# 全量单糖 SMILES 库及预编译参考分子现已统一定义在 glycan_reference_library.py 中。
# All monosaccharide SMILES and precompiled reference molecules are now
# centralized in glycan_reference_library.py (single source of truth).
try:
    import glycan_reference_library
    from glycan_reference_library import (
        RAW_MONOSACCHARIDE_SMILES,
        REFERENCE_MOLS,
        SPECIFICITY_ORDER,
        REFERENCE_RING_SIZE,
    )
except ImportError:
    from lib import glycan_reference_library
    from lib.glycan_reference_library import (
        RAW_MONOSACCHARIDE_SMILES,
        REFERENCE_MOLS,
        SPECIFICITY_ORDER,
        REFERENCE_RING_SIZE,
    )

# ---------- 1b. CIP+Exo Fingerprint Engine (New Tier 1) ----------
# CIP+Exo 指纹引擎: 虚拟脱修饰 → CIP 重算 → 指纹匹配 → 标签回填
# Virtual demod → CIP reassign → fingerprint match → label reassembly
try:
    from virtual_demodify import virtualDemodify, COMPILED_PATTERNS as DEMOD_PATTERNS
    from cip_exo_engine import (
        extractSugarFingerprint, matchSugarFingerprint,
        getReferenceFingerprintDb, walkSugarRing,
    )
    CIP_EXO_AVAILABLE = True
except ImportError:
    try:
        from lib.virtual_demodify import virtualDemodify, COMPILED_PATTERNS as DEMOD_PATTERNS
        from lib.cip_exo_engine import (
            extractSugarFingerprint, matchSugarFingerprint,
            getReferenceFingerprintDb, walkSugarRing,
        )
        CIP_EXO_AVAILABLE = True
    except ImportError:
        CIP_EXO_AVAILABLE = False

# =====================================================================
# 原子计数氧门控 (Atom-Counting Oxygen Gate)
# 解决 L-Col 等低 OH 数脱氧糖劫持高 OH 数单糖的问题
# Prevents 3,6-dideoxy sugars (2 OH) from hijacking 6-deoxy (3 OH)
# or normal hexoses (4 OH)
# =====================================================================
def _countExocyclicHeteroatoms(mol, ringAtoms: set) -> int:
    """统一拓扑遍历: 计算糖环外的杂原子取代基数量。
    Unified topology walk: count exocyclic heteroatom substituents.

    规则 (Rules, v3.0 — 2026-03-18):
      层 0: 直接连接在环碳上的 O/N → 每个直接计 1
      层 1 (一跳 C-C): 环碳→外环碳(C6 etc)→O/N → 每个外环碳最多贡献 1

    设计意图 (Design Intent):
      严格 1-hop 半径 + 每外环碳限 1: 防止修饰基团 (如 OAc 的 C=O)
      导致过度计数。-C(=O)OH (羧基) 在这里只贡献 1, 与 -CH2OH 一致。
      Strict 1-hop radius + max-1-per-exo-C: prevents modification groups
      (e.g. OAc's C=O) from overcounting. -C(=O)OH counts as 1, same as -CH2OH.

    此函数同时用于 参考糖 和 实际碎片 的 O 计数, 确保 ref/frag 完全一致。
    Used for BOTH reference and fragment O-counting to guarantee parity.
    """
    oCount = 0
    ringSet = set(ringAtoms)
    for rIdx in ringSet:
        atom = mol.GetAtomWithIdx(rIdx)
        if atom.GetAtomicNum() != 6:  # 只看环上的碳 (Only ring carbons)
            continue
        for nbr in atom.GetNeighbors():
            nIdx = nbr.GetIdx()
            if nIdx in ringSet:
                continue
            # 层 0: 直接杂原子 O(8) 或 N(7)
            # Layer 0: Direct heteroatom on ring carbon
            if nbr.GetAtomicNum() in (8, 7):
                oCount += 1
            # 层 1: 一跳 C-C 链路 → 每个外环碳最多贡献 1 个杂原子
            # Layer 1: One-hop C-C → max 1 heteroatom per exo-Carbon
            elif nbr.GetAtomicNum() == 6:
                foundHetero = False
                for nnbr in nbr.GetNeighbors():
                    if nnbr.GetIdx() == rIdx:
                        continue  # 不回头看环碳
                    if nnbr.GetAtomicNum() in (7, 8) and not foundHetero:
                        oCount += 1
                        foundHetero = True  # 每个 exo-C 限 1
    return oCount


def _countRefOxygens(smiles: str) -> int:
    """计算参考 SMILES 的外环杂原子数 (使用统一拓扑遍历)。
    Count exocyclic heteroatom substituents in reference sugar SMILES.
    Uses the same topology walk as fragment counting for parity.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return -1
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if len(ring) in (5, 6):
            ringO = [i for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 8]
            if len(ringO) == 1:
                return _countExocyclicHeteroatoms(mol, set(ring))
    return -1


def _countFragmentRingOxygens(mol, ringAtoms: set) -> int:
    """计算糖碎片的外环杂原子数 (使用统一拓扑遍历)。
    Count exocyclic heteroatom substituents in actual fragment.
    Delegates to unified function for ref/frag parity.
    """
    return _countExocyclicHeteroatoms(mol, ringAtoms)


# 预计算每个参考糖的外环氧原子数
REFERENCE_OXYGEN_COUNTS = {
    k: _countRefOxygens(v)
    for k, v in RAW_MONOSACCHARIDE_SMILES.items()
}


# =====================================================================
# C/N 原子计数门控 (Carbon/Nitrogen Atom-Count Gate)
# 防止氨基己糖 (6C, 1N) 被错误匹配到戊糖 (5C, 0N)
# Prevents amino hexoses from matching pentose templates
# =====================================================================
def _countRefCN(smiles: str) -> tuple:
    """计算参考糖 SMILES 中的碳原子和氮原子总数。
    Count total C and N atoms in a reference sugar SMILES.
    Returns: (carbonCount, nitrogenCount)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (-1, -1)
    cCount = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
    nCount = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 7)
    return (cCount, nCount)


# 预计算每个参考糖的碳/氮原子数
# Precompute C/N counts for all reference sugars
REFERENCE_CN_COUNTS = {
    k: _countRefCN(v)
    for k, v in RAW_MONOSACCHARIDE_SMILES.items()
}


def _countFragmentCN(mol, sugarUnitAtoms: set) -> tuple:
    """计算糖碎片中的碳原子和氮原子总数 (含环内+外环)。
    Count total C and N in fragment (ring + exocyclic atoms within the sugar unit).
    Returns: (carbonCount, nitrogenCount)
    """
    cCount = 0
    nCount = 0
    for idx in sugarUnitAtoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6:
            cCount += 1
        elif atom.GetAtomicNum() == 7:
            nCount += 1
    return (cCount, nCount)

def isolate_sugar_ring(mol, ring_atoms):
    bonds_to_cut = []
    exo_atoms = set()
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in ring_atoms: exo_atoms.add(nbr.GetIdx())
                
    c6_idxs = []
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in ring_atoms and nbr.GetSymbol() == 'C':
                for nnbr in nbr.GetNeighbors():
                    if nnbr.GetSymbol() in ('O','N','S'):
                        c6_idxs.append(nbr.GetIdx())
                        exo_atoms.add(nnbr.GetIdx())
                        break
                        
    for exo_idx in list(exo_atoms):
        exo_atom = mol.GetAtomWithIdx(exo_idx)
        for nbr in exo_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx not in ring_atoms and nbr_idx not in c6_idxs:
                bond = mol.GetBondBetweenAtoms(exo_idx, nbr_idx)
                if bond: bonds_to_cut.append(bond.GetIdx())
                    
    if not bonds_to_cut: return Chem.Mol(mol), dict(zip(ring_atoms, ring_atoms))
        
    frag_mol = Chem.FragmentOnBonds(mol, bonds_to_cut)
    frag_mols = Chem.GetMolFrags(frag_mol, asMols=True)
    frag_indices = Chem.GetMolFrags(frag_mol, asMols=False)
    
    for frag, indices in zip(frag_mols, frag_indices):
        if not all(r in indices for r in ring_atoms): continue
        ri = frag.GetRingInfo()
        valid = any(len(r) in (5,6) and [frag.GetAtomWithIdx(a).GetSymbol() for a in r].count('O') == 1 for r in ri.AtomRings())
        if valid:
            rwmol = Chem.RWMol(frag)
            atoms_to_process = [(atom.GetIdx(), atom.GetNeighbors()[0].GetSymbol()) for atom in rwmol.GetAtoms() if atom.GetAtomicNum() == 0 and atom.GetNeighbors()]
            for dummy_idx, nbr_sym in atoms_to_process:
                dummy_atom = rwmol.GetAtomWithIdx(dummy_idx)
                dummy_atom.SetAtomicNum(8 if nbr_sym == 'C' else 1)
                dummy_atom.SetFormalCharge(0); dummy_atom.SetIsotope(0)
            
            rwmol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(rwmol)
            clean_mol = Chem.RemoveHs(rwmol.GetMol())
            Chem.AssignStereochemistry(clean_mol, force=True, cleanIt=True)
            
            # CRITICAL FIX: useChirality=False prevents CIP inversion from breaking the mapping!
            matches = mol.GetSubstructMatches(clean_mol, useChirality=False)
            best_match = next((match for match in matches if sum(1 for t in match if t in ring_atoms) >= len(ring_atoms)), None)
            
            if best_match:
                clean_to_target = {clean_idx: target_idx for clean_idx, target_idx in enumerate(best_match)}
                target_to_clean = {target_idx: clean_idx for clean_idx, target_idx in clean_to_target.items()}
                mapped_ring = [target_to_clean.get(old_r, -1) for old_r in ring_atoms]
                if -1 not in mapped_ring:
                    return clean_mol, dict(zip(ring_atoms, mapped_ring))
            
            new_ri = clean_mol.GetRingInfo()
            if new_ri.AtomRings(): return clean_mol, dict(zip(ring_atoms, new_ri.AtomRings()[0]))
    return None, None

def check_modifications(mol, ring_atoms, ref_smiles=None):
    mods = []
    is_acid = False
    is_complex = False
    
    try:
        ref_mol = Chem.MolFromSmiles(ref_smiles) if ref_smiles else None
        
        if ref_smiles and ("A" in ref_smiles or "acid" in ref_smiles.lower()):
            is_acid = True
            
        def traverse(start_idx, avoid_atoms):
            visited = set([start_idx])
            q = [start_idx]
            while q:
                curr = q.pop(0)
                for nbr in mol.GetAtomWithIdx(curr).GetNeighbors():
                    idx = nbr.GetIdx()
                    if idx not in avoid_atoms and idx not in visited:
                        visited.add(idx)
                        q.append(idx)
            return visited

        for r_idx in ring_atoms:
            r_atom = mol.GetAtomWithIdx(r_idx)
            for nbr in r_atom.GetNeighbors():
                if nbr.GetIdx() not in ring_atoms:
                    branch = traverse(nbr.GetIdx(), set(ring_atoms))
                    if len(branch) > 15:
                        continue
                        
                    try:
                        sub_smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=list(branch))
                        if sub_smiles:
                            # 识别 NAc/Ac/NH2
                            if any(x in sub_smiles for x in ['NC(C)=O', 'CC(N)=O', 'CC(=O)N']): mods.append('NAc')
                            elif 'C(=O)C' in sub_smiles or sub_smiles == 'CC(=O)': mods.append('Ac')
                            elif sub_smiles == 'N': mods.append('N')
                            
                            # 识别脱氧 (直接连在环上的单碳)
                            elif sub_smiles == 'C': mods.append('deoxy')
                            
                            # 识别 O-甲基 (O-Me)
                            elif sub_smiles in ['CO', 'OC']:
                                # 注意: C6 的正常羟甲基 (CH2OH) 的 sub_smiles 也是 CO
                                # 我们之前在遍历字典时已经识别了基底(比如 Glc)，此时它不应该被当作修饰
                                # 真正的 O-Me 修饰通常在这个函数的更复杂逻辑中，或被基底忽略。
                                # 为了安全起见保留原样。
                                pass
                            elif 'S(=O)(=O)' in sub_smiles: mods.append('S')
                            
                            # 识别糖醛酸 (O=CO or C(=O)O)
                            elif 'O=CO' in sub_smiles or 'C(=O)O' in sub_smiles:
                                mods.append('A')
                                is_acid = True
                                
                            elif 'P(=O)' in sub_smiles: mods.append('P')
                            elif len(branch) >= 6: is_complex = True
                    except:
                        pass
                        
    except:
        pass
        
    final_mods = sorted(list(set(mods)))
    return final_mods, is_acid, is_complex

def create_virtual_standard_sugar(mol, ring_atoms):
    """
    创建分子的临时副本，用于 SubstructMatch 模板匹配。
    Create a virtual copy of the molecule for SubstructMatch template matching.

    1. N/S → O 突变: 使氨基糖匹配基础模板 (N/S mutation for amino sugars)
    2. 糖苷键切断 + OH 封端: 隔离糖片段 (Glycosidic bond cleavage + OH cap)
       - 设计意图: 大型苷元 (大环内酯/甾体/多萜) 会使 SubstructMatch 失败
       - 因为异头碳的 C-O-C 糖苷键不匹配模板的 C-OH (异头 OH)
       - Design: Large aglycones (macrolides/steroids/terpenes) cause SubstructMatch
         to fail because the anomeric C-O-C glycosidic bond doesn't match the
         template's C-OH (anomeric hydroxyl)
       - 此处将糖苷键 O 远端的苷元碳切断, 只保留 O 作为自由 OH
       - We cleave the aglycone carbon beyond the glycosidic O, preserving
         the O as a free hydroxyl

    注意: 脱氧位置（仅连有 H 或 C）不受影响，完美保护 Rhamnose 等脱氧糖。
    Note: Deoxy positions (bonded to H/C only) are unaffected, preserving deoxy sugars.
    """
    rwmol = Chem.RWMol(mol)
    ringSet = set(ring_atoms)

    # Step 1: N/S → O 突变 (原有逻辑)
    # Step 1: N/S → O mutation (original logic)
    for idx in ring_atoms:
        atom = rwmol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6:
            for nbr in atom.GetNeighbors():
                nbrIdx = nbr.GetIdx()
                if nbrIdx not in ringSet:
                    if nbr.GetAtomicNum() in [7, 16]:
                        nbr.SetAtomicNum(8)

    # Step 2: 糖苷键切断 — 隔离糖片段
    # Step 2: Glycosidic bond cleavage — isolate sugar fragment
    # 策略: 找到所有 "糖环 C → 环外 O → 非环 C" 链, 断开 O-苷元C 键
    # Strategy: find all "ring-C → exo-O → non-ring-C" chains, break the O-aglycone bond
    # 这样异头碳的 C-O-aglycone 变成 C-O-H (自由 OH), 原子索引不变
    # The anomeric C-O-aglycone becomes C-O-H (free OH), atom indices unchanged
    bondsToBreak = []
    for idx in ring_atoms:
        atom = rwmol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:
            continue
        for nbr in atom.GetNeighbors():
            nbrIdx = nbr.GetIdx()
            if nbrIdx in ringSet or nbr.GetAtomicNum() != 8:
                continue
            # 找到环外 O, 检查它是否连接到非环 C (糖苷键)
            for nbr2 in nbr.GetNeighbors():
                nbr2Idx = nbr2.GetIdx()
                if nbr2Idx == idx or nbr2Idx in ringSet:
                    continue
                if nbr2.GetAtomicNum() == 6:
                    # 苷元判定: 远端 C 是否有 ≥2 个重原子邻居 (不含桥连 O)?
                    # Aglycone test: does distal C have ≥2 heavy neighbors (excl. bridging O)?
                    # 如果只有 1 个重原子邻居 (= 桥连 O), 则是 -OCH3 (OMe), 不切
                    # If only 1 heavy neighbor (= bridging O), it's -OCH3 (OMe), don't cut
                    heavyNbrs = sum(1 for n in nbr2.GetNeighbors()
                                    if n.GetAtomicNum() > 1 and n.GetIdx() != nbrIdx)
                    if heavyNbrs >= 1:
                        # 远端 C 连接更大的结构 → 苷元, 断键
                        bondsToBreak.append((nbrIdx, nbr2Idx))

    for oIdx, cIdx in bondsToBreak:
        bondIdx = rwmol.GetBondBetweenAtoms(oIdx, cIdx)
        if bondIdx is not None:
            rwmol.RemoveBond(oIdx, cIdx)

    try:
        return rwmol.GetMol()
    except Exception:
        return Chem.RWMol(mol).GetMol()


def rescue_hexose_by_key_nodes(mol, ring_atoms):
    """CIP 立体化学救援引擎 v2.0 — 强化版 (Hardened CIP Rescue Engine v2.0)

    当 Tier 1-3 子结构匹配全部失败时启动。依靠图遍历读取 C2/C3/C4 的
    CIP (R/S) 构型进行靶向救援。

    v2.0 强化 (2026-03-19):
      1. C3 交叉验证: D-Glc/Gal/Man 的 C3 均应为 S, 若 C3=R → 置信度降级
      2. 脱氧糖扩展: 检测 O 计数偏低时推断 L-Rha/L-Fuc (天然产物高频糖)
      3. 置信度门控: 需 C2+C4 均有定义的 CIP 才执行 (否则直接 Hex)
      4. 修饰复杂度门控: 环外重原子 >25 时不信任 CIP (修饰可能翻转优先级)
      5. 置信度标注: (CIP:high) = C3 验证通过; (CIP:low) = C3 缺失或不一致

    Hardened CIP Rescue v2.0 improvements:
      1. C3 cross-validation: All three common hexoses (Glc/Gal/Man) have C3=S
      2. Deoxy sugar extension: Low O-count → L-Rha/L-Fuc candidates
      3. Confidence gate: Both C2 AND C4 must have defined CIP codes
      4. Modification complexity gate: >25 exocyclic heavy atoms → untrusted
      5. Confidence tag: (CIP:high) when C3 validated; (CIP:low) otherwise

    化学依据 (Chemical Basis):
      - D-Glc: C2=S, C3=S, C4=S (最常见己糖)
      - D-Gal: C2=S, C3=S, C4=R (第二常见)
      - D-Man: C2=R, C3=S, C4=S (第三常见)
      - L-Rha (6-deoxy-L-Man): C2=S, C3=R, C4=R + 低O计数
      - L-Fuc (6-deoxy-L-Gal): C2=S, C3=R, C4=S + 低O计数

    参考文献 (References):
      - Varki et al., Essentials of Glycobiology (4th ed.), Ch.2
      - McNaught, Nomenclature of Carbohydrates, Adv. Carbohydr. Chem. (1997)
    """
    if len(ring_atoms) != 6:
        return "Pen", "?"  # 仅救援六元糖环，五元保持 (pyranose only)

    # 强制 RDKit 计算分子的手性构型
    # Force RDKit to compute stereochemistry assignments
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    try:
        ring_set = set(ring_atoms)

        # ── Gate 1: 修饰复杂度门控 ──────────────────────────────
        # 环外重原子过多 → 修饰基团可能干扰 CIP 优先级排序
        # Too many exocyclic heavy atoms → modifications may invert CIP priorities
        exocyclicHeavy = 0
        for rIdx in ring_set:
            for nbr in mol.GetAtomWithIdx(rIdx).GetNeighbors():
                if nbr.GetIdx() not in ring_set and nbr.GetAtomicNum() > 1:
                    exocyclicHeavy += 1
        if exocyclicHeavy > 25:
            return "Hex", "?"

        # ── Step 1: 定位环氧 → C1 → C2 → C3 → C4 → C5 ─────────
        ring_oxygens = [idx for idx in ring_atoms
                        if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8]
        if not ring_oxygens:
            return "Hex", "?"
        ring_o_idx = ring_oxygens[0]

        # C1: 异头碳 (连环O + 环外杂原子)
        # Bug 6 修复: 增加 C-糖苷回退 (C-glycoside fallback)
        # C-糖苷 (如 Vitexin) 的 C1 连接环氧 + 外环碳, 不是外环氧
        # C-glycosides have C1 bonded to ring-O + exo-C, NOT exo-O
        c1_idx = None
        for nbr in mol.GetAtomWithIdx(ring_o_idx).GetNeighbors():
            if nbr.GetIdx() in ring_set:
                heteroCount = sum(1 for n in nbr.GetNeighbors()
                                 if n.GetAtomicNum() in (8, 7, 16))
                if heteroCount >= 2:
                    c1_idx = nbr.GetIdx()
                    break
        # Bug 6 回退: C-糖苷 C1 定位 — 选择与环氧相邻的碳中,
        # 外环邻居数最多的那个 (C1 通常连接苷元/OH/另一个糖)
        # Fallback for C-glycosides: choose the ring-O neighbor carbon
        # with the most non-ring neighbors (C1 connects to aglycon/OH/sugar)
        if c1_idx is None:
            candidates = []
            for nbr in mol.GetAtomWithIdx(ring_o_idx).GetNeighbors():
                if nbr.GetIdx() in ring_set and nbr.GetAtomicNum() == 6:
                    exoNbrCount = sum(1 for n in nbr.GetNeighbors()
                                     if n.GetIdx() not in ring_set)
                    # Bug 6 tiebreaker: C-糖苷的 C1 连接芳香性苷元碳
                    # 优先选择有芳香外环邻居的碳作为 C1
                    # C-glycoside tiebreaker: C1 bonds to aromatic aglycon C
                    hasAromaticExo = any(
                        n.GetIsAromatic() for n in nbr.GetNeighbors()
                        if n.GetIdx() not in ring_set
                    )
                    candidates.append((exoNbrCount, hasAromaticExo, nbr.GetIdx()))
            if candidates:
                # 排序: exo数降序 → 芳香性优先 (True > False) → 索引
                candidates.sort(key=lambda x: (-x[0], -x[1], x[2]))
                c1_idx = candidates[0][2]
        if c1_idx is None:
            return "Hex", "?"

        # 沿碳链遍历 C2→C3→C4→C5
        c2_idx = [n.GetIdx() for n in mol.GetAtomWithIdx(c1_idx).GetNeighbors()
                  if n.GetIdx() in ring_set and n.GetIdx() != ring_o_idx][0]
        c3_idx = [n.GetIdx() for n in mol.GetAtomWithIdx(c2_idx).GetNeighbors()
                  if n.GetIdx() in ring_set and n.GetIdx() != c1_idx][0]
        c4_idx = [n.GetIdx() for n in mol.GetAtomWithIdx(c3_idx).GetNeighbors()
                  if n.GetIdx() in ring_set and n.GetIdx() != c2_idx][0]
        c5_idx = [n.GetIdx() for n in mol.GetAtomWithIdx(c4_idx).GetNeighbors()
                  if n.GetIdx() in ring_set and n.GetIdx() != c3_idx][0]

        # ── Step 2: 读取 CIP 标签 ──────────────────────────────
        def _getCip(idx):
            a = mol.GetAtomWithIdx(idx)
            return a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else "?"

        c2Cip = _getCip(c2_idx)
        c3Cip = _getCip(c3_idx)
        c4Cip = _getCip(c4_idx)

        # ── Gate 2: 置信度门控 — C2+C4 必须均有定义 ────────────
        # Confidence gate: BOTH C2 AND C4 must have defined CIP
        if c2Cip == "?" or c4Cip == "?":
            return "Hex", "?"

        # ── Step 3: 检测氨基糖 ──────────────────────────────────
        hasNitrogenOnC2 = any(
            n.GetAtomicNum() == 7
            for n in mol.GetAtomWithIdx(c2_idx).GetNeighbors()
            if n.GetIdx() not in ring_set
        )
        hasAnyNitrogen = any(
            mol.GetAtomWithIdx(rIdx).GetAtomicNum() == 7
            or any(n.GetAtomicNum() == 7
                   for n in mol.GetAtomWithIdx(rIdx).GetNeighbors()
                   if n.GetIdx() not in ring_set)
            for rIdx in ring_set if mol.GetAtomWithIdx(rIdx).GetAtomicNum() == 6
        )
        aminoSuffix = "N" if (hasNitrogenOnC2 or hasAnyNitrogen) else ""

        # ── Step 4: 检测脱氧糖特征 ─────────────────────────────
        # 脱氧糖: C5 (或 C6 外环碳) 上仅连 C/H, 不连 O
        # Deoxy detection: C5 exocyclic neighbor is CH3, not CH2OH
        isDeoxy = False
        c5Atom = mol.GetAtomWithIdx(c5_idx)
        c5ExoO = sum(1 for n in c5Atom.GetNeighbors()
                     if n.GetIdx() not in ring_set and n.GetAtomicNum() == 8)
        c5ExoC = [n for n in c5Atom.GetNeighbors()
                  if n.GetIdx() not in ring_set and n.GetAtomicNum() == 6]
        if c5ExoO == 0 and len(c5ExoC) == 1:
            # C5 的外环碳上没有 O → 可能是 6-deoxy (CH3 而非 CH2OH)
            exoC = c5ExoC[0]
            exoCO = sum(1 for n in exoC.GetNeighbors()
                        if n.GetAtomicNum() == 8)
            if exoCO == 0:
                isDeoxy = True  # 确认: 外环碳无O = 6-deoxysugar

        # ── Step 5: C3 交叉验证 + 置信度判定 ────────────────────
        # D-Glc/Gal/Man 在 Fischer 投影中 C3-OH 均为右旋 → CIP=S
        # L-Rha/Fuc 为 L-构型 → C3-OH 左旋 → CIP=R
        # C3 cross-validation: common D-hexoses have C3=S; L-deoxysugars have C3=R
        c3Validated = False
        c3Confidence = "low"

        # ── Step 6: CIP 靶向判决 ────────────────────────────────
        # 仅输出常见糖, 禁止猜测罕见糖 (D-Tal, D-All 等)
        # COMMON sugars only. Rare sugars must earn labels via Tier-1.

        # ---- 6a. 脱氧糖优先检测 (Deoxy sugars first) ----
        if isDeoxy:
            # L-Rha (6-deoxy-L-Man): CIP 特征 C4=R, C2=S
            # L-Fuc (6-deoxy-L-Gal): CIP 特征 C4=S, C2=S
            # 注意: L-构型糖的 CIP 因 C5 构型反转而与 D-糖不同
            if c4Cip == "R" and c2Cip == "S":
                c3Confidence = "high" if c3Cip == "R" else "low"
                tag = f"(CIP:{c3Confidence})"
                return f"L-Rha{aminoSuffix}{tag}", "?"
            elif c4Cip == "S" and c2Cip == "S":
                c3Confidence = "high" if c3Cip == "R" else "low"
                tag = f"(CIP:{c3Confidence})"
                return f"L-Fuc{aminoSuffix}{tag}", "?"
            elif c4Cip == "S" and c2Cip == "R":
                c3Confidence = "high" if c3Cip == "R" else "low"
                tag = f"(CIP:{c3Confidence})"
                return f"D-Rha{aminoSuffix}{tag}", "?"
            else:
                # 其他脱氧糖组合 → dHex 泛标签
                return f"dHex{aminoSuffix}", "?"

        # ---- 6b. 标准己糖判决 (Standard hexoses) ----
        if c4Cip == "R" and c2Cip == "S":
            # D-Gal 特征: C2=S, C3=S, C4=R
            c3Confidence = "high" if c3Cip == "S" else "low"
            tag = f"(CIP:{c3Confidence})"
            return f"D-Gal{aminoSuffix}{tag}", "?"
        elif c4Cip == "S" and c2Cip == "S":
            # D-Glc 特征: C2=S, C3=S, C4=S
            c3Confidence = "high" if c3Cip == "S" else "low"
            tag = f"(CIP:{c3Confidence})"
            return f"D-Glc{aminoSuffix}{tag}", "?"
        elif c4Cip == "S" and c2Cip == "R":
            # D-Man 特征: C2=R, C3=S, C4=S
            c3Confidence = "high" if c3Cip == "S" else "low"
            tag = f"(CIP:{c3Confidence})"
            return f"D-Man{aminoSuffix}{tag}", "?"
        elif c4Cip == "R" and c2Cip == "R":
            # C4=R + C2=R → 指向 D-Tal (极罕见), 不猜测 → Hex
            return f"Hex{aminoSuffix}", "?"

        return "Hex", "?"
    except Exception:
        # 图遍历结构异常 → 安全退回 Hex (Safe fallback)
        return "Hex", "?"

def _countFragmentCarbons(mol, ring_atoms, maxHops: int = 4):
    """
    杂原子截断铁律: 严格计算糖的 C-C 连续骨架碳数。
    Heteroatom Blockade: count sugar backbone carbons via C-C bonds ONLY.

    化学定义: 糖的碳数分类 (Hex=6C, Hept=7C, Oct=8C, Non=9C) 仅由
    连续的碳-碳 (C-C) 骨架决定。杂原子 (O, N, S) 是硬截断屏障。
    Chemical rule: Sugar carbon classification is determined ONLY by the
    contiguous C-C backbone. Heteroatoms (O, N, S) are hard barriers.

    v2.1 修复: 增加 BFS 步数限制 (maxHops=4) 和分支度门控。
    v2.1 fix: Added BFS depth limit and branching-degree gate.
    - maxHops=4 覆盖: 环碳 → C6 → C7 → C8 → C9 (覆盖壬糖/唾液酸)
    - 分支度门控: 非环碳若有 ≥3 个碳邻居 → 已进入苷元骨架, 停止
    - maxHops=4 covers: ring-C → C6 → C7 → C8 → C9 (covers nonoses)
    - Branch gate: non-ring carbon with ≥3 carbon neighbors → in aglycon, stop

    示例 (Examples):
      - C-O-C(=O)CH3 (乙酰化): O 处截断, 外侧 2C 不计入
      - C-NH-C(=O)CH3 (NAc): N 处截断, 外侧 2C 不计入
      - C-C(OH)-C(OH)-CH2OH (甘油侧链): 沿 C-C 遍历, 全部计入 (≤4 步)
    """
    ringSet = set(ring_atoms)
    visited = set(ring_atoms)   # 环原子全部预标记已访问
    totalCarbons = 0

    # 首先计算环上的碳原子数
    carbonQueue = []
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6:
            totalCarbons += 1
            carbonQueue.append((idx, 0))  # (atomIdx, depth)
        # 环氧标记为 visited 但不扩展 (它是种子边界)

    # 从环上碳原子出发, 仅沿 C-C 键向外 BFS (带步数限制)
    # BFS from ring carbons, follow C-C bonds only (with depth limit)
    while carbonQueue:
        idx, depth = carbonQueue.pop(0)
        if depth >= maxHops:
            continue  # 步数限制: 不再扩展 (depth limit: stop expanding)
        for nbr in mol.GetAtomWithIdx(idx).GetNeighbors():
            nIdx = nbr.GetIdx()
            if nIdx in visited:
                continue
            visited.add(nIdx)

            # 杂原子截断铁律: O(8), N(7), S(16), P(15) → 硬停
            if nbr.GetAtomicNum() != 6:
                continue  # 遇到杂原子, 不再沿此分支前进

            # 分支度门控: 非环碳若有 ≥3 个碳邻居 → 已进入骨架, 停止
            # Branching gate: non-ring carbon with ≥3 C neighbors → in backbone
            if nIdx not in ringSet:
                carbonNbrCount = sum(
                    1 for n in nbr.GetNeighbors() if n.GetAtomicNum() == 6
                )
                if carbonNbrCount >= 3:
                    continue  # 进入了大型碳骨架, 不计入糖碳

            # 仅碳原子 → 计数并继续扩展
            totalCarbons += 1
            carbonQueue.append((nIdx, depth + 1))

    return totalCarbons


def _buildAnomericAmnestyRef(refMol):
    """法则 3 辅助: 去除参考糖异头碳的手性标记, 构建'宽恕版'SMARTS.
    Golden Rule 3 helper: remove chirality from the anomeric carbon in ref mol.
    异头碳定义: 环内碳, 同时连接环内氧和环外氧 (C bonded to ring-O and exo-O).
    Anomeric carbon: ring C bonded to both ring-O and exocyclic-O.
    """
    try:
        rwMol = Chem.RWMol(Chem.MolFromSmiles(Chem.MolToSmiles(refMol)))
        if rwMol is None:
            return None
        ri = rwMol.GetRingInfo()
        ringAtomSets = [set(r) for r in ri.AtomRings()]
        for ringSet in ringAtomSets:
            for idx in ringSet:
                atom = rwMol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() != 6:
                    continue
                # 检查是否连接环内 O 和环外 O
                hasRingO = False
                hasExoO = False
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 8 and nbr.GetIdx() in ringSet:
                        hasRingO = True
                    if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ringSet:
                        hasExoO = True
                if hasRingO and hasExoO:
                    atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        smartsStr = Chem.MolToSmarts(rwMol.GetMol())
        return Chem.MolFromSmarts(smartsStr)
    except Exception:
        return None


# 预编译异头碳宽恕版参考分子 (Pre-build anomeric amnesty refs)
AMNESTY_REFERENCE_MOLS = {}
for k, refMol in REFERENCE_MOLS.items():
    amnestyMol = _buildAnomericAmnestyRef(refMol)
    if amnestyMol is not None:
        AMNESTY_REFERENCE_MOLS[k] = amnestyMol


def identify_monosaccharide_v2(mol, ring_atoms):
    """核心匹配引擎 v5.0 — 罕见糖隔离 + 模糊手性打分
    Core Matching Engine v5.0 — Rare Sugar Quarantine + Fuzzy Chiral Scoring

    Tier 1: 绝对手性匹配 (所有糖, useChirality=True)
    Tier 2: 异头碳手性特赦 (仅 COMMON_HEXOSES, 罕见糖隔离)
    Tier 3: 模糊手性打分 (仅 COMMON_HEXOSES, useChirality=False + 手性重合度)

    Tier 1: Strict chirality (all sugars)
    Tier 2: Anomeric amnesty (COMMON only, rare sugars quarantined)
    Tier 3: Fuzzy chiral scoring (COMMON only, useChirality=False + chiral overlap)
    """
    fragmentRingSize = len(ring_atoms)
    base_name = "Hex" if fragmentRingSize == 6 else "Pen"
    matched_anomer = "?"

    # 虚拟分身: N→O 突变使氨基糖能匹配基础己糖模板
    # Virtual clone: N→O mutation enables amino sugars to match base hexose templates
    virtual_mol = create_virtual_standard_sugar(mol, ring_atoms)

    # 安全碳/氮计数 — 从 virtual_mol 计算, 与 SubstructMatch 目标一致
    # C/N counts from virtual_mol for consistency with SubstructMatch target
    # 修复 (2026-03-18): 旧代码从原始 mol 计算 fragN, 导致 GlcNAc 的
    # N-gate (fragN=1, refN=0) 阻止其匹配 D-Glc 模板。
    # Fix: compute from virtual_mol where N→O, so fragN=0 for amino sugars,
    # allowing GlcNAc virtual mol to match D-Glc template correctly.
    expandedAtoms = set(ring_atoms)
    for rIdx in ring_atoms:
        for nbr in virtual_mol.GetAtomWithIdx(rIdx).GetNeighbors():
            expandedAtoms.add(nbr.GetIdx())

    fragC = sum(1 for idx in expandedAtoms if virtual_mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    fragN = sum(1 for idx in expandedAtoms if virtual_mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)

    # O-Gate 预计算 — 同样从 virtual_mol 计算
    fragmentOxygens = _countFragmentRingOxygens(virtual_mol, set(ring_atoms))

    # ================================================================
    # 三层匹配 (Three-Tier Matching)
    # ================================================================

    # --- Tier 1: 绝对手性匹配 (所有糖参与) ---
    base_name, matched_anomer = _tierMatch(
        virtual_mol, mol, ring_atoms, fragmentRingSize, fragC, fragN, fragmentOxygens,
        REFERENCE_MOLS, useChirality=True, allowRare=True
    )
    if matched_anomer != "?":
        return _applyFuzzyAndMods(mol, ring_atoms, base_name, matched_anomer)

    # --- Tier 2: 异头碳特赦 (仅 COMMON, 罕见糖被隔离!) ---
    base_name, matched_anomer = _tierMatch(
        virtual_mol, mol, ring_atoms, fragmentRingSize, fragC, fragN, fragmentOxygens,
        AMNESTY_REFERENCE_MOLS, useChirality=True, allowRare=False
    )
    if matched_anomer != "?":
        return _applyFuzzyAndMods(mol, ring_atoms, base_name, matched_anomer)

    # --- Tier 3: 模糊手性打分 (仅 COMMON, useChirality=False + 手性重合度) ---
    base_name, matched_anomer = _tier3FuzzyChiralMatch(
        mol, virtual_mol, ring_atoms, fragmentRingSize, fragC, fragN, fragmentOxygens
    )
    if matched_anomer != "?":
        return _applyFuzzyAndMods(mol, ring_atoms, base_name, matched_anomer)

    # 全部失败 → CIP 救援
    base_name = "Hex" if fragmentRingSize == 6 else "Pen"
    matched_anomer = "?"
    base_name, matched_anomer = rescue_hexose_by_key_nodes(mol, ring_atoms)

    return _applyFuzzyAndMods(mol, ring_atoms, base_name, matched_anomer)


# ================================================================
# 罕见糖隔离协议 (Rare Sugar Quarantine Protocol)
# 天然产物中这些糖极其罕见, 必须 Tier 1 绝对手性全匹配才给标签
# ================================================================
COMMON_HEXOSES = {
    "D-Glc", "L-Glc", "D-Gal", "L-Gal", "D-Man", "L-Man",
    "L-Rha", "D-Rha", "L-Fuc", "D-Fuc", "D-Qui", "L-Qui",
}
RARE_HEXOSES = {
    "D-Tal", "D-All", "D-Alt", "D-Gul", "D-Ido", "L-Ido",
    "L-Tal", "L-All", "L-Alt", "L-Gul",
}


# ================================================================
# 贝叶斯先验权重 (Bayesian Prior Weights)
# 天然产物中的糖频率先验: Tier-S (绝对高频) >> Tier-C (罕见)
# ================================================================
BAYESIAN_PRIOR_WEIGHT = {}

# Tier-X: 特种精确氨基糖组 (权重 150)
# 当这些糖精确匹配自身特定配置 (含 N 位置等) 时, 必须压倒后续将 N->O 突变后匹配的由于字母排序优先的普通己糖 (如 D-Glc)。
for s in ("D-Kanosamine", "D-6aGlc"):
    BAYESIAN_PRIOR_WEIGHT[s] = 150

# Tier-S: 绝对高频组 (权重 100)
for s in ("D-Glc", "L-Glc", "D-Gal", "L-Gal", "D-Man", "L-Man",
          "L-Rha", "D-Rha", "L-Fuc", "D-Fuc",
          "D-Xyl", "L-Ara", "D-Ara", "D-Rib",
          "D-GlcA", "D-GalA", "D-ManA", "D-GulA", "D-IdoA",
          "D-GlcN", "D-GalN", "D-ManN",
          "D-GlcNAc", "D-GalNAc", "D-ManNAc",
          "D-Qui", "L-Qui", "D-Fru",
          "Neu5Ac", "Neu5Gc", "KDO", "D-Api"):
    BAYESIAN_PRIOR_WEIGHT[s] = 100
# Tier-C: 罕见组 (权重 1)
for s in ("D-Tal", "L-Tal", "D-All", "L-All", "D-Alt", "L-Alt",
          "D-Gul", "L-Gul", "D-Ido", "L-Ido",
          "D-Lyx", "L-Lyx",
          "D-6dTal", "D-6dGul", "D-6dAlt"):
    BAYESIAN_PRIOR_WEIGHT[s] = 1
# Tier-B: 特种脱氧+心苷糖 (权重 10, 低于标准糖但高于罕见)
for s in ("D-dGlc", "L-dGlc", "D-dGal", "L-dGal", "D-dMan",
          "D-Dig", "D-Boi", "L-Ole", "D-Cym", "L-The",
          "D-Dtx", "L-Dtx", "D-Sar", "D-Cma", "D-Din",
          "L-Ola", "D-Eva", "L-Eva"):
    BAYESIAN_PRIOR_WEIGHT[s] = 10


def _hasTargetChirality(mol, ring_atoms):
    """检测目标糖环是否有指定手性 (Detect if target sugar ring has chirality).

    如果环碳全部没有手性标记 (ChiralTag == CHI_UNSPECIFIED),
    说明原始 SMILES 是平面的, useChirality=True 必然失败。
    If all ring carbons have no chiral tag, the original SMILES is flat
    and useChirality=True will always fail for any chiral query.

    Returns:
        bool: True if at least one ring carbon has specified chirality.
    """
    ringSet = set(ring_atoms)
    for idx in ringSet:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6:
            if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                return True
    return False


def _tierMatch(virtual_mol, original_mol, ring_atoms, fragmentRingSize, fragC, fragN,
               fragmentOxygens, refDict, useChirality, allowRare):
    """全量竞争匹配层 v7.0: 收集所有合法匹配 + 贝叶斯先验裁决 + N-模板双轨匹配
    Full Competition Matching v7.0: All valid matches + Bayesian + dual-track N-template.

    v7.0 新增: 对含 N 的参考模板 (Neu5Ac, Neu5Gc), 使用 original_mol
    (未经 N→O 突变) 进行 SubstructMatch, 解决 N→O 突变导致的匹配失败。
    """
    validMatches = []  # list of (name, anomer, weight)

    for key in SPECIFICITY_ORDER:
        refMol = refDict.get(key)
        if refMol is None:
            continue
        name, anomer = key

        # 罕见糖隔离
        if not allowRare and name in RARE_HEXOSES:
            continue

        # 环大小门控
        refRingSize = REFERENCE_RING_SIZE.get(key, 6)
        if refRingSize != fragmentRingSize:
            continue

        # 碳数门控: 标准糖严格相等, 多碳糖放宽
        # Carbon gate: strict equality for standard sugars,
        # relaxed for multi-C sugars (Neu5Ac/KDO)
        refC, refN = REFERENCE_CN_COUNTS.get(key, (-1, -1))
        if refC >= 0:
            if refC <= 7:
                if fragC != refC:
                    continue
            else:
                if fragC > refC + 1:
                    continue
            # N-Gate: 对含 N 的参考糖 (Neu5Ac/Neu5Gc), 跳过 N 互斥检查
            # N-Gate: skip mutual exclusion for N-containing templates
            if refN == 0:
                if fragN > 0:
                    continue
            # refN > 0 时不检查 fragN (因为 virtual_mol 已将 N→O, fragN=0)

        # 子结构匹配 (Substructure Matching)
        # 异头碳手性已在字典加载阶段剥离 (_findAnomericCarbon),
        # 但若目标分子完全无手性标记 (缺少 @/@@ 的平面 SMILES),
        # useChirality=True 会因"查询有手性、目标没有"而必然失败。
        # 修复: 检测到无手性目标时, 自动降级为 useChirality=False。
        #
        # Achiral Target Fallback: If the target fragment has ZERO
        # specified stereocenters and useChirality=True is requested,
        # the match ALWAYS fails because query has @/@@ but target
        # has none. Fix: detect achiral targets and fallback to
        # useChirality=False to prevent false Hex results.
        # 选择匹配分子: 含 N 的参考模板用原始分子, 否则用 N→O 突变体
        # Dual-track: N-containing templates match against original_mol,
        # standard templates match against virtual_mol (N→O mutant)
        matchMol = original_mol if refN > 0 else virtual_mol
        effectiveChiral = useChirality
        if useChirality and not _hasTargetChirality(matchMol, ring_atoms):
            effectiveChiral = False
        matches = matchMol.GetSubstructMatches(refMol, useChirality=effectiveChiral)
        matched = False
        for match in matches:
            if not all(idx in match for idx in ring_atoms):
                continue

            # O-Gate
            refOCount = REFERENCE_OXYGEN_COUNTS.get(key, -1)
            if refOCount >= 0:
                diff = fragmentOxygens - refOCount
                if diff > 0 or diff < -1:
                    continue

            matched = True
            break  # 只需确认匹配存在, 不需多个 match

        if matched:
            weight = BAYESIAN_PRIOR_WEIGHT.get(name, 50)
            validMatches.append((name, anomer, weight))
            # 不 break! 继续遍历, 收集所有合法候选

    # ================================================================
    # 贝叶斯裁决: 最高权重者胜出
    # ================================================================
    if not validMatches:
        baseName = "Hex" if fragmentRingSize == 6 else "Pen"
        return baseName, "?"

    # 按权重排序 (降序), 权重相同时按字典序
    validMatches.sort(key=lambda x: (-x[2], x[0]))
    bestName, bestAnomer, _ = validMatches[0]
    return bestName, bestAnomer


def _chiralScore(mol, ring_atoms, refSmiles: str) -> float:
    """计算碎片与参考糖在环内碳 (C2-C5) 上的手性重合度。
    Compute chiral overlap score between fragment and ref sugar on ring carbons.
    返回 0.0 ~ 1.0:
      1.0 = 所有已定义的手性中心完全一致
      0.0 = 完全不匹配
    对于碎片上 CHI_UNSPECIFIED 的碳, 视为"不确定", 不扣分。
    For CHI_UNSPECIFIED atoms in fragment, count as "unknown" (neither match nor mismatch).
    """
    refMol = Chem.MolFromSmiles(refSmiles)
    if refMol is None:
        return 0.0

    # 获取参考分子的环内碳手性
    refRi = refMol.GetRingInfo()
    refRingAtoms = None
    for ring in refRi.AtomRings():
        if len(ring) in (5, 6):
            oInRing = sum(1 for i in ring if refMol.GetAtomWithIdx(i).GetAtomicNum() == 8)
            if oInRing >= 1:
                refRingAtoms = list(ring)
                break
    if refRingAtoms is None:
        return 0.0

    # 排除异头碳 (连环O + 环外O 的碳), 只对比 C2-C5
    refRingSet = set(refRingAtoms)
    refChirals = {}  # idx_in_ring -> ChiralType
    for idx in refRingAtoms:
        atom = refMol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:
            continue
        # 判断是否为异头碳
        hasRingO = any(n.GetAtomicNum() == 8 and n.GetIdx() in refRingSet for n in atom.GetNeighbors())
        hasExoO = any(n.GetAtomicNum() == 8 and n.GetIdx() not in refRingSet for n in atom.GetNeighbors())
        if hasRingO and hasExoO:
            continue  # 跳过异头碳
        tag = atom.GetChiralTag()
        if tag != Chem.ChiralType.CHI_UNSPECIFIED:
            refChirals[idx] = tag

    if not refChirals:
        return 0.0

    # 2D 骨架匹配 (useChirality=False) 获得原子映射
    refAsSmarts = Chem.MolFromSmarts(Chem.MolToSmarts(
        Chem.MolFromSmiles(refSmiles)))
    if refAsSmarts is None:
        return 0.0
    matches = mol.GetSubstructMatches(refAsSmarts, useChirality=False)
    if not matches:
        return 0.0

    # 选择覆盖最多 ring_atoms 的匹配
    bestMatch = None
    bestOverlap = -1
    ringSet = set(ring_atoms)
    for match in matches:
        overlap = sum(1 for idx in match if idx in ringSet)
        if overlap > bestOverlap:
            bestOverlap = overlap
            bestMatch = match

    if bestMatch is None:
        return 0.0

    # 对比手性
    totalChirals = len(refChirals)
    matchCount = 0
    unknownCount = 0

    for refIdx, refTag in refChirals.items():
        fragIdx = bestMatch[refIdx] if refIdx < len(bestMatch) else -1
        if fragIdx < 0 or fragIdx >= mol.GetNumAtoms():
            unknownCount += 1
            continue
        fragTag = mol.GetAtomWithIdx(fragIdx).GetChiralTag()
        if fragTag == Chem.ChiralType.CHI_UNSPECIFIED:
            unknownCount += 1  # 碎片手性缺失 → 不扣分
        elif fragTag == refTag:
            matchCount += 1
        # else: mismatch, 不计分

    # 得分: (匹配数 + 未知数*0.5) / 总手性中心数
    # 未知手性给 0.5 分 (有手性到来的可能)
    if totalChirals == 0:
        return 0.0
    score = (matchCount + unknownCount * 0.5) / totalChirals
    return score


def _tier3FuzzyChiralMatch(mol, virtual_mol, ring_atoms, fragmentRingSize,
                           fragC, fragN, fragmentOxygens):
    """Tier 3: 模糊手性打分引擎 — 仅限常见糖
    Tier 3: Fuzzy Chiral Scoring Engine — COMMON sugars only

    1. useChirality=False 确认骨架匹配
    2. 计算手性重合度 (chiralScore)
    3. score >= 0.75 → 接受匹配
    4. score < 0.75 → 返回 Hex (诚实退避)
    """
    base_name = "Hex" if fragmentRingSize == 6 else "Pen"
    matched_anomer = "?"

    bestScore = 0.0
    bestName = None
    bestAnomer = None

    for key in SPECIFICITY_ORDER:
        name, anomer = key

        # 仅限常见糖 (Rare quarantine)
        if name not in COMMON_HEXOSES and fragmentRingSize == 6:
            # 对非己糖 (戊糖等), 也允许常见的 D-Xyl, L-Ara 等
            if fragmentRingSize != 5:
                continue
        if name in RARE_HEXOSES:
            continue

        refMol = REFERENCE_MOLS.get(key)
        if refMol is None:
            continue

        # 环大小门控
        refRingSize = REFERENCE_RING_SIZE.get(key, 6)
        if refRingSize != fragmentRingSize:
            continue

        # 碳数绝对相等
        refC, refN = REFERENCE_CN_COUNTS.get(key, (-1, -1))
        if refC >= 0:
            if fragC != refC:
                continue
            if fragN > 0 and refN == 0:
                continue
            if fragN == 0 and refN > 0:
                continue

        # Step 1: useChirality=False 骨架匹配
        matches = virtual_mol.GetSubstructMatches(refMol, useChirality=False)
        hasRingMatch = False
        for match in matches:
            if all(idx in match for idx in ring_atoms):
                hasRingMatch = True
                break
        if not hasRingMatch:
            continue

        # O-Gate
        refOCount = REFERENCE_OXYGEN_COUNTS.get(key, -1)
        if refOCount >= 0:
            diff = fragmentOxygens - refOCount
            if diff > 0 or diff < -1:
                continue

        # Step 2: 计算手性重合度
        refSmiles = RAW_MONOSACCHARIDE_SMILES.get(key, "")
        if not refSmiles:
            continue
        score = _chiralScore(mol, ring_atoms, refSmiles)

        if score > bestScore:
            bestScore = score
            bestName = name
            bestAnomer = anomer

    # Step 3: 判断得分
    if bestScore >= 0.75 and bestName is not None:
        base_name = bestName
        matched_anomer = bestAnomer
    # else: 保持 Hex/Pen — 诚实退避

    return base_name, matched_anomer


def _applyFuzzyAndMods(mol, ring_atoms, base_name, matched_anomer):
    """碳计数退避 + 模糊分类 + 修饰检测 (从 identify_monosaccharide_v2 提取)。
    Carbon-count fallback + fuzzy classifier + modification detection.
    """
    # === 氨基糖 N 验证: 防止无 N 目标匹配含 N 参考 (Post-match N guard) ===
    # 设计意图: L-Vancosamine/D-Myc/D-Des 等含 N 糖模板在 N→O 虚拟突变后
    # 可能误匹配不含 N 的糖环 (如 L-Ole). 此处做最终验证: 如果 base_name
    # 是已知含 N 的氨基糖, 但目标环没有环外 N 原子, 则降级为 Hex/Pen.
    # Design: Amino sugar templates (with N) can accidentally match non-N rings
    # after N→O virtual mutation. Final guard: if name is amino sugar but ring
    # has no exocyclic N, degrade to generic Hex/Pen.
    _AMINO_SUGAR_NAMES = {
        "L-Vancosamine", "D-Myc", "D-Des", "L-Aco", "D-Bac",
        "D-GlcN", "L-GlcN", "D-GalN", "L-GalN", "D-ManN", "L-ManN",
        "D-GlcNAc", "L-GlcNAc", "D-GalNAc", "L-GalNAc", "D-ManNAc",
    }
    pureBase = base_name.split("(")[0] if "(" in base_name else base_name
    if pureBase in _AMINO_SUGAR_NAMES:
        ringSet = set(ring_atoms)
        hasExoN = False
        for rIdx in ring_atoms:
            rAtom = mol.GetAtomWithIdx(rIdx)
            if rAtom.GetAtomicNum() != 6:
                continue
            for nbr in rAtom.GetNeighbors():
                if nbr.GetIdx() not in ringSet and nbr.GetAtomicNum() == 7:
                    hasExoN = True
                    break
            if hasExoN:
                break
        if not hasExoN:
            # 目标环无 N → 氨基糖匹配无效, 降级 (No N on ring → invalid, degrade)
            base_name = "Hex" if len(ring_atoms) == 6 else "Pen"
            matched_anomer = "?"

    if base_name in ("Hex", "Pen") and matched_anomer == "?":
        totalC = _countFragmentCarbons(mol, ring_atoms)
        # 不再使用 Non (9C): 碳计数超标通常是苷元渗透, 不是真正的壬糖
        # No longer classify as Non (9C): excessive count is usually aglycon bleed
        if totalC == 8:
            base_name = "Oct"  # 辛糖 (8C, 如 KDO 类)
        elif totalC == 7:
            base_name = "Hept"  # 庚糖 (7C)
        # 5C, 6C, 9C+ 保持 Pen/Hex (Keep original Pen/Hex)

    # ---- 模糊分类器 3.0 (Fuzzy Classifier 3.0) ----
    # 替代旧的 Unknown_Sugar_MW_86 兜底逻辑
    # 使用碳/氮/氧原子计数对无法精确命名的糖进行化学分类
    # Replaces MW_86 fallback with chemistry-aware classification
    if base_name in ("Hex", "Pen") and matched_anomer == "?":
        # 收集完整糖单元原子 (环 + 直连外环原子)
        fullUnitAtoms = set(ring_atoms)
        for rIdx in ring_atoms:
            for nbr in mol.GetAtomWithIdx(rIdx).GetNeighbors():
                fullUnitAtoms.add(nbr.GetIdx())
        fragC, fragN = _countFragmentCN(mol, fullUnitAtoms)
        # 核苷糖修正: 如果 N 属于另一个环 (核碱基), 则不算氨基糖 N
        # Nucleoside fix: If N belongs to another ring (nucleobase), don't count as amino N
        # 设计意图: 核苷 (如胸腺嘧啶/胞嘧啶) 通过 N-糖苷键连接糖环,
        # 这个 N 是糖苷键 (类似 O-糖苷键的 O), 不是氨基修饰
        # Design: Nucleosides connect via N-glycosidic bond (analogous to O-glycosidic).
        # This N is a linkage, NOT an amino modification.
        ri = mol.GetRingInfo()
        allRingsSet = [set(r) for r in ri.AtomRings()]
        ringSet = set(ring_atoms)
        adjustedFragN = fragN
        for idx in fullUnitAtoms:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 7 and idx not in ringSet:
                # 检查这个 N 是否在任何其他环中 (Check if N is in another ring)
                inOtherRing = any(idx in rs for rs in allRingsSet if not rs.issubset(ringSet))
                if inOtherRing:
                    adjustedFragN -= 1
        fragN = max(0, adjustedFragN)
        # 计算含氧量 (O/C 比率判断脱氧程度)
        fragO = sum(1 for idx in fullUnitAtoms
                    if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        ringSize = len(ring_atoms)

        if ringSize == 6:  # 吡喃糖环
            if fragN > 0:
                base_name = "HexN"   # 氨基己糖 (有 N 原子)
            elif fragO <= 2:  # 环氧 + 最多 1 个 exo-O → 高度脱氧
                base_name = "dHex"   # 脱氧己糖
            else:
                base_name = "Hex"    # 普通己糖
        elif ringSize == 5:  # 呋喃糖环
            if fragN > 0:
                base_name = "PenN"   # 氨基戊糖 (罕见但化学完备)
            else:
                base_name = "Pen"    # 普通戊糖

    mods, is_acid, is_complex = check_modifications(mol, ring_atoms)

    if "GlcNAc" in base_name or "GalNAc" in base_name:
        if "NAc" in mods: mods.remove("NAc")
    if "GlcA" in base_name or "GalA" in base_name or "IdoA" in base_name:
        if "A" in mods: mods.remove("A")

    if "Glc" in base_name and "A" in mods:
        base_name = base_name.replace("Glc", "GlcA")
        mods.remove("A")
    elif "Gal" in base_name and "A" in mods:
        base_name = base_name.replace("Gal", "GalA")
        mods.remove("A")

    # === 智能 NAc 合并: 验证 N 直连糖环碳 (Smart NAc merge) ===
    # 设计意图: check_modifications 检测到 NAc 时, 验证 N 原子确实直连在糖环碳上
    # (真正的 GlcNAc: ring-C2-NH-CO-CH3), 而不是糖苷键上的 N (假阳性)
    # Design: verify N atom is directly bonded to ring carbon (true GlcNAc)
    # vs glycosidic bond N (false positive from cross-ring linkage)
    if "NAc" in mods and ("GlcNAc" not in base_name and "GalNAc" not in base_name
                          and "ManNAc" not in base_name):
        ringSet = set(ring_atoms)
        hasRealNAc = False
        for rIdx in ring_atoms:
            rAtom = mol.GetAtomWithIdx(rIdx)
            if rAtom.GetAtomicNum() != 6:  # 只看环碳 (only ring carbons)
                continue
            for nbr in rAtom.GetNeighbors():
                if nbr.GetIdx() in ringSet:
                    continue
                # 环外 N 原子: 检查是否连了乙酰基 -C(=O)CH3
                # Exocyclic N: check if it carries acetyl group -C(=O)CH3
                if nbr.GetAtomicNum() == 7:
                    # 防护: N 的重原子度应为 2 (连 ring-C + acetyl-C)
                    # 如果 N 度数 > 2 则可能是糖苷链接 N, 跳过
                    # Guard: N heavy-atom degree should be 2 (ring-C + acetyl-C)
                    # If degree > 2, likely glycosidic linker N → skip
                    nHeavyDeg = sum(1 for x in nbr.GetNeighbors() if x.GetAtomicNum() > 1)
                    if nHeavyDeg > 2:
                        continue
                    for nNbr in nbr.GetNeighbors():
                        if nNbr.GetIdx() == rIdx:
                            continue
                        if nNbr.GetAtomicNum() == 6:
                            # 检查 C=O 双键 (check for C=O)
                            hasDoubleBondO = False
                            for bond in nNbr.GetBonds():
                                otherAtom = mol.GetAtomWithIdx(bond.GetOtherAtomIdx(nNbr.GetIdx()))
                                if otherAtom.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                                    hasDoubleBondO = True
                            if hasDoubleBondO:
                                hasRealNAc = True
                                break
                if hasRealNAc:
                    break
            if hasRealNAc:
                break
        if hasRealNAc:
            _NAC_MERGE_MAP = {
                "D-Glc": "D-GlcNAc", "L-Glc": "L-GlcNAc",
                "D-Gal": "D-GalNAc", "L-Gal": "L-GalNAc",
                "D-Man": "D-ManNAc", "L-Man": "L-ManNAc",
            }
            pureBase = base_name.split("(")[0] if "(" in base_name else base_name
            if pureBase in _NAC_MERGE_MAP:
                base_name = _NAC_MERGE_MAP[pureBase]
                mods.remove("NAc")

    # === 氨基糖 base-name 提升: D-Glc(N) → D-GlcN ===
    # 设计意图: 如果检测到 N 修饰且 N 直连环碳, 提升为标准氨基糖名
    # Design: If N modification detected and N is directly on ring carbon,
    # promote to standard amino sugar name (D-GlcN, D-GalN)
    if "N" in mods and ("GlcN" not in base_name and "GalN" not in base_name
                        and "ManN" not in base_name):
        _AMINO_MERGE_MAP = {
            "D-Glc": "D-GlcN", "L-Glc": "L-GlcN",
            "D-Gal": "D-GalN", "L-Gal": "L-GalN",
            "D-Man": "D-ManN", "L-Man": "L-ManN",
        }
        pureBase = base_name.split("(")[0] if "(" in base_name else base_name
        if pureBase in _AMINO_MERGE_MAP:
            base_name = _AMINO_MERGE_MAP[pureBase]
            mods.remove("N")

    # === 直接 N 检测提升: N→O 突变隐藏的氨基糖 (Direct N detection promotion) ===
    # 设计意图: 当引擎通过 N→O 虚拟突变匹配了基础己糖 (D-Glc/D-Rha 等),
    # 但原始分子的糖环碳上确实有直连的 N 原子 (NH2/NHR), 则提升为氨基糖
    # Design: When engine matched base hexose via N→O mutation, but original
    # ring carbon has direct N atom (NH2/NHR), promote to amino sugar variant
    _DIRECT_N_PROMOTE_MAP = {
        "D-Glc": "D-GlcN", "L-Glc": "L-GlcN",
        "D-Gal": "D-GalN", "L-Gal": "L-GalN",
        "D-Man": "D-ManN", "L-Man": "L-ManN",
        "D-Rha": "D-GlcN",   # 3-amino-3-deoxy → kanosamine ≈ D-GlcN
        "D-Qui": "D-GlcN",   # 4-amino-4,6-dideoxy → amino-D-Glc variant
    }
    pureBase = base_name.split("(")[0] if "(" in base_name else base_name
    if pureBase in _DIRECT_N_PROMOTE_MAP and "GlcN" not in base_name and "GalN" not in base_name:
        ringSet = set(ring_atoms)
        hasDirectN = False
        for rIdx in ring_atoms:
            rAtom = mol.GetAtomWithIdx(rIdx)
            if rAtom.GetAtomicNum() != 6:
                continue
            for nbr in rAtom.GetNeighbors():
                if nbr.GetIdx() not in ringSet and nbr.GetAtomicNum() == 7:
                    hasDirectN = True
                    break
            if hasDirectN:
                break
        if hasDirectN:
            base_name = _DIRECT_N_PROMOTE_MAP[pureBase]
            # 移除冗余 N 修饰标签 (remove redundant N mod tag)
            if "N" in mods:
                mods.remove("N")

    # Bug 2 修复: 脱氧糖固有属性过滤 — L-Rha 不应标注 deoxy
    # Bug 2 fix: Remove redundant 'deoxy' tag for inherently-deoxy sugars.
    # 设计意图: 鼠李糖 (L-Rha = 6-deoxy-L-Man) 本身就是脱氧糖, deoxy 不是修饰
    # L-Rha IS a 6-deoxy sugar by definition; 'deoxy' is not an external modification.
    INHERENTLY_DEOXY = {
        "L-Rha", "D-Rha", "L-Fuc", "D-Fuc",
        "D-Qui", "L-Qui",
        "D-Dtx", "L-Dtx", "D-Dig", "L-Ole", "D-Cym", "L-The",
        "D-Boi", "D-Sar", "D-Cma", "D-Din", "L-Ola", "D-Eva", "L-Eva",
        "D-dGlc", "L-dGlc", "D-dGal", "L-dGal", "D-dMan",
        "dHex",  # 泛指脱氧己糖也不需要 deoxy 标签
        "L-Streptose",  # 分支脱氧醛戊糖 (branched deoxy aldopentose)
        "D-Myc",  # 3-amino-3,6-dideoxy-D-mannose
        "D-Des",  # desosamine — 脱氧氨基糖
    }
    pureBaseForDeoxy = base_name.split("(")[0] if "(" in base_name else base_name
    if pureBaseForDeoxy in INHERENTLY_DEOXY and "deoxy" in mods:
        mods.remove("deoxy")

    # 固有氨基糖过滤 (Inherently Amino Sugar Filter)
    # Kanosamine 和 6aGlc 本身就自带氨基，修饰扫描器产生的 'N' 需要被剔除
    INHERENTLY_AMINO = {
        "D-Kanosamine", "D-6aGlc", "D-Myc", "D-Des", "L-Vancosamine", "L-Aco", "D-Bac",
    }
    if pureBaseForDeoxy in INHERENTLY_AMINO and "N" in mods:
        mods.remove("N")

    final_name = f"{base_name}({','.join(mods)})" if mods else base_name
    return final_name, matched_anomer


# =====================================================================
# v10 工业级裸糖分离与匹配管线 (Industrial-Grade Bare Sugar Pipeline)
# =====================================================================

# 需要导入修饰剥离字典
try:
    from glycan_reference_library import STRIPPING_REACTIONS
except ImportError:
    from lib.glycan_reference_library import STRIPPING_REACTIONS

# O-酯/醚修饰: 剥离后还原为 -OH
# N-修饰 (N-Acetylated 等): 不剥离, 保留给 GlcNAc/GalNAc 精确匹配
# 关键安全规则 (Critical Safety Rules):
#   - "Lactylated" 已移除: 其 SMARTS 匹配糖苷键 (C-O-C), 导致多糖被错误剥离!
#   - "Methylated" 已移除: 在多糖中 O-Me 匹配糖苷键上的 O-C, 需要用直接扫描器处理
_O_STRIPPABLE_MODS = {
    "Sulfated", "Phosphated",
    "Acetylated", "Formylated", "Malonylated", "Succinylated",
    "Tigloylated",
    "Galloylated", "Benzoylated", "p-Coumaroylated",
    "Caffeoylated", "Feruloylated", "Sinapoylated",
}


def isolateBareSugarRing(mol, ringAtoms):
    """核心重构 1: 真·糖环切除与断口羟基化
    Core Refactor 1: True sugar ring cleavage with -OH capping.

    算法 (Algorithm):
    1. 收集 ring + 第一层 exo 原子 (包括 C6-CH₂OH)
    2. BFS 从 exo 原子向外遍历整个"尾巴"
    3. 用 FragmentOnBonds 切断尾巴键
    4. 提取含糖环的碎片, dummy→O (补 -OH)
    5. SanitizeMol + AssignStereochemistry

    Args:
        mol: 完整分子 RDKit Mol
        ringAtoms: 糖环原子索引列表 (list of int)

    Returns:
        (bareMol, ringAtomsNew): 裸糖分子 + 新原子索引映射, 或 (None, None)
    """
    ringSet = set(ringAtoms)

    # Step 1: 收集 ring + 直连 exo (包括 C6 和 C6-OH)
    exoAtoms = set()
    c6Candidates = set()
    for rIdx in ringSet:
        atom = mol.GetAtomWithIdx(rIdx)
        for nbr in atom.GetNeighbors():
            nIdx = nbr.GetIdx()
            if nIdx not in ringSet:
                exoAtoms.add(nIdx)
                # C6 检测: 环碳 → 非环碳 → 碳上连有 O/N
                if nbr.GetAtomicNum() == 6:
                    for nnbr in nbr.GetNeighbors():
                        if nnbr.GetIdx() not in ringSet and nnbr.GetIdx() != rIdx:
                            if nnbr.GetAtomicNum() in (8, 7, 16):
                                c6Candidates.add(nbr.GetIdx())
                                exoAtoms.add(nnbr.GetIdx())

    # 糖核心原子 = ring + exo 第一层
    sugarCoreAtoms = ringSet | exoAtoms

    # Step 2: 找到所有需要切断的键 (从 exo 原子向外的键)
    bondsToCut = []
    for exoIdx in exoAtoms:
        exoAtom = mol.GetAtomWithIdx(exoIdx)
        for nbr in exoAtom.GetNeighbors():
            nIdx = nbr.GetIdx()
            # 不切 ring 内的键, 不切 C6 与 환 的键
            if nIdx in ringSet:
                continue
            if nIdx in c6Candidates or exoIdx in c6Candidates:
                # C6 与其 OH 之间的键不切 (保留 C6-OH)
                if nIdx in exoAtoms:
                    continue
            if nIdx not in sugarCoreAtoms:
                bond = mol.GetBondBetweenAtoms(exoIdx, nIdx)
                if bond:
                    bondsToCut.append(bond.GetIdx())

    if not bondsToCut:
        # 没有外接键 → 整个分子就是糖 (或极小分子)
        return Chem.Mol(mol), list(range(len(ringAtoms)))

    # Step 3: FragmentOnBonds 切割
    try:
        fragMol = Chem.FragmentOnBonds(mol, bondsToCut)
        # Bug 1/24 Fix: sanitizeFrags=False prevents kekulization crash on nucleobases
        fragMols = Chem.GetMolFrags(fragMol, asMols=True, sanitizeFrags=False)
        fragIndices = Chem.GetMolFrags(fragMol, asMols=False)
    except Exception:
        return Chem.Mol(mol), list(range(len(ringAtoms)))

    # Step 4: 找到含糖环的碎片
    for frag, indices in zip(fragMols, fragIndices):
        # 检查是否包含所有 ring 原子
        if not all(r in indices for r in ringAtoms):
            continue
            
        try:
            # 仅对糖碎片执行 sanitize
            Chem.SanitizeMol(frag)
        except Exception:
            pass

        # 验证含有有效糖环
        ri = frag.GetRingInfo()
        hasValidRing = any(
            len(r) in (5, 6) and
            sum(1 for a in r if frag.GetAtomWithIdx(a).GetAtomicNum() == 8) >= 1
            for r in ri.AtomRings()
        )
        if not hasValidRing:
            continue

        # Step 5: dummy→OH (补 -OH, 避免 O-O 过氧化物)
        # When dummy's neighbour is C: replace dummy with O → creates C-OH ✓
        # When dummy's neighbour is O: replace dummy with H → keeps O-H (the glycosidic O IS the OH) ✓
        rwmol = Chem.RWMol(frag)
        for atom in rwmol.GetAtoms():
            if atom.GetAtomicNum() == 0 and atom.GetNeighbors():
                nbr = atom.GetNeighbors()[0]
                nbrNum = nbr.GetAtomicNum()
                if nbrNum == 6:
                    # C-[dummy] → C-OH: dummy becomes O
                    atom.SetAtomicNum(8)
                else:
                    # O-[dummy] or N-[dummy] → O-H or N-H: dummy becomes H
                    atom.SetAtomicNum(1)
                atom.SetFormalCharge(0)
                atom.SetIsotope(0)

        # Step 5b: 清除异头碳手性 (切断糖苷键后 C1 的 CIP 不再有意义)
        # Clear anomeric carbon chirality: after cleavage, CIP at C1 may flip
        # due to substituent change → must unspecify to match SMARTS [C:1]
        fragRi = rwmol.GetRingInfo()
        for ring in fragRi.AtomRings():
            if len(ring) not in (5, 6):
                continue
            ringS = set(ring)
            for aIdx in ring:
                a = rwmol.GetAtomWithIdx(aIdx)
                if a.GetAtomicNum() != 6:
                    continue
                hasRingO = any(
                    n.GetAtomicNum() == 8 and n.GetIdx() in ringS
                    for n in a.GetNeighbors()
                )
                hasExoO = any(
                    n.GetAtomicNum() == 8 and n.GetIdx() not in ringS
                    for n in a.GetNeighbors()
                )
                if hasRingO and hasExoO:
                    a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)

        # Step 6: SanitizeMol + 手性重建
        try:
            rwmol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(rwmol)
            cleanMol = Chem.RemoveHs(rwmol.GetMol())
            Chem.AssignStereochemistry(cleanMol, force=True, cleanIt=True)

            # Step 6b: 清除异头碳手性 (必须在 AssignStereochemistry 之后)
            # 切断糖苷键后 C1 的 CIP 可能翻转 → 清除以匹配 [C:1] 通配符
            editMol = Chem.RWMol(cleanMol)
            fragRi2 = editMol.GetRingInfo()
            for ring2 in fragRi2.AtomRings():
                if len(ring2) not in (5, 6):
                    continue
                ringS2 = set(ring2)
                for aIdx2 in ring2:
                    a2 = editMol.GetAtomWithIdx(aIdx2)
                    if a2.GetAtomicNum() != 6:
                        continue
                    hasRingO2 = any(
                        n.GetAtomicNum() == 8 and n.GetIdx() in ringS2
                        for n in a2.GetNeighbors()
                    )
                    hasExoO2 = any(
                        n.GetAtomicNum() == 8 and n.GetIdx() not in ringS2
                        for n in a2.GetNeighbors()
                    )
                    if hasRingO2 and hasExoO2:
                        a2.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
            cleanMol = editMol.GetMol()
        except Exception:
            cleanMol = rwmol.GetMol()

        # 构建新旧原子映射 (用于后续 ring_atoms 追踪)
        # indices[newIdx] = oldIdx_in_fragMol
        oldToNew = {}
        for newIdx, oldIdx in enumerate(indices):
            if oldIdx < mol.GetNumAtoms():
                oldToNew[oldIdx] = newIdx

        newRingAtoms = [oldToNew.get(r, -1) for r in ringAtoms]
        if -1 in newRingAtoms:
            # 映射失败 → 用 substructure match 恢复
            matches = mol.GetSubstructMatches(cleanMol, useChirality=False)
            bestMatch = next((m for m in matches
                              if sum(1 for t in m if t in ringSet) >= len(ringSet)),
                             None)
            if bestMatch:
                newRingAtoms = list(range(len(ringAtoms)))

        return cleanMol, newRingAtoms

    return None, None


def stripModificationsToBareSugar(fragMol, ringAtoms):
    """核心重构 2: 真·虚拟水解酶
    Core Refactor 2: Virtual Hydrolase — strip modifications, restore -OH.

    扫描 STRIPPING_REACTIONS 中的 O-酯/醚修饰 (不剥离 N-Acetyl 等 N-修饰)。
    对匹配到的修饰: 找到连接点 → 删除修饰基团原子 → 保留 O (还原为 -OH)。

    Args:
        fragMol: 隔离后的糖碎片 (isolateBareSugarRing 的输出)
        ringAtoms: 糖环原子索引列表

    Returns:
        (bareMol, modList): 裸糖分子 + 检测到的修饰列表
    """
    modList = []
    currentMol = Chem.RWMol(fragMol)

    for modName, modPat in STRIPPING_REACTIONS.items():
        if modName not in _O_STRIPPABLE_MODS:
            continue
        if modPat is None:
            continue

        # 反复扫描直到没有本地匹配 (处理多位点修饰)
        maxIter = 5
        for _ in range(maxIter):
            matches = currentMol.GetSubstructMatches(modPat)
            if not matches:
                break

            valid_match = None
            ringSet = set(ringAtoms) if ringAtoms else set()

            for match in matches:
                anchorIdx = match[0]
                anchorAtom = currentMol.GetAtomWithIdx(anchorIdx)
                
                belongs_to_ring = False
                if anchorIdx in ringSet:
                    belongs_to_ring = True
                else:
                    for nbr in anchorAtom.GetNeighbors():
                        if nbr.GetIdx() in ringSet:
                            belongs_to_ring = True
                            break
                        if nbr.GetAtomicNum() == 6: # e.g. C6
                            for c6nbr in nbr.GetNeighbors():
                                if c6nbr.GetIdx() in ringSet:
                                    belongs_to_ring = True
                                    break
                        if belongs_to_ring:
                            break
                            
                if belongs_to_ring:
                    valid_match = match
                    break
                    
            if not valid_match:
                break # No local matches found for this ring

            match = valid_match
            anchorIdx = match[0]

            # 找到修饰基团的"远端"原子 (不是环原子的那一侧)
            modAtoms = set()
            for mIdx in match:
                if mIdx != anchorIdx and mIdx not in ringSet:
                    modAtoms.add(mIdx)

            if not modAtoms:
                break

            # BFS 从 modAtoms 扩展, 不越过 anchor 和 ring
            barrier = ringSet | {anchorIdx}
            queue = list(modAtoms)
            visited = set(modAtoms)
            while queue:
                curr = queue.pop(0)
                for nbr in currentMol.GetAtomWithIdx(curr).GetNeighbors():
                    nIdx = nbr.GetIdx()
                    if nIdx not in visited and nIdx not in barrier:
                        visited.add(nIdx)
                        queue.append(nIdx)

            allModAtoms = visited

            # 删除修饰基团原子 (从大到小删, 避免索引偏移)
            atomsToRemove = sorted(allModAtoms, reverse=True)
            for aIdx in atomsToRemove:
                # 先删除所有连接键
                atom = currentMol.GetAtomWithIdx(aIdx)
                bondsToRemove = []
                for bond in atom.GetBonds():
                    otherIdx = bond.GetOtherAtomIdx(aIdx)
                    if otherIdx not in allModAtoms:
                        bondsToRemove.append((aIdx, otherIdx))
                for a1, a2 in bondsToRemove:
                    currentMol.RemoveBond(a1, a2)

            # 批量删除原子 (从大到小)
            for aIdx in atomsToRemove:
                currentMol.RemoveAtom(aIdx)

            # 更新 ringAtoms (原子索引可能偏移)
            # 每次删除原子后, 比它大的索引全部 -1
            newRingAtoms = list(ringAtoms) if ringAtoms else []
            for removedIdx in sorted(allModAtoms):
                newRingAtoms = [r - 1 if r > removedIdx else r for r in newRingAtoms]
            ringAtoms = newRingAtoms

            modList.append(modName)

            # 重新 Sanitize
            try:
                currentMol.UpdatePropertyCache(strict=False)
                Chem.SanitizeMol(currentMol)
            except Exception:
                pass

    # 最终清理
    try:
        finalMol = currentMol.GetMol()
        Chem.AssignStereochemistry(finalMol, force=True, cleanIt=True)
        return finalMol, modList
    except Exception:
        return currentMol.GetMol(), modList


def matchBareSugar(bareMol, ringAtoms):
    """核心重构 3: 极简裸糖精确匹配
    Core Refactor 3: Bare sugar precise matching — no O-Gate, no C-Gate.

    因为裸糖已经完美隔离, 碳数/氧数天然正确, 只需简单子结构匹配。
    保留贝叶斯裁决防止残缺手性撞车。

    Tier 1: 严格手性 (所有糖)
    Tier 2: 异头碳特赦 (仅 COMMON)
    Tier 3: 模糊手性打分 (仅 COMMON, score >= 0.75)

    Args:
        bareMol: 裸糖 Mol
        ringAtoms: 裸糖中的环原子索引

    Returns:
        (sugarName, anomer)
    """
    fragmentRingSize = len(ringAtoms) if ringAtoms else 6

    # 异头碳手性清除 (在匹配前统一执行)
    # Anomeric carbon chirality clearing (unified for all pathways)
    # 切断糖苷键后 C1 的 CIP 可能翻转 → 必须清除
    ringSet = set(ringAtoms)
    try:
        rwBare = Chem.RWMol(bareMol)
        for aIdx in ringAtoms:
            a = rwBare.GetAtomWithIdx(aIdx)
            if a.GetAtomicNum() != 6:
                continue
            hasRingO = any(
                n.GetAtomicNum() == 8 and n.GetIdx() in ringSet
                for n in a.GetNeighbors()
            )
            hasExoO = any(
                n.GetAtomicNum() == 8 and n.GetIdx() not in ringSet
                for n in a.GetNeighbors()
            )
            if hasRingO and hasExoO:
                a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        bareMol = rwBare.GetMol()
    except Exception:
        pass

    # 计算裸糖的实际 Gate 值 (真实 O/C/N 计数)
    # Real O/C/N counts from bare sugar — these should perfectly match references
    fragmentOxygens = _countFragmentRingOxygens(bareMol, ringSet)
    fullUnitAtoms = set(ringAtoms)
    for rIdx in ringAtoms:
        for nbr in bareMol.GetAtomWithIdx(rIdx).GetNeighbors():
            fullUnitAtoms.add(nbr.GetIdx())
    fragC, fragN = _countFragmentCN(bareMol, fullUnitAtoms)

    # 创建虚拟标准糖 (N/S→O 突变, 为了匹配 GlcN→Glc 模板)
    virtualMol = create_virtual_standard_sugar(bareMol, ringAtoms)

    # --- Tier 1: 严格手性 (全量竞争 + 贝叶斯裁决) ---
    name, anomer = _tierMatch(
        virtualMol, ringAtoms, fragmentRingSize, fragC, fragN, fragmentOxygens,
        REFERENCE_MOLS, useChirality=True, allowRare=True
    )
    if anomer != "?":
        return name, anomer

    # --- Tier 2: 异头碳特赦 (仅 COMMON) ---
    name, anomer = _tierMatch(
        virtualMol, ringAtoms, fragmentRingSize, fragC, fragN, fragmentOxygens,
        AMNESTY_REFERENCE_MOLS, useChirality=True, allowRare=False
    )
    if anomer != "?":
        return name, anomer

    # --- Tier 3: 模糊手性打分 ---
    name, anomer = _tier3FuzzyChiralMatch(
        bareMol, virtualMol, ringAtoms, fragmentRingSize, fragC, fragN, fragmentOxygens
    )
    if anomer != "?":
        return name, anomer

    # --- CIP 救援 ---
    name, anomer = rescue_hexose_by_key_nodes(bareMol, ringAtoms)
    if anomer != "?":
        return name, anomer

    # 退化
    return ("Hex" if fragmentRingSize == 6 else "Pen"), "?"


# 修饰标签简化映射 (Human-readable modification labels)
_MOD_LABEL_MAP = {
    "Sulfated": "S", "Phosphated": "P",
    "Acetylated": "Ac", "Formylated": "Fo",
    "Malonylated": "Mal", "Succinylated": "Suc",
    "Lactylated": "Lac", "Tigloylated": "Tig",
    "Galloylated": "Gal_ester", "Benzoylated": "Bz",
    "p-Coumaroylated": "pCou", "Caffeoylated": "Caf",
    "Feruloylated": "Fer", "Sinapoylated": "Sin",
    "Methylated": "Me",
}


def identify_monosaccharide_v10(mol, ringAtoms):
    import sys
    # print("DEBUG v10 START: length", len(ringAtoms), file=sys.stderr)
    GENERIC_LABELS = {"Hex", "Pen", "dHex", "HexN", "PenN", "HexA"}

    # ================================================================
    # Step 1: 无损基线 — 直接在原始分子上用 v2 匹配
    # Lossless baseline: run v2 directly on original molecule
    # ================================================================
    baselineResult = identify_monosaccharide_v2(mol, ringAtoms)
    # print(f"DEBUG v10 V2 Result: {baselineResult}", file=sys.stderr)
    baselineName = baselineResult[0] if isinstance(baselineResult, tuple) else baselineResult
    # 去除修饰后缀以提取纯糖名, 例如 "D-Glc(deoxy)" → "D-Glc"
    # Strip mod suffix for comparison, e.g. "D-Glc(deoxy)" → "D-Glc"
    pureBaseName = baselineName.split("(")[0] if "(" in baselineName else baselineName
    # print(f"DEBUG v10 pureBaseName: {pureBaseName}", file=sys.stderr)

    # ================================================================
    # Step 2: 检测修饰基团 (仅用于标签装饰, 不影响匹配结果)
    # Detect modifications (for labeling only, does NOT affect matching)
    # ================================================================
    try:
        taggedMol = Chem.RWMol(mol)
        for atom in taggedMol.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx() + 1)
        taggedMol = taggedMol.GetMol()
        _, modList = stripModificationsToBareSugar(taggedMol, list(ringAtoms))
    except Exception:
        modList = []

    # ================================================================
    # Step 2.1: 显式 Sulfate 检测 (Bug 3 修复)
    # Explicit Sulfate Detection (Bug 3 Fix)
    # stripModificationsToBareSugar 使用 SMARTS [O:1]S(=O)(=O)[O-,OH] 匹配
    # 直连环碳的 O-Sulfate, 但无法匹配通过 C6 间接连接的 6-Sulfate
    # (Ring-C → C6 → O → S(=O)(=O)). 此处显式扫描这种模式。
    # SMARTS-based detection misses C6-linked sulfates. Scan explicitly.
    # ================================================================
    if "Sulfated" not in modList:
        ringSet = set(ringAtoms)
        for rIdx in ringAtoms:
            rAtom = mol.GetAtomWithIdx(rIdx)
            if rAtom.GetAtomicNum() != 6:
                continue
            for nbr in rAtom.GetNeighbors():
                nIdx = nbr.GetIdx()
                if nIdx in ringSet:
                    continue
                # 检查两层: exo-C(C6) → O → S, 或 exo-O → S
                # Check two hops: exo-C(C6) → O → S, or exo-O → S
                atomsToCheck = [nbr]
                if nbr.GetAtomicNum() == 6:  # C6 碳
                    atomsToCheck.extend([
                        n for n in nbr.GetNeighbors()
                        if n.GetIdx() != rIdx and n.GetIdx() not in ringSet
                    ])
                for checkAtom in atomsToCheck:
                    if checkAtom.GetAtomicNum() == 8:  # O
                        for sNbr in checkAtom.GetNeighbors():
                            if sNbr.GetAtomicNum() == 16:  # S
                                # 验证 S 连接 ≥2 个 =O (sulfate/sulfonate)
                                # Verify S has ≥2 double-bonded O (sulfate)
                                dblO = sum(
                                    1 for so in sNbr.GetNeighbors()
                                    if so.GetAtomicNum() == 8
                                    and mol.GetBondBetweenAtoms(
                                        sNbr.GetIdx(), so.GetIdx()
                                    ).GetBondTypeAsDouble() == 2.0
                                )
                                if dblO >= 2:
                                    modList.append("Sulfated")
                                    break
                    if "Sulfated" in modList:
                        break
                if "Sulfated" in modList:
                    break
            if "Sulfated" in modList:
                break

    # ================================================================
    # Step 2.5: 氨基糖 N 验证 (Amino Sugar N-Validation Guard)
    # 防止含 N 糖模板 (L-Vancosamine, D-Myc 等) 误匹配无 N 的糖环
    # Prevents amino sugar templates from falsely matching non-N rings
    # (e.g. L-Vancosamine template matching L-Ole after N→O mutation)
    # ================================================================
    _AMINO_SUGAR_NAMES_V10 = {
        "L-Vancosamine", "D-Myc", "L-Aco",
    }
    if pureBaseName in _AMINO_SUGAR_NAMES_V10:
        ringSet = set(ringAtoms)
        hasExoN = any(
            any(nbr.GetAtomicNum() == 7 and nbr.GetIdx() not in ringSet
                for nbr in mol.GetAtomWithIdx(rIdx).GetNeighbors())
            for rIdx in ringAtoms
            if mol.GetAtomWithIdx(rIdx).GetAtomicNum() == 6
        )
        if not hasExoN:
            # 目标环无 N → 降级为泛指, 让后续 CIP+Exo 引擎重新匹配
            # No N on ring → degrade to generic, let CIP+Exo re-match
            pureBaseName = "Hex" if len(ringAtoms) == 6 else "Pen"
            baselineName = pureBaseName
            baselineResult = (pureBaseName, "?")

    # ================================================================
    # Step 3: 如果基线已经是具体糖名 → 直接采用, 附加修饰标签
    # If baseline already specific → use directly, append mod labels
    # ================================================================
    if pureBaseName not in GENERIC_LABELS:
        result = baselineResult
        if modList:
            shortMods = sorted(set(_MOD_LABEL_MAP.get(m, m) for m in modList))
            if shortMods and not any(s in baselineName for s in shortMods):
                finalName = f"{baselineName}({','.join(shortMods)})"
                if isinstance(result, tuple):
                    return finalName, result[1]
                else:
                    return finalName
        return result

    # ================================================================
    # Step 3.2: 核苷呋喃糖专用匹配 (Bug 1 修复)
    # Nucleoside Furanose Matching (Bug 1 Fix)
    # 核苷 (Thymidine/Cytarabine 等) 中的 2-deoxy 呋喃糖因缺少 C2-OH
    # 导致 CIP 引擎手性信息不足 → Pen 退化。此处用 OH 计数 + N 邻居检测。
    # Nucleosides with 2-deoxy furanose lack C2-OH → CIP degrades to Pen.
    # Use OH count + N-neighbor detection for specific identification.
    # ================================================================
    if len(ringAtoms) == 5 and pureBaseName == "Pen":
        ringSet = set(ringAtoms)
        # 检查 C1 是否邻接 N (核碱基标志)
        # Check if C1 has an N neighbor (nucleobase signature)
        ringOIdx = None
        for rIdx in ringAtoms:
            if mol.GetAtomWithIdx(rIdx).GetAtomicNum() == 8:
                ringOIdx = rIdx
                break
        hasNucleobaseN = False
        if ringOIdx is not None:
            ringO = mol.GetAtomWithIdx(ringOIdx)
            for cNbr in ringO.GetNeighbors():
                if cNbr.GetIdx() in ringSet and cNbr.GetAtomicNum() == 6:
                    for exoNbr in cNbr.GetNeighbors():
                        if exoNbr.GetIdx() not in ringSet and exoNbr.GetAtomicNum() == 7:
                            hasNucleobaseN = True
                            break
                if hasNucleobaseN:
                    break

        if hasNucleobaseN:
            # 计算环碳上的环外 OH 数 (不含环 O, 不含 N)
            # Count exocyclic OH on ring carbons (exclude ring O, exclude N)
            exoOhCount = 0
            for rIdx in ringAtoms:
                rAtom = mol.GetAtomWithIdx(rIdx)
                if rAtom.GetAtomicNum() != 6:
                    continue
                for nbr in rAtom.GetNeighbors():
                    nIdx = nbr.GetIdx()
                    if nIdx in ringSet:
                        continue
                    if nbr.GetAtomicNum() == 8:
                        exoOhCount += 1
            # 包括 C5-CH2OH (通过 C5 → C6 → OH)
            # Include C5-CH2OH (via C5 → C6 → OH)
            for rIdx in ringAtoms:
                rAtom = mol.GetAtomWithIdx(rIdx)
                if rAtom.GetAtomicNum() != 6:
                    continue
                for nbr in rAtom.GetNeighbors():
                    nIdx = nbr.GetIdx()
                    if nIdx in ringSet:
                        continue
                    if nbr.GetAtomicNum() == 6:  # C6
                        for nnbr in nbr.GetNeighbors():
                            if nnbr.GetIdx() != rIdx and nnbr.GetAtomicNum() == 8:
                                exoOhCount += 1

            # 判定: 2 OH → dRib (2-脱氧核糖)
            # 3+ OH → D-Rib (核糖) 或 D-Ara (阿拉伯糖)
            # 根据生物学常见核苷，取决于 C2 碳的手性（D-Rib: C2为R; D-Ara: C2为S）
            nucleosideName = "Pen"
            if exoOhCount <= 2:
                nucleosideName = "dRib"
            else:
                # 尝试区分 D-Rib 和 D-Ara (C2 手性)
                c1Idx = None
                for rIdx in ringAtoms:
                    rAtom = mol.GetAtomWithIdx(rIdx)
                    if rAtom.GetAtomicNum() == 6:
                        for nbr in rAtom.GetNeighbors():
                            if nbr.GetIdx() not in ringSet and nbr.GetAtomicNum() == 7:
                                c1Idx = rIdx
                                break
                    if c1Idx is not None: break
                
                if c1Idx is not None:
                    c1Atom = mol.GetAtomWithIdx(c1Idx)
                    c2Idx = None
                    for nbr in c1Atom.GetNeighbors():
                        if nbr.GetIdx() in ringSet and nbr.GetAtomicNum() == 6:
                            c2Idx = nbr.GetIdx()
                            break
                    if c2Idx is not None:
                        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
                        c2Atom = mol.GetAtomWithIdx(c2Idx)
                        cip = c2Atom.GetProp("_CIPCode") if c2Atom.HasProp("_CIPCode") else "?"
                        if cip == "R":
                            nucleosideName = "D-Rib"
                        elif cip == "S":
                            nucleosideName = "D-Ara"
            
            anomer = baselineResult[1] if isinstance(baselineResult, tuple) else "?"
            if nucleosideName != "Pen":
                if modList:
                    shortMods = sorted(set(_MOD_LABEL_MAP.get(m, m) for m in modList))
                    if shortMods:
                        nucleosideName = f"{nucleosideName}({','.join(shortMods)})"
                return nucleosideName, anomer

    # ================================================================
    # Step 3.5: CIP+Exo 指纹引擎 (New Tier 1)
    # 当 v2 返回泛指标签时, 启动 CIP+Exo 引擎:
    #   1. 用 isolateBareSugarRing 切除苷元/相邻糖环, 只保留单独糖环
    #   2. 虚拟剥离所有修饰 → 裸糖
    #   3. 强制 CIP 重算 (Chem.AssignStereochemistry)
    #   4. 环行走提取 CIP+Exo 指纹
    #   5. 查表匹配
    #   6. 回填修饰标签
    # ================================================================
    if CIP_EXO_AVAILABLE:
        try:
            # Step 3.5.0: 在运行指纹提取前, 首先必须彻底切除所有苷元和相邻糖连接
            # Isolation is absolutely required so walkSugarRing doesn't get confused by macrocycles
            from lib.monosaccharide_identifier import isolateBareSugarRing
            bareSugarIsolated, adjustedRingAtoms = isolateBareSugarRing(mol, ringAtoms)
            
            # Step 3.5.1: 虚拟剥离 (Virtual Demodification)
            allDemodPatterns = list(DEMOD_PATTERNS.keys())
            bareMol, detectedMods = virtualDemodify(bareSugarIsolated, modNames=allDemodPatterns)

            # Step 3.5.2: 强制 CIP 重算 (Force CIP reassignment on bare sugar)
            try:
                Chem.AssignStereochemistry(bareMol, force=True, cleanIt=True)
            except Exception:
                pass

            # Step 3.5.3: 提取指纹 (Extract fingerprint)
            # CRITICAL: DO NOT PASS adjustedRingAtoms! bareMol atom indices have shifted!
            fp = extractSugarFingerprint(bareMol)

            if fp is not None:
                # Step 3.5.4: 查表匹配 (Fingerprint matching)
                matches = matchSugarFingerprint(fp)

                if matches:
                    topName, topAnomer, topScore = matches[0][0], matches[0][1], matches[0][2]

                    # N-guard: CIP+Exo 匹配到氨基糖时验证原始环上有 N
                    # N-guard: validate original isolated ring has N for amino sugar CIP match
                    _pureCipName = topName.split("(")[0] if "(" in topName else topName
                    if _pureCipName in _AMINO_SUGAR_NAMES_V10:
                        _ringSetCip = set(adjustedRingAtoms)
                        _hasExoNCip = any(
                            any(nbr.GetAtomicNum() == 7 and nbr.GetIdx() not in _ringSetCip
                                for nbr in bareSugarIsolated.GetAtomWithIdx(rIdx).GetNeighbors())
                            for rIdx in adjustedRingAtoms
                            if bareSugarIsolated.GetAtomWithIdx(rIdx).GetAtomicNum() == 6
                        )
                        if not _hasExoNCip:
                            topName = None  # 无效 → 跳过 CIP+Exo 结果

                    if topName is not None:
                        # Step 3.5.5: 回填修饰标签 (Reassemble labels)
                        # 检查是否有氨基修饰需要回填 NAc/NS 等
                        finalCipName = topName
                        if "N-Acetyl" in detectedMods:
                            if finalCipName.endswith("N"):
                                finalCipName = finalCipName + "Ac"
                            elif not finalCipName.endswith("NAc"):
                                finalCipName = finalCipName + "NAc"

                        # 附加其他修饰标签 (Append other mod labels)
                        otherMods = [m for m in detectedMods
                                     if m not in ("N-Acetyl",)]
                        if otherMods:
                            shortMods = sorted(set(
                                _MOD_LABEL_MAP.get(m, m) for m in otherMods
                            ))
                            if shortMods:
                                finalCipName = f"{finalCipName}({','.join(shortMods)})"

                        return finalCipName, topAnomer
        except Exception:
            pass  # CIP+Exo 失败 → 回退到 Step 4

    # ================================================================
    # Step 4: 基线是泛指标签 → 尝试剥离修饰后重新匹配 (升级路径)
    # Baseline is generic → try stripping mods and re-matching (escalation)
    # ================================================================
    if modList:
        try:
            taggedMol2 = Chem.RWMol(mol)
            for atom in taggedMol2.GetAtoms():
                atom.SetAtomMapNum(atom.GetIdx() + 1)
            taggedMol2 = taggedMol2.GetMol()
            ringMapNums = set(idx + 1 for idx in ringAtoms)

            strippedMol, _ = stripModificationsToBareSugar(taggedMol2, list(ringAtoms))

            # 通过 AtomMapNum 回映射找回糖环原子
            # Recover ring atoms via AtomMapNum
            workRingAtoms = []
            for atom in strippedMol.GetAtoms():
                if atom.GetAtomMapNum() in ringMapNums:
                    workRingAtoms.append(atom.GetIdx())

            if len(workRingAtoms) >= len(ringAtoms) - 1:
                # 清除 AtomMapNum
                cleanMol = Chem.RWMol(strippedMol)
                for atom in cleanMol.GetAtoms():
                    atom.SetAtomMapNum(0)
                try:
                    Chem.SanitizeMol(cleanMol)
                    Chem.AssignStereochemistry(cleanMol.GetMol(), force=True, cleanIt=True)
                    workingMol = cleanMol.GetMol()
                except Exception:
                    workingMol = cleanMol.GetMol()

                escalatedResult = identify_monosaccharide_v2(workingMol, workRingAtoms)
                escalatedName = escalatedResult[0] if isinstance(escalatedResult, tuple) else escalatedResult
                pureEscName = escalatedName.split("(")[0] if "(" in escalatedName else escalatedName

                # N-guard: 氨基糖升级结果需验证原始分子的环上有 N
                # N-guard: Amino sugar escalation result requires N on original ring
                if pureEscName in _AMINO_SUGAR_NAMES_V10:
                    ringSetOrig = set(ringAtoms)
                    hasExoNOrig = any(
                        any(nbr.GetAtomicNum() == 7 and nbr.GetIdx() not in ringSetOrig
                            for nbr in mol.GetAtomWithIdx(rIdx).GetNeighbors())
                        for rIdx in ringAtoms
                        if mol.GetAtomWithIdx(rIdx).GetAtomicNum() == 6
                    )
                    if not hasExoNOrig:
                        pureEscName = "Hex"  # 降级 → 不采用升级结果

                # 仅当升级结果更具体时才采用 (防止降级)
                # Only adopt if escalation is MORE specific (prevent degradation)
                if pureEscName not in GENERIC_LABELS:
                    shortMods = sorted(set(_MOD_LABEL_MAP.get(m, m) for m in modList))
                    if shortMods and not any(s in escalatedName for s in shortMods):
                        finalName = f"{escalatedName}({','.join(shortMods)})"
                        if isinstance(escalatedResult, tuple):
                            return finalName, escalatedResult[1]
                        else:
                            return finalName
                    return escalatedResult
        except Exception:
            pass

    # ================================================================
    # 兜底: 返回基线结果 + 修饰标签
    # Fallback: return baseline with mod labels
    # ================================================================
    if modList:
        shortMods = sorted(set(_MOD_LABEL_MAP.get(m, m) for m in modList))
        if shortMods and not any(s in baselineName for s in shortMods):
            finalName = f"{baselineName}({','.join(shortMods)})"
            if isinstance(baselineResult, tuple):
                return finalName, baselineResult[1]
            else:
                return finalName

    return baselineResult


# =====================================================================
# LLM 文献回溯占位函数 (LLM Literature Retrieval Placeholder)
# 用于未来接入 Data Mining LLM 模块, 通过文献自动纠正泛指糖标签
# Placeholder for future Data Mining LLM module integration
# =====================================================================
def requestLlmSugarIdentification(
    compoundId: str,
    doi: str = "",
    pubchemCid: str = "",
    genericLabel: str = ""
) -> str:
    """
    [占位函数 / Placeholder] 请求 LLM 模块从文献中识别具体糖单元。
    Request LLM module to identify specific sugars from literature.

    触发条件 (Trigger): 化合物含有模糊标签 (Hex, dHex, HexN, Pen, PenN)。
    工作流 (Workflow):
      1. 通过 PubChem CID 或 DOI 获取文献摘要
      2. 提问 LLM: "What are the specific sugar units (monosaccharides)
         present in this compound?"
      3. 解析 LLM 输出, 映射到标准糖名

    Args:
        compoundId: 化合物编号 (如 CNP0420046)
        doi: 文献 DOI (如 10.1021/...)
        pubchemCid: PubChem CID
        genericLabel: 当前模糊标签 (如 "Hex", "dHex")

    Returns:
        具体糖名 (如 "D-Glc") 或空字符串 (无法解析)

    TODO: 接入实际 LLM API 后实现
    """
    # [PLACEHOLDER] 当前返回空字符串, 表示无法解析
    # 实际实现时将调用 OpenAI/Claude API + PubChem/OpenAlex
    return ""


def get_linkage_type(mol, unit1_atoms, unit2_atoms):
    u1_set = set(unit1_atoms)
    u2_set = set(unit2_atoms)
    for idx1 in unit1_atoms:
        atom1 = mol.GetAtomWithIdx(idx1)
        for nbr in atom1.GetNeighbors():
            if nbr.GetIdx() in u1_set: continue
            if nbr.GetIdx() in u2_set: return "-"
            if nbr.GetSymbol() == "O":
                for nnbr in nbr.GetNeighbors():
                    if nnbr.GetIdx() == idx1: continue
                    if nnbr.GetIdx() in u2_set: return "-"
    return "/"

def _hasValidSugarRing(mol, ring_atoms):
    """
    拓扑铁律: 开环糖杀手。
    Strict Ring Enforcement: verify fragment has a valid 5 or 6-membered
    oxygen-containing heterocycle.

    返回 False → 该碎片是开环结构 (如糖醇), 应标记为 Non_Cyclic_Invalid。
    """
    from rdkit.Chem import rdMolDescriptors
    if rdMolDescriptors.CalcNumRings(mol) == 0:
        return False

    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if len(ring) not in (5, 6):
            continue
        # 检查是否含有至少一个氧原子 (含氧杂环)
        oxygenInRing = sum(
            1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8
        )
        if oxygenInRing >= 1:
            return True
    return False


def generate_refined_sequence(mol):
    if not mol: return "", ""
    # 使用 glycan_topology 获取糖单元
    # Use glycan_topology to get sugar units (was sugar_utils, which doesn't exist)
    sugar_units, atom_to_sugar = glycan_topology.get_sugar_units(mol)
    if not sugar_units: return "", ""

    # ---- 拓扑铁律: 开环糖过滤 (Ring Enforcement) ----
    # 标记不含有效含氧杂环的碎片为 Non_Cyclic_Invalid
    validUnits = []
    for unit in sugar_units:
        ringAtoms = unit.get('ring_atoms', [])
        if ringAtoms:
            validUnits.append(unit)
        else:
            # 如果 glycan_topology 没提供 ring_atoms, 仍保留
            validUnits.append(unit)

    sugar_units = validUnits
    if not sugar_units:
        return "Non_Cyclic_Invalid", ""

    G = nx.DiGraph()

    for i, unit in enumerate(sugar_units):
        unit['order_id'] = i + 1
        name = unit.get('name', 'Unknown')
        anomer = unit.get('anomeric_config', '?')
        mods = unit.get('modifications', [])
        
        G.add_node(unit['id'], name=name, anomer=anomer, mods=mods, order_id=unit['order_id'])
        
    linkages = glycan_topology.find_glycosidic_linkages(mol, sugar_units)
    
    for l in linkages:
        src = l['sugar_donor']
        dst = l['sugar_acceptor']
        link_pos = l['linkage']
        G.add_edge(src, dst, link=link_pos)
        
    roots = [n for n in G.nodes() if G.out_degree(n) == 0]
    
    if not roots:
        seq_parts = []
        mod_parts = []
        for n in G.nodes():
            nd = G.nodes[n]
            seq_parts.append(nd['name'])
            mod_str = ",".join(f"*{m}" for m in nd['mods']) if nd['mods'] else ""
            mod_parts.append(f"{nd['name']}_{nd['order_id']}({mod_str})")
        return "Cyclic", " ; ".join(mod_parts)
    
    def traverse_two_pass(node):
        nd = G.nodes[node]
        name = nd['name']
        mods = nd['mods']
        order_id = nd['order_id']
        
        mod_label = f"{name}_{order_id}({','.join(f'*{m}' for m in mods)})" if mods else f"{name}_{order_id}()"
        
        children = list(G.predecessors(node))
        children.sort()
        
        if not children:
            return name, mod_label
            
        main_child = children[0]
        branches = children[1:]
        
        m_link = G.edges[main_child, node].get('link', '?')
        m_anomer = G.nodes[main_child].get('anomer', '?')
        m_link_str = f"({m_anomer}{m_link})"
        
        m_seq, m_mod = traverse_two_pass(main_child)
        
        res_seq = f"{m_seq}-{m_link_str}-"
        res_mod = f"{m_mod}-"
        
        for b in branches:
            b_link = G.edges[b, node].get('link', '?')
            b_anomer = G.nodes[b].get('anomer', '?')
            b_link_str = f"({b_anomer}{b_link})"
            
            b_seq, b_mod = traverse_two_pass(b)
            res_seq = f"[{b_seq}-{b_link_str}]-" + res_seq
            res_mod = f"[{b_mod}]-" + res_mod
            
        res_seq += name
        res_mod += mod_label
        
        return res_seq, res_mod

    seq_list = []
    mod_list = []
    for root in roots:
        s, m = traverse_two_pass(root)
        seq_list.append(s)
        mod_list.append(m)
        
    return " ; ".join(seq_list), " ; ".join(mod_list)

def analyze_glycan(smiles):
    if not smiles or pd.isna(smiles) or smiles == "nan": return "", ""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if not mol: return "Invalid", ""
        return generate_refined_sequence(mol)
    except Exception as e:
        return "Error", ""
