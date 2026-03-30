from lib.glycan_reference_library import (
    NUCLEOBASE_LIBRARY, 
    PHOSPHATE_SMARTS,
    AMINO_ACID_LIBRARY,
    STRICT_PEPTIDE_BOND_SMARTS,
    STRIPPING_REACTIONS
)
from rdkit import Chem


def strip_and_record_modifications(sugar_mol):
    """
    Identifies and records modifications on the sugar ring by matching the reactant side
    of STRIPPING_REACTIONS. It returns a list of modifications found.
    We don't actually 'strip' the atoms out of the RDKit molecule because that ruins 
    atom indexing required for DAG linkage. We instead identify the modifications 
    attached to the sugar ring atoms.
    """
    modifications = []
    
    for mod_name, pat in STRIPPING_REACTIONS.items():
        if not pat: continue
        matches = sugar_mol.GetSubstructMatches(pat)
        for match in matches:
            # We don't strictly care which atom it is attached to for the label right now,
            # just that the modification exists on this sugar molecule.
            modifications.append(mod_name)
            
    return modifications

def _classifyAcylGroup(mol, carbonylCarbonIdx, esterOxygenIdx):
    """
    从酯键/酰胺键的羰基碳出发, BFS 沿 C-C 键计数碳原子数 → 精确分类酰基。
    Classify acyl group by counting carbons via C-C BFS from carbonyl carbon.
    Heteroatoms (O, N, S, P) act as hard barriers — only contiguous C-C bonds
    are followed.

    设计意图: 区分乙酰基 (2C), 丙酰 (3C), 丁酰 (4C), 异戊酰 (5C), 己酰 (6C)
    Design: distinguish Ac(2C), Propionyl(3C), Butyryl(4C), Isovaleryl(5C), Hexanoyl(6C)

    Args:
        mol: RDKit 分子对象
        carbonylCarbonIdx: 羰基碳的原子索引 (atom index of carbonyl carbon C(=O))
        esterOxygenIdx: 酯氧的原子索引 (atom index of ester oxygen, used as barrier)

    Returns:
        str: 修饰标签 (e.g. "O-Ac", "Propionyl", "Butyryl", "Isovaleryl", "Hexanoyl")
    """
    visited = {esterOxygenIdx}  # 酯氧作为屏障, 不向回走 (ester O acts as barrier)
    queue = [carbonylCarbonIdx]
    visited.add(carbonylCarbonIdx)
    carbonCount = 0
    aromaticCount = 0  # Bug 1 修复: 统计芳香碳数 (Count aromatic carbons)

    while queue:
        idx = queue.pop(0)
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6:
            carbonCount += 1
            if atom.GetIsAromatic():
                aromaticCount += 1
        for nbr in atom.GetNeighbors():
            nIdx = nbr.GetIdx()
            if nIdx in visited:
                continue
            visited.add(nIdx)
            # 杂原子截断铁律: O(8), N(7), S(16), P(15) → 硬停
            # Heteroatom barrier: stop at O, N, S, P (don't follow)
            if nbr.GetAtomicNum() != 6:
                continue
            queue.append(nIdx)

    # Bug 1 修复: 芳香性优先判断 — 芳香碳 ≥ 5 则为芳香酰基
    # Bug 1 fix: Aromaticity-first classification.
    # Benzoyl (C6H5CO-) has 6 aromatic C + 1 carbonyl C = 7C total.
    # Galloyl has similar aromatic ring. Must NOT be classified as Heptanoyl.
    if aromaticCount >= 5:
        # 已知芳香酰基映射 (Known aromatic acyl mapping)
        # 设计意图: 区分苯甲酰基 (7C), 没食子酰基 (7C+3OH), 桂皮酰基 (9C) 等
        # Design: distinguish Benzoyl (7C), Galloyl (7C+3OH), cinnamoyl (9C) etc.
        if carbonCount <= 7:
            return "Bz"          # 苯甲酰基 (Benzoyl, C₆H₅CO-)
        elif carbonCount == 9:
            return "pCou"        # 对香豆酰基 (p-Coumaroyl, C₆H₄(OH)CH=CHCO-)
        elif carbonCount == 10:
            return "Fer"         # 阿魏酰基 (Feruloyl, 10C)
        elif carbonCount == 11:
            return "Sin"         # 芥子酰基 (Sinapoyl, 11C)
        return f"Aroyl-{carbonCount}C"  # 未知芳香酰基 (Unknown aromatic acyl)

    # 脂肪酰基: 碳数 → 酰基名称映射 (Aliphatic acyl: carbon count → name)
    ACYL_NAMES = {
        2: "O-Ac",          # 乙酰基 (Acetyl, CH₃CO-)
        3: "Propionyl",     # 丙酰基 (Propionyl, CH₃CH₂CO-)
        4: "Butyryl",       # 丁酰基 (Butyryl, CH₃(CH₂)₂CO-)
        5: "Tig",           # 当归酰/替格酰基 (Tigloyl/Angeloyl, common in NP)
        6: "Hexanoyl",      # 己酰基 (Hexanoyl)
        7: "Heptanoyl",     # 庚酰基 (Heptanoyl) — 仅脂肪链才到达这里
    }
    return ACYL_NAMES.get(carbonCount, f"Acyl-{carbonCount}C")

def find_mapped_sugar_units(mol):
    """
    1. Finds all potential sugar rings topologically.
    2. Identifies C1 -> C_n and maps positions.
    3. Uses lib.monosaccharide_identifier to identify basic sugar names.
    4. Extracts unit modifications.
    """
    ri = mol.GetRingInfo()
    matched_units = []
    
    potential_rings = []
    for ring in ri.AtomRings():
        if len(ring) in (5, 6):
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if [a.GetSymbol() for a in atoms].count("O") == 1:
                potential_rings.append(list(ring))

    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    conf_map = {idx: conf for idx, conf in centers}

    nucleo_pats = list(NUCLEOBASE_LIBRARY.values())
    
    valid_rings = []
    ring_set_cache = set()
    for ring in potential_rings:
        ring_set = set(ring)
        is_fused = any(ri.NumAtomRings(idx) > 1 for idx in ring)
        if is_fused:
            # 颁布"桥连/宏环糖"拓扑豁免法案 (Stapled Sugar Exemption)
            # 评估所有与当前糖环共享原子的外部环 (Evaluate exemptions for all externally shared rings)
            all_shared_rings_exempt = True
            for other_ring in ri.AtomRings():
                other_ring_set = set(other_ring)
                if other_ring_set == ring_set:
                    continue
                
                # 如果这个外部环和我们共享了原子 (If they share atoms: Fused/Bridged/Spiro)
                shared_atoms = ring_set.intersection(other_ring_set)
                if not shared_atoms:
                    continue
                    
                # === 拓扑铁律：多环含氧聚醚检测 (Polycyclic Oxygenated Polyether Detection) ===
                # 如果两个含氧环共享一条边 (≥2 原子)，且外部环也是含氧五/六元环，
                # 说明这极其可能是聚醚类有机多环体系 (如莫能菌素)，绝对不是糖！
                # If two oxygen-containing rings share an edge (≥2 atoms), and the other
                # ring is also a 5/6-membered oxygenated ring, this is a polyether, NOT sugar.
                if len(shared_atoms) >= 2 and len(other_ring) in (5, 6):
                    other_symbols = [mol.GetAtomWithIdx(i).GetSymbol() for i in other_ring]
                    if other_symbols.count('O') >= 1:
                        # 两个含氧环共享一条边 → 聚醚拓扑，强制拒绝
                        # Two oxygen-containing rings sharing an edge → polyether topology, reject
                        all_shared_rings_exempt = False
                        break

                # 豁免条件 A (宏环豁免): 外部环是一个大环 (环原子数 >= 8)，说明这是大环交联绝对允许
                is_macrocycle = len(other_ring) >= 8
                
                # 豁免条件 B (缩醛/脱水桥豁免): 外部环是一个纯由碳和氧组成的小环 (如 1,6-脱水桥或 4,6-缩酮) 且必须包含氧！
                # 且必须只共享 1 个原子 (Spiro junction, 如螺环缩酮)
                other_atoms = [mol.GetAtomWithIdx(i) for i in other_ring]
                symbols = [a.GetSymbol() for a in other_atoms]
                is_acetal_bridge = all(sym in ('C', 'O') for sym in symbols) and ('O' in symbols)
                
                if not (is_macrocycle or is_acetal_bridge):
                    all_shared_rings_exempt = False
                    break
                        
            if not all_shared_rings_exempt:
                continue

        # === 不饱和环检测 (Unsaturated Ring Bond Detection) ===
        # 仅检查 **环内键** (ring-internal bonds) 是否为 SINGLE。
        # 环外键 (如 C=O 羧基, C=O 乙酰基) 不影响糖环本身的饱和度判定。
        # Only check bonds BETWEEN ring atoms for saturation. Exocyclic bonds
        # (e.g., C=O of uronic acids, acetyl groups) do NOT disqualify a sugar ring.
        has_unsaturated_ring_bond = False
        for i, idx1 in enumerate(ring):
            for idx2 in ring[i+1:]:
                bond = mol.GetBondBetweenAtoms(idx1, idx2)
                if bond and bond.GetBondType() != Chem.BondType.SINGLE:
                    has_unsaturated_ring_bond = True
                    break
            if has_unsaturated_ring_bond:
                break
        if has_unsaturated_ring_bond:
            continue

        # === 内酯/酯环检测 (Lactone / Ester Ring Detection) ===
        # Bug1 修复: 糖环上的碳原子必须是 SP3。
        # 如果环上的碳有环外双键 (C=O, C=N), 则该碳是 SP2 → 内酯/酯环, 非糖。
        # Bug1 fix: Sugar ring carbons must be SP3. If any ring carbon has an
        # exocyclic double bond (C=O, C=N → SP2 hybridization), this is a
        # lactone/ester ring, NOT a sugar ring.
        has_sp2_ring_carbon = False
        ring_set = set(ring)
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() == 'C':
                # 检查该碳是否有环外不饱和键
                # Check if this ring carbon has exocyclic unsaturated bond
                for bond in atom.GetBonds():
                    otherIdx = bond.GetOtherAtomIdx(idx)
                    if otherIdx not in ring_set:
                        if bond.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE):
                            has_sp2_ring_carbon = True
                            break
                if has_sp2_ring_carbon:
                    break
        if has_sp2_ring_carbon:
            continue
            
        is_nucleoside = False
        if len(ring) == 5:
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() not in ring and nbr.GetSymbol() == 'N':
                        for npat in nucleo_pats:
                            if npat and mol.HasSubstructMatch(npat):
                                for b_match in mol.GetSubstructMatches(npat):
                                    if nbr.GetIdx() in b_match:
                                        is_nucleoside = True; break
                            if is_nucleoside: break
                    if is_nucleoside: break
                if is_nucleoside: break
        # 2026-03-25: 不再排除核苷糖环 — 呋喃糖鉴定引擎已就绪
        # No longer exclude nucleoside sugar rings — furanose engine is ready
        # 仅标记, 不跳过 (Label only, do not skip)
        pass

        # Topo oxygenation check: standard sugar requires >= 2 valid exocyclic oxygen/heteroatoms
        # Typical glycosides / pure sugars will easily pass this. 
        # C-O (like C6) also counts as a valid substituent path.
        exo_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() == 'C':
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() not in ring:
                        sym = nbr.GetSymbol()
                        if sym in ('O', 'N', 'S', 'P'):
                            exo_count += 1
                        elif sym == 'C':
                            for nnbr in nbr.GetNeighbors():
                                if nnbr.GetSymbol() in ('O', 'N', 'S', 'P') and nnbr.GetIdx() != atom.GetIdx():
                                    exo_count += 1
        if exo_count < 2:
            continue

        # === 羟基密度门控 (Hydroxyl Density Gate) ===
        # 真正的糖环碳上必须有足够的 -OH 基团 (如 Glc=4, Xyl=3, dHex=3)
        # 非糖含氧杂环 (如四氢吡喃类) 通常 -OH 很少
        # True sugar ring carbons must have sufficient -OH groups.
        # Non-sugar oxacycles (tetrahydropyrans, cyclitol derivatives) typically have fewer.
        hydroxylCount = 0
        ring_set_local = set(ring)
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != 'C':
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set_local:
                    continue
                # 直接 -OH: O 原子且度为 1 (仅连一个重原子) 或度为 2 (含 H)
                # Direct -OH: oxygen neighbor outside ring with degree ≤ 2 (terminal or OH)
                if nbr.GetSymbol() == 'O' and nbr.GetDegree() <= 2:
                    # 排除酯键和缩醛键的 O (这些 O 度通常为 2 且连着 C=O 或另一个 C)
                    # Exclude ester/acetal O (connected to C=O or another C)
                    isTerminalOH = True
                    for onnbr in nbr.GetNeighbors():
                        if onnbr.GetIdx() == idx:
                            continue
                        # 如果 O 连着另一个碳 → 可能是糖苷键 O、甲氧基等，仍算有 O
                        # If O connects to another C → glycosidic O, still valid
                    hydroxylCount += 1
                elif nbr.GetSymbol() == 'N':
                    # 氨基也算作等价的极性取代基 (Amino counts as equivalent polar substituent)
                    hydroxylCount += 1
        # 呋喃糖放宽: 5元环 deoxyribose 仅 1 个 OH, 阈值降为 1
        # Furanose relaxation: deoxyribose has only 1 OH, lower threshold to 1
        minOH = 1 if len(ring) == 5 else 2
        if hydroxylCount < minOH:
            continue

        # === 手性中心门控 (Chiral Center Gate) ===
        # 真正的糖环有多个手性碳 (Glc=4, Xyl=3, 最少的如 glyceraldehyde 也有 1)
        # 非糖杂环通常 0-1 个手性中心
        # True sugar rings have multiple chiral carbons (Glc=4, Xyl=3).
        # Non-sugar heterocycles typically have 0-1 chiral centers.
        chiralRingCarbons = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != 'C':
                continue
            # 手性碳代理指标: 连接 ≥3 个不同的重原子邻居 (含环内 + 环外)
            # Chiral carbon proxy: ≥3 distinct heavy-atom neighbors (in-ring + exocyclic)
            heavyNeighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() > 1]
            if len(heavyNeighbors) >= 3:
                # 进一步检查: 至少有一个环外杂原子 (O/N/S) 邻居
                # Further check: at least 1 exocyclic heteroatom neighbor
                hasExoHetero = any(
                    n.GetIdx() not in ring_set_local and n.GetAtomicNum() in (7, 8, 16)
                    for n in heavyNeighbors
                )
                if hasExoHetero:
                    chiralRingCarbons += 1
        if chiralRingCarbons < 2:
            continue

        # === 聚醚/大环排除: 环外长碳链检测 (Polyether Exclusion: Exocyclic C-chain) ===
        # 真正的糖环碳最多连 CH2OH (2 个重原子的短链)
        # 聚醚/大环内酯在环碳上有长碳链延伸 (C-C-C-C...)
        # True sugar ring C has at most CH2OH (2 heavy atom short chain).
        # Polyether/macrolides have long carbon chains extending from ring C.
        # Gate: if any ring C → exo C → follow C-C chain ≥ 4 → not sugar
        hasLongCarbonChain = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != 'C':
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set_local:
                    continue
                # 仅追踪 sp3 非芳香碳 (排除 C-糖苷的芳香链)
                # Only follow sp3 non-aromatic C (exclude C-glycoside aromatic chain)
                if nbr.GetSymbol() != 'C' or nbr.GetIsAromatic():
                    continue
                # 从这个环外 C 开始, 沿纯 sp3 碳链走 (BFS)
                # Walk pure sp3 carbon chain from this exocyclic C
                chainLen = 1
                visited = {idx, nbr.GetIdx()}
                frontier = [nbr.GetIdx()]
                while frontier and chainLen < 4:
                    nextFrontier = []
                    for fIdx in frontier:
                        fAtom = mol.GetAtomWithIdx(fIdx)
                        for fNbr in fAtom.GetNeighbors():
                            if fNbr.GetIdx() in visited:
                                continue
                            # 只跟踪非芳香 C (Only follow non-aromatic C)
                            if fNbr.GetSymbol() == 'C' and not fNbr.GetIsAromatic():
                                chainLen += 1
                                visited.add(fNbr.GetIdx())
                                nextFrontier.append(fNbr.GetIdx())
                    frontier = nextFrontier
                if chainLen >= 4:
                    hasLongCarbonChain = True
                    break
            if hasLongCarbonChain:
                break
        if hasLongCarbonChain:
            continue

        # === 糖键连方式约束 (Glycosidic Linkage Type Constraint) ===
        # 真正的糖环与外部骨架仅通过以下桥原子类型连接:
        # C-O-C, C-N-C, C-S-C, C-C (直接碳碳键)
        # 禁止异常连接如 C-P, C-B, C-Si 等无化学意义的键
        # Real sugar rings connect to external structures only through: 
        # C-O-C, C-N-C, C-S-C, or direct C-C bonds.
        has_invalid_linkage = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != 'C':
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set:
                    continue
                sym = nbr.GetSymbol()
                # 允许 O, N, S, C, H (以及 P 在磷酸酯情况下也允许)
                # Allow O, N, S, C, H (and P for phosphate esters)
                if sym not in ('O', 'N', 'S', 'C', 'H', 'P'):
                    has_invalid_linkage = True
                    break
            if has_invalid_linkage:
                break
        if has_invalid_linkage:
            continue
            
        # === 大环共碳检测 (Macrolide Co-Ring Detection) ===
        # 如果环中任何碳原子同时属于 7+ 元大环 → 这不是独立糖环
        # 而是大环内酯 (如 Amphotericin B) 折叠形成的假糖环
        # If any ring carbon also belongs to a 7+ membered ring → not an
        # independent sugar ring, but a macrolide fold artifact
        isMacrolideFold = False
        for idx in ring:
            if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6:
                continue
            for otherRing in ri.AtomRings():
                if len(otherRing) >= 7 and idx in otherRing:
                    isMacrolideFold = True
                    break
            if isMacrolideFold:
                break
        if isMacrolideFold:
            continue

        valid_rings.append(ring)
        
    import lib.monosaccharide_identifier as seq
        
    for sid, ring in enumerate(valid_rings):
        ring_o_idx = None
        for idx in ring:
            if mol.GetAtomWithIdx(idx).GetSymbol() == 'O':
                ring_o_idx = idx
                break
                
        o_atom = mol.GetAtomWithIdx(ring_o_idx)
        ring_neighbors = [n.GetIdx() for n in o_atom.GetNeighbors() if n.GetIdx() in ring]
        
        if len(ring_neighbors) != 2: continue
        n1, n2 = ring_neighbors
        
        c1_idx = None
        for n in ring_neighbors:
            atom = mol.GetAtomWithIdx(n)
            oxygens_attached = 0
            for nn in atom.GetNeighbors():
                if nn.GetSymbol() == 'O':
                    oxygens_attached += 1
            if oxygens_attached >= 2:
                c1_idx = n
                break
                
        if c1_idx is None:
            # Fallback if no anomeric oxygen is found (e.g. reduced sugar)
            c1_idx = n1
            
        path = [c1_idx]
        curr = c1_idx
        prev = ring_o_idx
        for _ in range(len(ring)-2):
            for n in mol.GetAtomWithIdx(curr).GetNeighbors():
                n_idx = n.GetIdx()
                if n_idx in ring and n_idx != prev:
                    prev = curr
                    curr = n_idx
                    path.append(curr)
                    break
                    
        pos_map = {}
        is_ketose = False
        c1_exo_c = None
        c1_atom = mol.GetAtomWithIdx(c1_idx)
        for nbr in c1_atom.GetNeighbors():
            if nbr.GetIdx() not in ring and nbr.GetSymbol() == 'C':
                is_ketose = True
                c1_exo_c = nbr.GetIdx()
                break
                
        offset = 2 if is_ketose else 1
        for i, idx in enumerate(path):
            pos_map[idx] = i + offset
            
        if is_ketose:
            pos_map[c1_exo_c] = 1
            
        oxygen_map = {}
        c6_idx = None
        for idx, pos in list(pos_map.items()):
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx not in ring:
                    sym = nbr.GetSymbol()
                    if sym in ('O', 'N', 'S'):
                        oxygen_map[nbr_idx] = pos
                    elif sym == 'C' and pos == len(ring) - 1:
                        c6_idx = nbr_idx
                        pos_map[c6_idx] = pos + 1
                        for nnbr in nbr.GetNeighbors():
                            if nnbr.GetSymbol() in ('O', 'N', 'S') and nnbr.GetIdx() != idx:
                                oxygen_map[nnbr.GetIdx()] = pos + 1
                                break
                                
        # Extract stereochemistry strictly on the isolated monomer to avoid polymer-induced CIP inversion
        anomeric_idx = c1_idx
        ref_idx = path[-1]
            
        ring_set = frozenset(ring)
        if any(frozenset(ex['ring_atoms']) == ring_set for ex in matched_units):
            continue
            
        # 接收基于拓扑匹配的绝对结果
        identify_result = seq.identify_monosaccharide_v10(mol, ring)
        if isinstance(identify_result, tuple):
            base_name, anomer_config = identify_result
        else:
            base_name = identify_result
            anomer_config = "?"

        # === CIP-based α/β 异头构型强制检测 (Forced Anomeric Config via CIP) ===
        # =================================================================
        # 设计修正 (2026-03-18): 字典模板已在加载时剥离异头碳手性 →
        # identify_monosaccharide_v2 永远返回 "a" → 旧条件 "if == ?" 永远不触发!
        # 修复: 始终使用 CIP 推算 α/β, 字典返回值仅作为 fallback。
        #
        # Design Fix: Dictionary templates have anomeric chirality stripped at
        # load time → v2 always returns "a" → old gate "if == ?" never fired!
        # Fix: ALWAYS compute α/β from CIP codes. Dict value is unreliable.
        #
        # 哈沃斯法则 (Haworth Convention):
        #   D-糖: C1 CIP == Cref CIP → α (axial); C1 CIP != Cref CIP → β (equatorial)
        #   L-糖: 规则相反 (reversed)
        # =================================================================
        # 手术 1: 断口清洗 — 清除虚拟原子对 CIP 的污染
        # Surgery 1: Dummy atom cleanup — prevent CIP poisoning from cleavage scars
        try:
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
            c1_atom_chiral = mol.GetAtomWithIdx(anomeric_idx)
            cref_atom_chiral = mol.GetAtomWithIdx(ref_idx)
            c1_cip = c1_atom_chiral.GetProp('_CIPCode') if c1_atom_chiral.HasProp('_CIPCode') else None
            cref_cip = cref_atom_chiral.GetProp('_CIPCode') if cref_atom_chiral.HasProp('_CIPCode') else None
            if c1_cip and cref_cip:
                isLSugar = base_name.startswith("L-")
                if isLSugar:
                    # L-系列: C1 与 Cref CIP 不同 → α; 相同 → β
                    # L-sugar: C1 CIP != Cref CIP → α; same → β
                    anomer_config = "a" if c1_cip != cref_cip else "b"
                else:
                    # D-系列及其他: C1 与 Cref CIP 相同 → α; 不同 → β
                    # D-sugar: C1 CIP == Cref CIP → α; different → β
                    anomer_config = "a" if c1_cip == cref_cip else "b"
        except Exception:
            pass  # 保留字典返回值 (Keep dict value as fallback)
            
        # === 清除遗留修饰标签 (Strip Legacy Modification Labels) ===
        # identify_monosaccharide_v2 中的 check_modifications() 使用全分子 SMARTS 匹配
        # 可能将其他位置的修饰误加给当前糖环. v2.0 直接扫描器是唯一权威.
        # check_modifications() in identify_monosaccharide_v2 uses whole-molecule SMARTS
        # matching which can misattribute modifications from other parts of the molecule.
        # The v2.0 direct scanner is now the sole authority for modification detection.
        if "(" in base_name:
            name = base_name.split("(")[0]
        else:
            name = base_name
        seq_mods = []  # 清空遗留, 由 v2.0 直接扫描重新填充
                       # Clear legacy mods; v2.0 direct scanner will repopulate

        # === 直接修饰基团扫描 v2.0 (Direct Modification Scanning v2.0) ===
        # 扫描糖环碳的外环取代基: O-Me, O-Acyl(精确碳数), Sulfate, Phosphate, NAc, NH2
        # Scan exocyclic substituents: O-Me, O-Acyl(precise carbon count),
        # Sulfate (2-hop), Phosphate (2-hop), NAc, NH2
        ring_set_scan = set(ring)
        for rIdx in ring:
            rAtom = mol.GetAtomWithIdx(rIdx)
            if rAtom.GetSymbol() != 'C':
                continue
            for exoNbr in rAtom.GetNeighbors():
                exoIdx = exoNbr.GetIdx()
                if exoIdx in ring_set_scan:
                    continue
                exoSym = exoNbr.GetSymbol()

                if exoSym == 'O':
                    for oNbr in exoNbr.GetNeighbors():
                        if oNbr.GetIdx() == rIdx:
                            continue

                        # --- O-CH₃ (O-Methyl) ---
                        if (oNbr.GetSymbol() == 'C'
                            and oNbr.GetDegree() == 1
                            and oNbr.GetTotalNumHs() == 3):
                            if "O-Me" not in seq_mods:
                                seq_mods.append("O-Me")

                        # --- O-Acyl (酯键酰基, 精确碳数分类) ---
                        # O-C(=O)-R: 检查羰基, 然后 BFS 计碳数
                        # O-Acyl (ester): detect carbonyl, then BFS count carbons
                        elif oNbr.GetSymbol() == 'C' and oNbr.GetDegree() >= 2:
                            hasCarbonyl = any(
                                mol.GetBondBetweenAtoms(oNbr.GetIdx(), acNbr.GetIdx()).GetBondTypeAsDouble() == 2.0
                                for acNbr in oNbr.GetNeighbors()
                                if acNbr.GetSymbol() == 'O' and acNbr.GetIdx() != exoIdx
                                and mol.GetBondBetweenAtoms(oNbr.GetIdx(), acNbr.GetIdx()) is not None
                            )
                            if hasCarbonyl:
                                acylLabel = _classifyAcylGroup(mol, oNbr.GetIdx(), exoIdx)
                                if acylLabel not in seq_mods:
                                    seq_mods.append(acylLabel)

                        # --- Sulfate: O-S(=O)(=O)-OH (二级邻居检测) ---
                        # Sulfate: Ring-C → O → S(=O)₂(OH) (two-hop detection)
                        elif oNbr.GetSymbol() == 'S':
                            doubleBondO = 0
                            for sNbr in oNbr.GetNeighbors():
                                if sNbr.GetIdx() == exoIdx:
                                    continue
                                sBond = mol.GetBondBetweenAtoms(oNbr.GetIdx(), sNbr.GetIdx())
                                if sNbr.GetSymbol() == 'O' and sBond:
                                    if sBond.GetBondTypeAsDouble() == 2.0:
                                        doubleBondO += 1
                            if doubleBondO >= 2:
                                if "Sulfate" not in seq_mods:
                                    seq_mods.append("Sulfate")

                        # --- Phosphate: O-P(=O)(OH)(OH) (二级邻居检测) ---
                        # Phosphate: Ring-C → O → P(=O)(OH)₂ (two-hop detection)
                        elif oNbr.GetSymbol() == 'P':
                            doubleBondO_P = 0
                            for pNbr in oNbr.GetNeighbors():
                                if pNbr.GetIdx() == exoIdx:
                                    continue
                                pBond = mol.GetBondBetweenAtoms(oNbr.GetIdx(), pNbr.GetIdx())
                                if pNbr.GetSymbol() == 'O' and pBond:
                                    if pBond.GetBondTypeAsDouble() == 2.0:
                                        doubleBondO_P += 1
                            if doubleBondO_P >= 1:
                                if "Phosphate" not in seq_mods:
                                    seq_mods.append("Phosphate")

                elif exoSym == 'N':
                    # 检查 NAc: N-C(=O)-CH₃
                    # Check for N-Acetyl: N-C(=O)-CH₃
                    isNAc = False
                    for nNbr in exoNbr.GetNeighbors():
                        if nNbr.GetIdx() == rIdx:
                            continue
                        if nNbr.GetSymbol() == 'C':
                            hasCarbonylN = False
                            hasMethylN = False
                            for nacNbr in nNbr.GetNeighbors():
                                if nacNbr.GetIdx() == exoIdx:
                                    continue
                                bond = mol.GetBondBetweenAtoms(nNbr.GetIdx(), nacNbr.GetIdx())
                                if nacNbr.GetSymbol() == 'O' and bond and bond.GetBondTypeAsDouble() == 2.0:
                                    hasCarbonylN = True
                                if (nacNbr.GetSymbol() == 'C'
                                    and nacNbr.GetDegree() == 1
                                    and nacNbr.GetTotalNumHs() == 3):
                                    hasMethylN = True
                            if hasCarbonylN and hasMethylN:
                                isNAc = True
                    if isNAc and "NAc" not in seq_mods:
                        seq_mods.append("NAc")
                    elif not isNAc:
                        # 游离氨基 (Free amino): N 未被乙酰化
                        if exoNbr.GetTotalNumHs() >= 2 and "NH2" not in seq_mods:
                            seq_mods.append("NH2")

        final_mods = sorted(list(set(seq_mods)))
            
        matched_units.append({
            'name': name,
            'ring_atoms': list(ring),
            'ring_oxygen': ring_o_idx,
            'anomeric_idx': anomeric_idx,
            'anomeric_config': anomer_config,
            'position_map': pos_map,
            'oxygen_map': oxygen_map,
            'modifications': final_mods
        })
        
    return matched_units

def get_sugar_units(mol):
    units_info = find_mapped_sugar_units(mol)
    sugar_units = []
    atom_to_sugar = {}
    
    for sid, unit in enumerate(units_info):
        unit['id'] = sid
        sugar_units.append(unit)
        for a_idx in unit['position_map'].keys():
            atom_to_sugar[a_idx] = sid
        for a_idx in unit['oxygen_map'].keys():
            atom_to_sugar[a_idx] = sid
            
    return sugar_units, atom_to_sugar


def find_nucleotides(mol):
    """
    Returns a set of atom indices that belong to an intact Nucleotide 
    (Ribose/Deoxyribose + Nucleobase + Phosphate group).
    It colors the sugar ring, the base, and the phosphate.
    """
    nucleotide_atoms = set()
    ri = mol.GetRingInfo()
    
    # Pre-find 5-membered oxygenated rings (Ribose candidates)
    ribose_rings = []
    for ring in ri.AtomRings():
        if len(ring) == 5:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if [a.GetSymbol() for a in atoms].count("O") == 1:
                ribose_rings.append(list(ring))

    nucleo_pats = list(NUCLEOBASE_LIBRARY.values())
    
    # 1. Has Phosphate?
    has_phosphate = False
    phosphate_atoms = set()
    if PHOSPHATE_SMARTS and mol.HasSubstructMatch(PHOSPHATE_SMARTS):
        matches = mol.GetSubstructMatches(PHOSPHATE_SMARTS)
        for match in matches:
            phosphate_atoms.update(match)
            has_phosphate = True
            
    if not has_phosphate:
        return nucleotide_atoms # Strict validation: Must have phosphate to be a Nucleotide
    
    for ring in ribose_rings:
        is_nucleoside = False
        matched_base_atoms = set()
        
        # 2. Check for Nucleobase connection
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in ring and nbr.GetSymbol() == 'N':
                    for pat in nucleo_pats:
                        if pat and mol.HasSubstructMatch(pat):
                            matches = mol.GetSubstructMatches(pat)
                            for match in matches:
                                if nbr.GetIdx() in match:
                                    is_nucleoside = True
                                    matched_base_atoms.update(match)
                                    break
                        if is_nucleoside: break
                if is_nucleoside: break
            if is_nucleoside: break
            
        if is_nucleoside:
            # 3. Verify it connects to the Phosphate via C5
            connected_to_phosphate = False
            extended_atoms = set(ring).union(matched_base_atoms)
            
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() not in ring and nbr.GetSymbol() == "C" and nbr.GetHybridization() == Chem.HybridizationType.SP3:
                        # Exocyclic C5. Check if C5 is bonded to an O that is in the phosphate group
                        extended_atoms.add(nbr.GetIdx())
                        for c5_nbr in nbr.GetNeighbors():
                            if c5_nbr.GetIdx() in phosphate_atoms:
                                connected_to_phosphate = True
                                extended_atoms.update(phosphate_atoms)
                                
            if connected_to_phosphate:
                nucleotide_atoms.update(extended_atoms)
                        
    return nucleotide_atoms
                        
    return nucleoside_atoms


# ---------- 3. Build sugar units ----------
# This section is now handled by the new get_sugar_units above.

# ---------- 4. Find glycosidic oxygens (sugar–sugar bridge) ----------

def is_bridge_oxygen_between_sugars(mol, atom, atom_to_sugar):
    """
    Check if:
      - atom is O
      - it is not in a ring
      - it has 2 neighbors (C–C) with single bonds
      - and both C neighbors belong to sugars (atom_to_sugar)
    """
    if atom.GetSymbol() != "O":
        return False
    if atom.IsInRing():
        return False

    bonds = list(atom.GetBonds())
    if len(bonds) != 2:
        return False

    neighbors = []
    for b in bonds:
        if b.GetBondType() != Chem.BondType.SINGLE:
            return False
        other = b.GetOtherAtom(atom)
        if other.GetSymbol() != "C":
            return False
        neighbors.append(other)

    # must be attached to TWO sugar carbons
    sugar_ids = []
    for n in neighbors:
        sid = atom_to_sugar.get(n.GetIdx(), None)
        sugar_ids.append(sid)

    if sugar_ids[0] is None or sugar_ids[1] is None:
        return False

    # if both are sugars, treat as a bridge
    return True


# The original find_glycosidic_linkages is removed as per instruction.
# The new implementation of sugar unit detection (find_mapped_sugar_units)
# provides position_map directly, making the old numbering function obsolete.
# A new linkage detection function would be needed if this functionality is still desired.


# ---------- 4. Find glycosidic oxygens (sugar–sugar bridge) ----------

def is_anomeric(mol, idx):
    """判断一个碳原子是否为异头碳（连接>=2个氧原子）"""
    atom = mol.GetAtomWithIdx(idx)
    return sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 8) >= 2

def find_glycosidic_linkages(mol, units):
    linkages = []
    # 1. 构建扩展环 (包含环外 C6)
    for u in units:
        ring = set(u['ring_atoms'])
        exo_carbons = set()
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetSymbol() == 'C' and nbr.GetIdx() not in ring:
                    exo_carbons.add(nbr.GetIdx())
        u['extended_ring'] = ring.union(exo_carbons)
        
    # 2. 严格的单向寻径
    for i, u1 in enumerate(units):
        for j, u2 in enumerate(units):
            if i == j: continue
            
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == 'O':
                    nbrs = [n.GetIdx() for n in atom.GetNeighbors()]
                    if len(nbrs) == 2:
                        n0, n1 = nbrs[0], nbrs[1]
                        donor_idx, acceptor_idx = None, None
                        
                        # 绝对法则：u1 必须是提供异头碳(Donor)的一方，u2 必须提供受体碳(Acceptor)
                        if (n0 in u1['ring_atoms'] and is_anomeric(mol, n0)) and (n1 in u2['extended_ring'] and n1 not in u1['ring_atoms']):
                            donor_idx, acceptor_idx = n0, n1
                        elif (n1 in u1['ring_atoms'] and is_anomeric(mol, n1)) and (n0 in u2['extended_ring'] and n0 not in u1['ring_atoms']):
                            donor_idx, acceptor_idx = n1, n0
                            
                        if donor_idx is not None:
                            # 对于头对头的非还原键(如 Sucrose)，防止双向重复记录，仅保留单向
                            if is_anomeric(mol, acceptor_idx) and u1['id'] > u2['id']:
                                continue
                                
                            pos_donor = u1.get('position_map', {}).get(donor_idx, "?")
                            pos_acc_num = u2.get('position_map', {}).get(acceptor_idx, 6 if acceptor_idx not in u2['ring_atoms'] else "?")
                            
                            # 动态生成受体后缀（非还原糖带上 b2/a2 等）
                            if is_anomeric(mol, acceptor_idx):
                                acc_config = u2.get('anomeric_config', '?')
                                pos_acceptor = f"{acc_config}{pos_acc_num}"
                            else:
                                pos_acceptor = str(pos_acc_num)
                                
                            linkages.append({
                                'sugar_donor': u1['id'],
                                'sugar_acceptor': u2['id'],
                                'linkage': f"{pos_donor}-{pos_acceptor}"
                            })
    return linkages

# ---------- 5. Main function from SMILES ----------

def get_sugar_linkages_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")

    sugar_units, atom_to_sugar = get_sugar_units(mol)
    linkages = find_glycosidic_linkages(mol, sugar_units)
    return sugar_units, linkages



# ---------- 6. High-level Validation ----------

def validate_sugar_structure(smiles: str):
    """
    Check if a molecule contains at least one sugar unit.
    Returns: (is_sugar: bool, reason: str)
    """
    if not smiles or not isinstance(smiles, str):
        return False, "Invalid SMILES"

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Parse Error"
        
        sugar_units, _ = get_sugar_units(mol)
        
        if len(sugar_units) > 0:
            return True, f"Found {len(sugar_units)} sugar units"
        else:
            return False, "No sugar rings found"
            
    except Exception as e:
        return False, f"Error: {str(e)}"

if __name__ == "__main__":
    smiles = "CO[C@@H]1O[C@H](CO[C@@H]2O[C@H](COC(=O)OCC3c4ccccc4-c4ccccc43)[C@@H](OCc3ccccc3)[C@H](OCc3ccccc3)[C@H]2OC(=O)c2ccccc2)[C@@H](OCc2ccccc2)[C@H](OCc2ccccc2)[C@H]1OC(=O)c1ccccc1"

    sugars, linkages = get_sugar_linkages_from_smiles(smiles)

    print("Sugars found:")
    for su in sugars:
        print(su['id'], su['position_map'])

    print("\nGlycosidic linkages:")
    for l in linkages:
        print(l)
    
    valid, reason = validate_sugar_structure(smiles)
    print(f"\nValidation: {valid}, {reason}")

# ---------- 7. Advanced Classification & Splitting ----------

# Common Substituents (SMARTS)
SUBSTITUENTS_SMARTS = {
    # Simple Acyl
    "Acetyl": "CC(=O)O",  # Matches the ester oxygen? No, usually attach to O of sugar.
                          # Better to match the group itself attached.
                          # For simplicity, we match the group structure.
    "Malonyl": "OC(=O)CC(=O)O",
    "Succinyl": "OC(=O)CCC(=O)O",
    
    # Phenylpropanoids
    "Caffeoyl": "OC(=O)/C=C/c1ccc(O)c(O)c1",
    "Feruloyl": "OC(=O)/C=C/c1ccc(O)c(OC)c1",
    "Coumaroyl": "OC(=O)/C=C/c1ccc(O)cc1",
    "Sinapoyl": "OC(=O)/C=C/c1cc(OC)c(O)c(OC)c1",
    
    # Others
    "Galloyl": "OC(=O)c1cc(O)c(O)c(O)c1",
    "Benzoyl": "OC(=O)c1ccccc1"
} 

def classify_sugar_parts(mol):
    """
    Classifies atoms into:
    - 'sugar_ring_atoms': Atoms in the pyranose/furanose rings
    - 'sugar_substituent_atoms': Atoms in small substituents attached to sugars (including Acetyl, etc.)
    - 'aglycone_atoms': Atoms in the Aglycone (non-sugar part)
    
    Returns a dict with 3 sets of atom indices.
    """
    sugar_units, atom_to_sugar = get_sugar_units(mol)
    
    all_sugar_unit_atoms = set()
    sugar_ring_atoms = set()
    
    for unit in sugar_units:
        all_sugar_unit_atoms.update(unit['position_map'].keys())
        sugar_ring_atoms.update(unit['ring_atoms'])
        
    sugar_substituent_atoms = set()
    aglycone_atoms = set()
    
    # Pre-calculate partial matches for known substituents to assist classification?
    # Actually, the traversal heuristic (size) + SMARTS is safer.
    
    # Set of atoms already accounted for (Ring + Exocyclic Sugar Carbons)
    processed_atoms = set(all_sugar_unit_atoms) 
    
    # Helper: Traverse a branch
    def get_branch(start_node, origin_node):
        branch_nodes = set()
        queue = [start_node]
        branch_nodes.add(start_node)
        is_linkage = False # If it connects to another sugar
        
        while queue:
            curr = queue.pop(0)
            atom = mol.GetAtomWithIdx(curr)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx == origin_node:
                    continue
                
                # If we hit another sugar unit, it's a linkage/bridge
                if nbr_idx in all_sugar_unit_atoms:
                    is_linkage = True
                    continue # Do not traverse INTO the other sugar
                
                if nbr_idx not in branch_nodes:
                    branch_nodes.add(nbr_idx)
                    queue.append(nbr_idx)
        return branch_nodes, is_linkage

    # Identify "Subs vs Aglycone" by traversing out from Sugar Atoms
    # Sort for determinism
    sorted_sugar_atoms = sorted(list(all_sugar_unit_atoms))
    
    # We need to track which atoms are assigned to avoid double counting
    assigned_atoms = set()
    
    for s_idx in sorted_sugar_atoms:
        atom = mol.GetAtomWithIdx(s_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            
            # Internal sugar bond
            if nbr_idx in all_sugar_unit_atoms:
                continue
            
            # Already assigned (e.g. from another traversal)
            if nbr_idx in sugar_substituent_atoms or nbr_idx in aglycone_atoms:
                continue
                
            # Traverse
            branch, is_linkage = get_branch(nbr_idx, s_idx)
            
            # CLASSIFICATION LOGIC
            branch_size = len(branch)
            # 阈值: 10 原子 — 与 molecular_visualizer.py 的 minAglyconHeavyAtoms=10 保持一致
            # Threshold: 10 atoms — unified with minAglyconHeavyAtoms=10 in molecular_visualizer.py
            # E.g. Acetyl (C2H3O) ~ 4 heavy atoms. 
            # Caffeoyl (C9H7O3) ~ 12 heavy atoms → now classified as Aglycone.
            # Long Carbon chains (>10 atoms) → definitely Aglycone.
            MAX_SUBSTITUENT_SIZE = 10
            
            if is_linkage:
                # If it connects two sugars, it acts as a bridge.
                # If the bridge is huge, it might be aglycone-like, but usually it's just O or O-CH2-O.
                # For visualization, we can color it yellow (substituent-like).
                if branch_size > MAX_SUBSTITUENT_SIZE:
                    # Rare case: A large linker? Treat as aglycone
                    aglycone_atoms.update(branch)
                else:
                    sugar_substituent_atoms.update(branch)
            else:
                # Terminal branch
                if branch_size <= MAX_SUBSTITUENT_SIZE:
                    sugar_substituent_atoms.update(branch)
                else:
                    aglycone_atoms.update(branch)
    
    # Anything not in sugar or substitutents is aglycone (e.g. detached salts or unconnected parts, though rare in single mol)
    # But for a connected graph, the traversal above covers all non-sugar atoms attached to sugar.
    # What if the Aglycone is the "root" and sugars are attached?
    # The loop iterates ALL sugar atoms. So any attachment to ANY sugar is checked.
    # If there are atoms NOT attached to any sugar (e.g. Aglycone atoms deep in the core),
    # they won't be in 'branch' if we only traverse from sugar?
    # Wait, 'get_branch' traverses the *entire* connected component until it hits another sugar or ends.
    # So if we start from a sugar O attached to Aglycone C, we traverse the WHOLE Aglycone.
    
    return {
        'sugar_ring_atoms': sugar_ring_atoms,
        'sugar_substituent_atoms': sugar_substituent_atoms,
        'aglycone_atoms': aglycone_atoms,
        'all_sugar_framework_atoms': all_sugar_unit_atoms # Ring + C6 etc.
    }

def get_split_smiles(mol):
    """
    Splits the molecule into Aglycone and Glycan parts.
    Returns: (aglycone_smiles, glycan_smiles)
    
    Strategy:
    - Identify bonds bridging aglycone and glycan indices.
    - Cleave them using Chem.FragmentOnBonds to generate explicit '*' markers on both sides.
    - Extract disjoint fragments into corresponding SMILES.
    """
    if not mol:
        return "", ""
        
    classification = classify_sugar_parts(mol)
    
    aglycone_indices = classification['aglycone_atoms']
    sugar_part_indices = classification['sugar_ring_atoms'].union(classification['sugar_substituent_atoms']).union(classification['all_sugar_framework_atoms'])
    
    # Identify bonds bridging aglycone and sugar parts
    bonds_to_break = []
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if (a1 in aglycone_indices and a2 in sugar_part_indices) or \
           (a2 in aglycone_indices and a1 in sugar_part_indices):
            bonds_to_break.append(bond.GetIdx())
            
    if bonds_to_break:
        try:
            frag_split = Chem.FragmentOnBonds(mol, bonds_to_break, dummyLabels=[(0, 0)] * len(bonds_to_break))
            frags_indices = list(Chem.GetMolFrags(frag_split))
            
            # Combine all fragments that intersect with aglycone indices into aglycone_smiles
            ag_combined = []
            sug_combined = []
            for indices in frags_indices:
                if any(idx in aglycone_indices for idx in indices):
                    ag_combined.extend(indices)
                elif any(idx in sugar_part_indices for idx in indices):
                    sug_combined.extend(indices)
                    
            aglycone_smiles = Chem.MolFragmentToSmiles(frag_split, atomsToUse=ag_combined, isomericSmiles=True) if ag_combined else ""
            glycan_smiles = Chem.MolFragmentToSmiles(frag_split, atomsToUse=sug_combined, isomericSmiles=True) if sug_combined else ""
        except Exception as e:
            aglycone_smiles = f"Error: {e}"
            glycan_smiles = f"Error: {e}"
    else:
        # 1. Get Aglycone SMILES
        if not aglycone_indices:
            aglycone_smiles = ""
        else:
            try:
                aglycone_smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=list(aglycone_indices), canonical=True, isomericSmiles=True)
            except Exception as e:
                aglycone_smiles = f"Error: {e}"

        # 2. Get Glycan SMILES
        if not sugar_part_indices:
            glycan_smiles = ""
        else:
            try:
                glycan_smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=list(sugar_part_indices), canonical=True, isomericSmiles=True)
            except Exception as e:
                glycan_smiles = f"Error: {e}"
            
    return aglycone_smiles, glycan_smiles

def strict_classify_sugar_parts(mol):
    """
    Enforces Phase 2 Strict Topological Rules.
    1. Identifies Polycyclic Oxygenated structures (False Positives) and rejects them.
    2. Uses classify_sugar_parts if purely valid.
    Returns: (classification_dict, is_false_positive)
    """
    if not mol:
        return None, True
            
    classification = classify_sugar_parts(mol)
    
    # 假阳性剔除：若在这个严格规则下未识别到任何糖环，
    # 说明原数据库标注错误，直接将该化合物从当前处理流中丢弃。
    if not classification['sugar_ring_atoms']:
        return None, True
        
    return classification, False


def strict_cleavage(mol):
    """
    Phase 2 Core Pipeline cleavage.
    Breaks bonds between Aglycan and Glycan Parts and uses Dummy Atom [14*], [15*] to protect topology.
    Returns: (aglycan_smiles_str, glycan_smiles_str, is_false_positive_bool)
    """
    if not mol:
        return "", "", True
        
    try:
        matched_units = find_mapped_sugar_units(mol)
        
        from lib.cleavage_engine import cleave_glycan_aglycan
        glycan_smiles, aglycan_smiles = cleave_glycan_aglycan(mol, matched_units)
        
        if glycan_smiles is None and aglycan_smiles is None:
            return "", "", True
            
        return aglycan_smiles, glycan_smiles, False
        
    except Exception as e:
        print(f"Strict cleavage error: {e}")
        return "", "", True

def extract_amino_acids_and_peptides(mol):
    """
    Extracts 20 standard amino acids and generic amino acids.
    Returns a list of dictionaries for each isolated peptide/amino acid cluster:
    [
      {
        'type': '多肽' or '氨基酸修饰',
        'sequence': 'Ala-Gly-Ser', # or 'Xaa'
        'smiles': '...',
        'atoms': set(...)
      }, ...
    ]
    """
    if not mol: return []
    
    # 按照天然的 alpha-氨基酸骨架:
    # 必须有 N(氨基或酰胺) - C(alpha) - C(=O)(羧基或酰胺) 结构。
    # 核心骨架定义： N - C(alpha) - C(=O) - [O or N]
    # N可以是游离氨基(-NH2, -NH3+)，或者是多肽里的酰胺(-NH-)。
    # C=O必须连接到一个氧(羧基)或氮(下一个氨基酸的酰胺)。
    # 最重要的：Alpha碳绝对禁止存在于普通的3~8元环中！（彻底排除各种生物碱、香豆素、大环内酰胺等假阳性）。
    # 唯有脯氨酸(Pro)的Alpha碳必须在5元环中。
    
    # 1. 查找所有的氨基酸实体
    found_aa_nodes = []
    matched_atoms = set()
    
    # 精确匹配 20 种天然氨基酸 (使用 strict SMARTS)
    for name, pat in AMINO_ACID_LIBRARY.items():
        if pat:
            matches = mol.GetSubstructMatches(pat)
            for m in matches:
                # 检查是否已经被大片段匹配过
                if not any(idx in matched_atoms for idx in m):
                    found_aa_nodes.append({'name': name, 'atoms': set(m)})
                    matched_atoms.update(m)
                    
    # We remove GENERIC_PAT to enforce the strict Whitelist rule requested by the user.
    # If it's not one of the 20 standard AAs matching the strict amide pattern, it is Aglycan.
                
    if not found_aa_nodes:
        return []
        
    # 2. 将氨基酸通过酰胺键(-CO-NH-)聚类拼装成肽链
    clusters = []
    visited_nodes = set()
    
    def are_connected_by_amide(aa1, aa2):
        # We enforce strict peptide bond tracking [NX3][CX3](=O)
        # However, because our dictionary already requires the C=O to be present
        # we just need to verify that an N from one AA connects to the carboxyl C of the other AA.
        for idx1 in aa1['atoms']:
            atom1 = mol.GetAtomWithIdx(idx1)
            for nbr in atom1.GetNeighbors():
                idx2 = nbr.GetIdx()
                if idx2 in aa2['atoms']:
                    # 只允许 C=O 与 N 连接的酰胺键 (或是 C-N 直接相连)
                    if (atom1.GetSymbol() == 'C' and nbr.GetSymbol() == 'N') or \
                       (atom1.GetSymbol() == 'N' and nbr.GetSymbol() == 'C'):
                        # Optional: Could further verify STRICT_PEPTIDE_BOND_SMARTS overlaps these two atoms
                        return True
        return False

    for i in range(len(found_aa_nodes)):
        if i in visited_nodes: continue
        cluster = [found_aa_nodes[i]]
        visited_nodes.add(i)
        
        # 拓展寻找相连的氨基酸 (BFS)
        queue = [i]
        while queue:
            curr_idx = queue.pop(0)
            curr_node = found_aa_nodes[curr_idx]
            for j in range(len(found_aa_nodes)):
                if j not in visited_nodes:
                    if are_connected_by_amide(curr_node, found_aa_nodes[j]):
                        cluster.append(found_aa_nodes[j])
                        visited_nodes.add(j)
                        queue.append(j)
        clusters.append(cluster)
        
    # 3. 提取结果
    results = []
    for c in clusters:
        seq_names = [n['name'] for n in c]
        # 合并所有原子，并通过 BFS 获取它完整的侧链碎片
        core_atoms = set()
        for idx in range(mol.GetNumAtoms()):
            # Fallback for full extraction if needed? 
            # SubstructMatches for the full SMARTS (like Trp, Tyr) already capture the entire side chain.
            pass
            
        full_atoms = set()
        for n in c:
            full_atoms.update(n['atoms'])
            
        is_peptide = len(c) >= 2
        smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=list(full_atoms))
        
        results.append({
            'type': '多肽' if is_peptide else '氨基酸修饰',
            'sequence': "-".join(seq_names),
            'smiles': smiles,
            'atoms': full_atoms
        })
        
    return results
