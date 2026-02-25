from lib.chemical_dictionaries import (
    MONOSACCHARIDE_LIBRARY, 
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

def find_mapped_sugar_units(mol):
    """
    1. Finds all potential sugar rings topologically.
    2. Identifies C1 -> C_n and maps positions.
    3. Uses lib.sugar_sequence to identify basic sugar names.
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
    for ring in potential_rings:
        is_fused = any(ri.NumAtomRings(idx) > 1 for idx in ring)
        if is_fused: continue
            
        invalid_connection = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for bond in atom.GetBonds():
                if bond.GetOtherAtomIdx(idx) not in ring and bond.GetBondType() != Chem.BondType.SINGLE:
                    invalid_connection = True
                    break
            if invalid_connection: break
        if invalid_connection: continue
            
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
        if is_nucleoside: continue

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
            
        valid_rings.append(ring)
        
    import lib.sugar_sequence as seq
        
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
        identify_result = seq.identify_monosaccharide_v2(mol, ring)
        if isinstance(identify_result, tuple):
            base_name, anomer_config = identify_result
        else:
            base_name = identify_result
            anomer_config = "?"
            
        if "(" in base_name:
            name = base_name.split("(")[0]
            seq_mods = [m.strip() for m in base_name.split("(")[1].replace(")", "").split(",")]
        else:
            name = base_name
            seq_mods = []
            
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
            # Threshold: 15 atoms. 
            # E.g. Acetyl (C2H3O) ~ 4 heavy atoms. 
            # Caffeoyl (C9H7O3) ~ 12 heavy atoms.
            # Aglycone usually Steroid (C17+) or Triterpene (C30+).
            MAX_SUBSTITUENT_SIZE = 15 
            
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
    Breaks bonds between Aglycan and Glycan Parts and uses Dummy Atom (*) to protect topology.
    Returns: (aglycan_smiles_str, glycan_smiles_str, is_false_positive_bool)
    """
    class_dict, is_fp = strict_classify_sugar_parts(mol)
    if is_fp:
        return "", "", True
        
    aglycone_indices = class_dict['aglycone_atoms']
    sugar_part_indices = class_dict['sugar_ring_atoms'].union(class_dict['sugar_substituent_atoms']).union(class_dict['all_sugar_framework_atoms'])
    
    bonds_to_break = []
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if (a1 in aglycone_indices and a2 in sugar_part_indices) or \
           (a2 in aglycone_indices and a1 in sugar_part_indices):
            bonds_to_break.append(bond.GetIdx())
            
    if not bonds_to_break:
        if not aglycone_indices: return "", Chem.MolToSmiles(mol), False
        if not sugar_part_indices: return Chem.MolToSmiles(mol), "", False
        
    try:
        # Generate dummy markers (*) at broken connections
        frag_split = Chem.FragmentOnBonds(mol, bonds_to_break, addDummies=True) # default adds *
        frags_indices = list(Chem.GetMolFrags(frag_split))
        
        ag_combined = []
        sug_combined = []
        
        # We need to map the original indices to the new split fragments.
        # FragmentOnBonds adds dummy atoms which will have new indices > mol.GetNumAtoms().
        # So we check if ANY of the original indices match our sub-parts to classify the fragment.
        
        for indices in frags_indices:
            # Check for intersection with original indices
            if any(idx in aglycone_indices for idx in indices if idx < mol.GetNumAtoms()):
                ag_combined.extend(indices)
            elif any(idx in sugar_part_indices for idx in indices if idx < mol.GetNumAtoms()):
                sug_combined.extend(indices)
                
        aglycone_smiles = Chem.MolFragmentToSmiles(frag_split, atomsToUse=ag_combined, isomericSmiles=True) if ag_combined else ""
        glycan_smiles = Chem.MolFragmentToSmiles(frag_split, atomsToUse=sug_combined, isomericSmiles=True) if sug_combined else ""
        
        return aglycone_smiles, glycan_smiles, False
        
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
