from rdkit import Chem

# ---------- 1. Find sugar rings (furanose & pyranose) ----------

def find_sugar_rings(mol):
    """
    Returns a list of sugar-like rings.
    Here: 5- or 6-membered ring with exactly 1 oxygen (furanose or pyranose).
    """
    ri = mol.GetRingInfo()
    sugar_rings = []
    for ring in ri.AtomRings():
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        symbols = [a.GetSymbol() for a in atoms]
        # Accept 5- or 6-member rings with exactly one O
        if len(ring) in (5, 6) and symbols.count("O") == 1:
            sugar_rings.append(list(ring))
    return sugar_rings


# ---------- 2. Number a sugar ring (C1–C5/6 depending on size) ----------

def number_sugar_ring(mol, ring_atoms):
    """
    Given a sugar-like ring (5 or 6 members, 1 O),
    returns:
      - position_map: {atom_idx: position (1..5 or 1..6)}
      - ring_oxygen_idx

    Conventions:
      - The ring is O–C1–C2–... (or the inverse order).
      - C1 is the neighbor of the ring oxygen that has an exocyclic oxygen.
      - For a 6-membered ring: ring carbons are C1–C5, exocyclic CH2 is C6.
      - For a 5-membered ring: ring carbons are C1–C4, exocyclic CH2 is C5.
    """
    ring_atoms = list(ring_atoms)
    ring_size = len(ring_atoms)
    atom_objs = {i: mol.GetAtomWithIdx(i) for i in ring_atoms}

    # sanity check: we only handle 5- or 6-member rings here
    if ring_size not in (5, 6):
        return None, None

    # find the ring oxygen
    ring_oxygen_idx = None
    for idx in ring_atoms:
        if atom_objs[idx].GetSymbol() == "O":
            ring_oxygen_idx = idx
            break
    if ring_oxygen_idx is None:
        return None, None  # not a sugar ring

    # rotate list so it starts at the ring oxygen
    pos_O = ring_atoms.index(ring_oxygen_idx)
    ring_order = ring_atoms[pos_O:] + ring_atoms[:pos_O]
    # now ring_order[0] = ring oxygen;
    # neighbors in the ring: ring_order[1] and ring_order[-1]
    neighbor1 = ring_order[1]
    neighbor2 = ring_order[-1]

    def has_exocyclic_oxygen(atom_idx):
        atom = mol.GetAtomWithIdx(atom_idx)
        for n in atom.GetNeighbors():
            if n.GetIdx() == ring_oxygen_idx:
                continue
            if n.GetIdx() in ring_order:
                continue
            if n.GetSymbol() == "O":
                return True
        return False

    # choose C1: neighbor of O that has an exocyclic oxygen
    c1_candidate = None
    if has_exocyclic_oxygen(neighbor1):
        c1_candidate = neighbor1
    if c1_candidate is None and has_exocyclic_oxygen(neighbor2):
        c1_candidate = neighbor2

    # coarse fallback: if none found, assume neighbor1 is C1
    if c1_candidate is None:
        c1_candidate = neighbor1

    # list of ring carbons in order
    # ring_order = [O, C?, C?, C?, ...]
    carbons_in_ring = ring_order[1:]  # all atoms except the ring O

    # orient the ring so that carbons_path = [C1, C2, ...]
    if c1_candidate == neighbor1:
        carbons_path = carbons_in_ring[:]  # already in the correct direction
    else:
        carbons_path = list(reversed(carbons_in_ring))

    # map ring carbon positions:
    # for 6-member ring: C1..C5 (5 carbons in ring)
    # for 5-member ring: C1..C4 (4 carbons in ring)
    position_map = {}
    for i, idx in enumerate(carbons_path):
        # there will be ring_size-1 carbons:
        #  - if ring_size=6: 5 carbons -> positions 1..5
        #  - if ring_size=5: 4 carbons -> positions 1..4
        position_map[idx] = i + 1

    # Find the "terminal" ring carbon (last in the path) to search for an exocyclic CH2 group
    c_last_idx = carbons_path[-1]
    c_last_atom = mol.GetAtomWithIdx(c_last_idx)
    exocyclic_idx = None
    for n in c_last_atom.GetNeighbors():
        if n.GetIdx() in ring_order:
            continue
        if n.GetSymbol() == "C":
            if n.GetHybridization() == Chem.HybridizationType.SP3:
                exocyclic_idx = n.GetIdx()
                break

    # If we find an exocyclic carbon, assign it the next position:
    #  - for ring_size=6 (pyranose): position 6
    #  - for ring_size=5 (furanose): position 5
    if exocyclic_idx is not None:
        next_pos = len(carbons_path) + 1
        position_map[exocyclic_idx] = next_pos

    return position_map, ring_oxygen_idx


# ---------- 3. Build sugar units ----------

def get_sugar_units(mol):
    """
    Steps:
      - detect sugar-like rings (5- or 6-member)
      - number each one (C1..C5/6)
    Returns a list of dicts:
      {
        'id': int,
        'ring_atoms': [...],
        'ring_oxygen': int,
        'position_map': {atom_idx: position}
      }
    And also a mapping atom_idx -> sugar_id.
    """
    sugar_rings = find_sugar_rings(mol)
    sugar_units = []
    atom_to_sugar = {}

    for sid, ring in enumerate(sugar_rings):
        pos_map, ring_O = number_sugar_ring(mol, ring)
        if pos_map is None:
            continue
        sugar_units.append({
            'id': sid,
            'ring_atoms': list(ring),
            'ring_oxygen': ring_O,
            'position_map': pos_map
        })
        for a_idx in pos_map.keys():
            atom_to_sugar[a_idx] = sid

    return sugar_units, atom_to_sugar


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


def find_glycosidic_linkages(mol, sugar_units, atom_to_sugar):
    """
    Returns a list of glycosidic linkages between sugar units:
      [
        {
          'sugar_donor': id,
          'sugar_acceptor': id,
          'pos_donor': 1,
          'pos_acceptor': 6,
          'linkage': '1-6'
        },
        ...
      ]
    Convention: 'donor' is the sugar where the bond originates from C1.
    """
    linkages = []
    for atom in mol.GetAtoms():
        if not is_bridge_oxygen_between_sugars(mol, atom, atom_to_sugar):
            continue

        # carbon neighbors of the bridge oxygen
        neighbors = [b.GetOtherAtom(atom) for b in atom.GetBonds()]
        c1, c2 = neighbors[0], neighbors[1]
        idx1, idx2 = c1.GetIdx(), c2.GetIdx()

        s1 = atom_to_sugar[idx1]
        s2 = atom_to_sugar[idx2]
        # get position maps for the two sugars
        su1 = next(su for su in sugar_units if su['id'] == s1)
        su2 = next(su for su in sugar_units if su['id'] == s2)

        pos1 = su1['position_map'].get(idx1, None)
        pos2 = su2['position_map'].get(idx2, None)

        if pos1 is None or pos2 is None:
            # could be a bridge at a non-numbered position (e.g. weird substituent), ignore
            continue

        # decide who is "donor" (who has C1)
        if pos1 == 1 and pos2 != 1:
            sugar_donor, sugar_acceptor = s1, s2
            pos_donor, pos_acceptor = pos1, pos2
        elif pos2 == 1 and pos1 != 1:
            sugar_donor, sugar_acceptor = s2, s1
            pos_donor, pos_acceptor = pos2, pos1
        else:
            # if neither side is C1 (or both are 1, e.g. 1–1), just record symmetrically
            sugar_donor, sugar_acceptor = s1, s2
            pos_donor, pos_acceptor = pos1, pos2

        linkage_str = f"{pos_donor}-{pos_acceptor}"

        linkages.append({
            'sugar_donor': sugar_donor,
            'sugar_acceptor': sugar_acceptor,
            #'pos_donor': pos_donor,
            #'pos_acceptor': pos_acceptor,
            'linkage': linkage_str,
            #'bridge_oxygen_idx': atom.GetIdx()
        })

    return linkages


# ---------- 5. Main function from SMILES ----------

def get_sugar_linkages_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES")

    sugar_units, atom_to_sugar = get_sugar_units(mol)
    linkages = find_glycosidic_linkages(mol, sugar_units, atom_to_sugar)
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
    - Aglycone SMILES: Mask Sugar+Substituents, or extract Aglycone atoms?
      - Better: Delete Sugar+Subst atoms to get Aglycone.
    - Glycan SMILES: Delete Aglycone atoms to get Glycan.
    """
    if not mol:
        return "", ""
        
    classification = classify_sugar_parts(mol)
    
    aglycone_indices = classification['aglycone_atoms']
    sugar_part_indices = classification['sugar_ring_atoms'].union(classification['sugar_substituent_atoms']).union(classification['all_sugar_framework_atoms'])
    
    # 1. Get Aglycone SMILES
    # If no aglycone atoms, return empty
    if not aglycone_indices:
        aglycone_smiles = ""
    else:
        # Create a new editable mol
        try:
            # We want to keep Aglycone atoms.
            # Easiest way in RDKit is to determine atoms TO DELETE.
            # Delete all atoms NOT in aglycone_indices?
            # Wait, RDKit atom indices change after deletion. Use specific function provided by generic RDKit flow or EditableMol.
            
            # But wait, we want to show attachment points?
            # User request: "Split SMILES". Usually means the isolated structure.
            # If we just cut bonds, we get radicals or hydrogens.
            # RDKit's PathToSubmol or similar might accept a list of atom indices.
            
            # APPROACH: Use MolFragmentToSmiles
            # It allows specifying atomsToUse.
            aglycone_smiles = Chem.MolFragmentToSmiles(mol, atomsToUse=list(aglycone_indices), canonical=True, isomericSmiles=True)
            
            # Note: This might create disjoint fragments if the aglycone is effectively split by sugars (unlikely for saponins, usually 1 aglycone core).
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
