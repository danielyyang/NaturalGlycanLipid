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
    import chemical_dictionaries
    from chemical_dictionaries import MONOSACCHARIDE_LIBRARY
except ImportError:
    from lib import chemical_dictionaries
    from lib.chemical_dictionaries import MONOSACCHARIDE_LIBRARY

try:
    import sugar_utils
except ImportError:
    from lib import sugar_utils

# ---------- 1. Reference Library Definitions ----------

# [THE COMPREHENSIVE NATURAL PRODUCTS DICTIONARY] D & L Enantiomers
RAW_MONOSACCHARIDE_SMILES = {
    # === Hexopyranoses (6-membered Hexoses) ===
    # D-Glucose
    ("D-Glc", "a"): "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    ("D-Glc", "b"): "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    # L-Glucose (Rare, but included for completeness)
    ("L-Glc", "a"): "O[C@@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](CO)O1",
    ("L-Glc", "b"): "O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](CO)O1",
    
    # D-Galactose
    ("D-Gal", "a"): "O[C@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@@H](CO)O1",
    ("D-Gal", "b"): "O[C@@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@@H](CO)O1",
    # L-Galactose (Found in agar and some plant cell walls)
    ("L-Gal", "a"): "O[C@@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H](CO)O1",
    ("L-Gal", "b"): "O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H](CO)O1",
    
    # D-Mannose
    ("D-Man", "a"): "O[C@H]1[C@@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    ("D-Man", "b"): "O[C@@H]1[C@@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    # L-Mannose (Base skeleton for L-Rhamnose which is 6-deoxy-L-Man)
    ("L-Man", "a"): "O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](CO)O1",
    ("L-Man", "b"): "O[C@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](CO)O1",
    
    # === Furanoses (5-membered Ketoses) ===
    ("D-Fru", "a"): "OC[C@]1(O)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    ("D-Fru", "b"): "OC[C@@]1(O)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    
    # === Pentopyranoses (6-membered Pentoses, vital for Natural Products) ===
    # D-Xylose
    ("D-Xyl", "a"): "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)CO1",
    ("D-Xyl", "b"): "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)CO1",
    
    # L-Arabinose (Extremely common in saponins and flavonoids)
    ("L-Ara", "a"): "O[C@H]1[C@@H](O)[C@@H](O)[C@@H](O)CO1",
    ("L-Ara", "b"): "O[C@@H]1[C@@H](O)[C@@H](O)[C@@H](O)CO1",
    # D-Arabinose
    ("D-Ara", "a"): "O[C@@H]1[C@H](O)[C@H](O)[C@H](O)CO1",
    ("D-Ara", "b"): "O[C@H]1[C@H](O)[C@H](O)[C@H](O)CO1",
}
def create_query_mol(smiles):
    """将标准 SMILES 转化为无视隐式氢的通用 Query 分子"""
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return None
    qmol = Chem.RWMol(mol)
    for atom in qmol.GetAtoms():
        if atom.GetAtomicNum() == 8:
            # 将所有氧原子替换为通用氧查询，允许其在聚合物中连接其他重原子
            q_atom = Chem.rdqueries.AtomNumEqualsQueryAtom(8)
            qmol.ReplaceAtom(atom.GetIdx(), q_atom)
    return qmol.GetMol()

# 使用改造后的通用查询图预编译字典
REFERENCE_MOLS = {k: create_query_mol(v) for k, v in RAW_MONOSACCHARIDE_SMILES.items() if create_query_mol(v) is not None}

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
                            if 'C(=O)C' in sub_smiles or sub_smiles == 'CC(=O)': mods.append('Ac')
                            elif 'S(=O)(=O)' in sub_smiles: mods.append('S')
                            elif sub_smiles == 'C': mods.append('Me')
                            elif 'P(=O)' in sub_smiles: mods.append('P')
                            elif len(branch) >= 6: is_complex = True
                    except:
                        pass
                        
    except Exception as e:
        pass
        
    final_mods = sorted(list(set(mods)))
    return final_mods, is_acid, is_complex

def identify_monosaccharide_v2(mol, ring_atoms):
    base_name = "Hex" if len(ring_atoms) == 6 else "Pen"
    matched_anomer = "?"
    
    # 直接在完整的 polymer (mol) 上进行子图透视匹配，绝不切割分子！
    for (name, anomer), ref_mol in REFERENCE_MOLS.items():
        matches = mol.GetSubstructMatches(ref_mol, useChirality=True)
        for match in matches:
            # 如果这个匹配的子图完美覆盖了我们当前的 ring_atoms，说明就是它！
            if all(idx in match for idx in ring_atoms):
                base_name = name
                matched_anomer = anomer
                break
        if matched_anomer != "?":
            break
                
    mods, is_acid, is_complex = check_modifications(mol, ring_atoms)
    
    if "GlcNAc" in base_name or "GalNAc" in base_name:
        if "NAc" in mods: mods.remove("NAc")
    if "GlcA" in base_name or "GalA" in base_name or "IdoA" in base_name:
        if "A" in mods: mods.remove("A")
        
    if base_name == "Glc" and "A" in mods:
        base_name = "GlcA"; mods.remove("A")
    elif base_name == "Gal" and "A" in mods:
         base_name = "GalA"; mods.remove("A")
         
    final_name = f"{base_name}({','.join(mods)})" if mods else base_name
    return final_name, matched_anomer

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

def generate_refined_sequence(mol):
    if not mol: return "", ""
    # Inject identify back to sugar_utils to ensure compatibility
    sugar_units, atom_to_sugar = sugar_utils.get_sugar_units(mol)
    if not sugar_units: return "", ""
    
    G = nx.DiGraph()
    
    for i, unit in enumerate(sugar_units):
        unit['order_id'] = i + 1
        name = unit.get('name', 'Unknown')
        anomer = unit.get('anomeric_config', '?')
        mods = unit.get('modifications', [])
        
        G.add_node(unit['id'], name=name, anomer=anomer, mods=mods, order_id=unit['order_id'])
        
    linkages = sugar_utils.find_glycosidic_linkages(mol, sugar_units)
    
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
            mod_str = ",".join(nd['mods']) if nd['mods'] else ""
            mod_parts.append(f"{nd['name']}_{nd['order_id']}({mod_str})")
        return "Cyclic", " ; ".join(mod_parts)
    
    def traverse_two_pass(node):
        nd = G.nodes[node]
        name = nd['name']
        mods = nd['mods']
        order_id = nd['order_id']
        
        mod_label = f"{name}_{order_id}({','.join(mods)})" if mods else f"{name}_{order_id}()"
        
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
