"""
Diagnostic: Print the RS_LIBRARY contents and compare against what 
get_rs_signature_core produces when parsing maltose SMILES.
"""
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__))))
from rdkit import Chem
from lib.monosaccharide_identifier import get_rs_signature_core, SUGAR_SMILES_LIB

print("\n=== Testing on Maltose ===")
maltose = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O[C@H]1OC(CO)[C@@H](O)[C@H](O)[C@H]1O"
mol = Chem.MolFromSmiles(maltose)
Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
print(f"  Chiral centers in maltose: {centers}")

def debug_sig(mol, ring_atoms):
    conf_map = {idx: conf for idx, conf in centers}
    
    o_idx = None
    for idx in ring_atoms:
        if mol.GetAtomWithIdx(idx).GetSymbol() == "O": o_idx=idx; break
        
    o_atom = mol.GetAtomWithIdx(o_idx)
    neighbors = [n for n in o_atom.GetNeighbors() if n.GetIdx() in ring_atoms]
    n1, n2 = neighbors
    
    def score_c1(atom_idx):
        atom = mol.GetAtomWithIdx(atom_idx)
        score = 0
        for n in atom.GetNeighbors():
            if n.GetIdx() == o_idx: continue
            sym = n.GetSymbol()
            if sym in ["O", "N", "S", "P"]: score += 2
            if sym == "C": score += 0.5 
        return score
        
    s1 = score_c1(n1.GetIdx())
    s2 = score_c1(n2.GetIdx())
    print(f"    Nodes: n1={n1.GetIdx()}(score={s1}), n2={n2.GetIdx()}(score={s2})")
    
    c1 = None
    if s1 > s2: c1 = n1
    elif s2 > s1: c1 = n2
    else: c1 = n1

    print(f"    Selected C1: {c1.GetIdx()}")
    path = [c1.GetIdx()]
    curr = c1
    prev = o_idx
    for _ in range(len(ring_atoms)-2):
        found = False
        for n in curr.GetNeighbors():
            if n.GetIdx() in ring_atoms and n.GetIdx() != prev:
                prev = curr.GetIdx()
                curr = n
                path.append(curr.GetIdx())
                found = True
                break
        if not found: break
        
    print(f"    Path generated: {path}")
    sig = []
    for idx in path[1:]:
         sig.append(conf_map.get(idx, '?'))
    print(f"    Signature: {tuple(sig)}")

ri = mol.GetRingInfo()
for i, ring in enumerate(ri.AtomRings()):
    atoms = [mol.GetAtomWithIdx(a) for a in ring]
    has_o = any(a.GetSymbol() == 'O' for a in atoms)
    if not has_o: continue
    print(f"\n  Ring {i}: atoms={list(ring)}")
    debug_sig(mol, list(ring))

