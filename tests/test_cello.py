import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__))))
from rdkit import Chem
from lib.monosaccharide_identifier import isolate_sugar_ring, get_rs_signature_core, RS_LIBRARY

cellobiose = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1OC1OC(CO)[C@@H](O)[C@H](O)[C@H]1O"
mol = Chem.MolFromSmiles(cellobiose)
ri = mol.GetRingInfo()

for i, ring in enumerate(ri.AtomRings()):
    atoms = [mol.GetAtomWithIdx(a).GetSymbol() for a in ring]
    if atoms.count('O') == 1:
        clean, mapping = isolate_sugar_ring(mol, list(ring))
        if clean and mapping:
            mapped_ring = [mapping[idx] for idx in ring]
            sig = get_rs_signature_core(clean, mapped_ring)
            print(f"Ring {i} cleaned centers: {Chem.FindMolChiralCenters(clean)}")
            print(f"SMILES: {Chem.MolToSmiles(clean, isomericSmiles=True)}")
            print(f"Signature: {sig}")
            if sig in RS_LIBRARY:
                print(f"  Matched: {RS_LIBRARY[sig]}")
            else:
                print(f"  NO MATCH")
