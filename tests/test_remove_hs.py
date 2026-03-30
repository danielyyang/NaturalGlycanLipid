import sys, os
from rdkit import Chem
from lib.monosaccharide_identifier import get_rs_signature_core

maltose = "O[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)[C@@H](CO)O1"
mol = Chem.MolFromSmiles(maltose)
Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
ri = mol.GetRingInfo()

ring0 = ri.AtomRings()[0]
from lib.monosaccharide_identifier import isolate_sugar_ring
clean, mapping = isolate_sugar_ring(mol, list(ring0))
print("Before RemoveHs:")
print(Chem.MolToSmiles(clean, isomericSmiles=True))
print(f"Centers: {Chem.FindMolChiralCenters(clean)}")
mapped_ring = [mapping[idx] for idx in ring0]
print(f"Signature: {get_rs_signature_core(clean, mapped_ring)}")

print("\nAfter RemoveHs:")
clean_no_h = Chem.RemoveHs(clean)
Chem.AssignStereochemistry(clean_no_h, force=True, cleanIt=True)
print(Chem.MolToSmiles(clean_no_h, isomericSmiles=True))
print(f"Centers: {Chem.FindMolChiralCenters(clean_no_h)}")

# Wait, we must remap the ring atoms because RemoveHs changes atom indices!
match = clean.GetSubstructMatch(clean_no_h)
# However, RemoveHs simply deletes atoms. The remaining atoms usually keep their old indices relative to each other? No, they shift.
# The safest way is to just find the ring in clean_no_h.
new_ri = clean_no_h.GetRingInfo()
new_ring = new_ri.AtomRings()[0]
print(f"Signature: {get_rs_signature_core(clean_no_h, list(new_ring))}")

