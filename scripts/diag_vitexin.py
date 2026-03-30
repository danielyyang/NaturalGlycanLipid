"""Quick Vitexin CIP diagnosis"""
import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import json
from rdkit import Chem

with open("data/benchmark_200.json", "r", encoding="utf-8") as f:
    entries = json.load(f)

for entry in entries:
    if entry["id"] == 10:
        smi = entry["smiles"]
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print("SMILES PARSE FAILED for Vitexin!")
            print("SMILES: %s" % smi)
            break
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
        print("Chiral centers: %s" % str(centers))

        # SubstructMatch with D-Glc
        glcSmi = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"
        glcMol = Chem.MolFromSmiles(glcSmi)
        m1 = mol.HasSubstructMatch(glcMol, useChirality=True)
        m2 = mol.HasSubstructMatch(glcMol, useChirality=False)
        print("SubstructMatch D-Glc (chiral): %s" % m1)
        print("SubstructMatch D-Glc (no chiral): %s" % m2)

        # Check CIP on sugar ring atoms
        ri = mol.GetRingInfo()
        for ring in ri.AtomRings():
            if len(ring) == 6:
                oCount = sum(1 for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 8)
                if oCount == 1:
                    print("Sugar-like ring:")
                    for i in ring:
                        a = mol.GetAtomWithIdx(i)
                        cip = a.GetProp("_CIPCode") if a.HasProp("_CIPCode") else "?"
                        print("  Atom %d (%s): CIP=%s" % (i, a.GetSymbol(), cip))
        break
