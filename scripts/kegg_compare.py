"""Compare current dictionary SMILES vs KEGG authoritative SMILES."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from rdkit import Chem

# KEGG authoritative (fetched from KEGG REST API)
KEGG_SMILES = {
    "D-Glc": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    "D-Gal": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O",
    "L-Gal": "OC[C@@H]1OC(O)[C@@H](O)[C@H](O)[C@@H]1O",
    "D-Man": "OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O",
    "D-All": "OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@@H]1O",
    "D-Tal": "OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@H]1O",
    "D-Gul": "OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@H]1O",
    "L-Alt": "OC[C@@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O",
    "L-Ido": "OC[C@@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    "L-Rha": "C[C@@H]1OC(O)[C@H](O)[C@H](O)[C@H]1O",
    "D-Rha": "C[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",
    "L-Fuc": "C[C@@H]1OC(O)[C@@H](O)[C@H](O)[C@@H]1O",
    "D-Fuc": "C[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O",
    "D-Qui": "C[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
    "D-Xyl": "OC1OC[C@@H](O)[C@H](O)[C@H]1O",
    "D-Rib": "OC[C@H]1OC(O)[C@H](O)[C@@H]1O",
    "D-Lyx": "OC1OC[C@@H](O)[C@H](O)[C@@H]1O",
}

# Load our current dictionary
from lib.glycan_reference_library import RAW_MONOSACCHARIDE_SMILES

print(f"{'Sugar':<12} {'Match?':<8} {'Ours (canonical)':<50} {'KEGG (canonical)':<50}")
print("=" * 120)

mismatches = []
for name, keggSmi in sorted(KEGG_SMILES.items()):
    # Find our entry (alpha anomer)
    ourKey = (name, "a")
    ourSmi = RAW_MONOSACCHARIDE_SMILES.get(ourKey, None)
    
    if ourSmi is None:
        print(f"{name:<12} {'MISS':<8} {'NOT IN DICT':<50} {keggSmi}")
        mismatches.append(name)
        continue
    
    # Canonicalize both (RDKit canonical is the gold standard)
    ourMol = Chem.MolFromSmiles(ourSmi)
    keggMol = Chem.MolFromSmiles(keggSmi)
    
    if ourMol is None:
        print(f"{name:<12} {'ERR':<8} {'INVALID SMILES':<50} {keggSmi}")
        mismatches.append(name)
        continue
    if keggMol is None:
        print(f"{name:<12} {'K-ERR':<8} {ourSmi:<50} {'INVALID KEGG SMILES'}")
        continue
    
    ourCanon = Chem.MolToSmiles(ourMol, isomericSmiles=True)
    keggCanon = Chem.MolToSmiles(keggMol, isomericSmiles=True)
    
    # KEGG uses anomeric-unspecified C1, ours uses specific alpha
    # Compare InChI instead for stereochemistry comparison
    ourInchi = Chem.MolToInchi(ourMol)
    keggInchi = Chem.MolToInchi(keggMol)
    
    match = "OK" if ourInchi == keggInchi else "DIFF"
    if match == "DIFF":
        mismatches.append(name)
    
    print(f"{name:<12} {match:<8} {ourCanon:<50} {keggCanon}")
    if match == "DIFF":
        print(f"{'':12} {'':8} InChI ours: {ourInchi}")
        print(f"{'':12} {'':8} InChI kegg: {keggInchi}")

print(f"\n{'='*120}")
print(f"Mismatches: {len(mismatches)} / {len(KEGG_SMILES)}")
if mismatches:
    print(f"  {', '.join(mismatches)}")
