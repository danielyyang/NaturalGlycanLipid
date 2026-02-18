
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from lib import sugar_utils
from rdkit import Chem

# A SMILES from the user's Ginsenoside.csv (Row 2, Notoginsenoside G)
smiles = "CC(C)=CCC[C@](C)(O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H]1CC[C@]2(C)[C@@H]1[C@H](O)C[C@@H]1[C@@]2(C)[C@@H](O)C=C2C(C)(C)[C@@H](O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)CC[C@@]21C"

print(f"Testing SMILES: {smiles}")

try:
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        print("RDKit Mol: Valid")
        print(f"Num Atoms: {mol.GetNumAtoms()}")
        
        sugar_units, _ = sugar_utils.get_sugar_units(mol)
        print(f"Sugar Units Found: {len(sugar_units)}")
        for i, su in enumerate(sugar_units):
            print(f"  Sugar {i+1}: {su}")
            
        valid, reason = sugar_utils.validate_sugar_structure(smiles)
        print(f"Validation Result: {valid}, Reason: {reason}")
        
        # Test Visualization
        from lib.visualizer import StructureVisualizer
        viz = StructureVisualizer()
        output_img = "images/debug_ginsenoside.png"
        viz.analyze_glycolipid(smiles, output_img)
        print(f"Visualization saved to: {output_img}")
        
    else:
        print("RDKit Mol: Invalid")
except Exception as e:
    print(f"Error: {e}")
