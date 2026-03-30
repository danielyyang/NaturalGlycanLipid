import os
import sys
from rdkit import Chem

# 绝对引入根目录库 (Absolute import root libraries)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from lib.glycan_topology import find_mapped_sugar_units

def test_stapled_macrocycle():
    """
    盲测复杂的宏环糖脂 (Macrocyclic Glycolipid) SMILES
    Blind test for complex macrocyclic glycolipid SMILES
    """
    smiles = "CCCCCCCCCCCCCCCC(=O)OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O[C@@H]3O[C@@H](C)[C@H](OC(=O)CCCCCCCCCCCCCCC)[C@@H](O)[C@H]3O)[C@H](C)O[C@@H](O[C@@H]3[C@H]4OC(=O)CCCCCCCCC[C@H](CCCCC)O[C@@H]5O[C@H](C)[C@H](O)[C@H](O)[C@H]5O[C@H](O[C@H]3C)[C@@H]4O)[C@@H]2OC(=O)CCCCCCCCCCCCCCC)[C@H](O)[C@@H](O)[C@@H]1O"
    print("Testing Stapled Macrocyclic Glycolipid SMILES:")
    print(smiles)
    print("-" * 50)
    
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print("Failed to parse SMILES!")
        return
        
    try:
        mapped_units = find_mapped_sugar_units(mol)
        
        if not mapped_units:
            print("[FAIL] Test Failed: Returned NO sugar units (FALSE_POSITIVE)")
        else:
            print(f"[SUCCESS] Test Passed: Found {len(mapped_units)} sugar units!")
            for idx, unit in enumerate(mapped_units):
                print(f"  Unit {idx + 1}: {unit.get('name', 'Unknown')}")
                print(f"    Ring Atoms: {unit.get('ring_atoms', [])}")
                print(f"    Modifications: {unit.get('modifications', [])}")
                
    except Exception as e:
        print(f"[ERROR] Exception occurred during testing: {e}")

if __name__ == "__main__":
    test_stapled_macrocycle()
