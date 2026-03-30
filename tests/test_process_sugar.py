import sys
import os
import pytest
from rdkit import Chem

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../scripts")))

from lib import glycan_topology
import process_sugar

def test_identify_extended_sugars():
    # Quinovose (6-deoxy-Glc)
    qui_smiles = "C[C@@H]1[C@H]([C@@H]([C@H]([C@@H](O1)O)O)O)O"
    mol = Chem.MolFromSmiles(qui_smiles)
    sugar_units, _ = sugar_utils.get_sugar_units(mol)
    assert process_sugar.identify_monosaccharide_v2(mol, sugar_units[0]['ring_atoms']) == "Qui"
    
    # Talose
    tal_smiles = "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@H]1O"
    mol = Chem.MolFromSmiles(tal_smiles)
    sugar_units, _ = sugar_utils.get_sugar_units(mol)
    assert process_sugar.identify_monosaccharide_v2(mol, sugar_units[0]['ring_atoms']) == "Tal"

def test_modifications():
    # Acetylated Glucose (at C6)
    # Glc base: OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O
    # At C6: CC(=O)OC...
    # SMILES: CC(=O)OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O
    glc_oac_smiles = "CC(=O)OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O"
    mol = Chem.MolFromSmiles(glc_oac_smiles)
    sugar_units, _ = sugar_utils.get_sugar_units(mol)
    
    # Should contain OAc
    name = process_sugar.identify_monosaccharide_v2(mol, sugar_units[0]['ring_atoms'])
    assert "Glc" in name
    assert "OAc" in name
    
    # GlcNAc (Amino N-Acetyl)
    # Should identify as GlcNAc or GlcN(NAc) -> GlcNAc
    glcnac_smiles = "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O"
    mol = Chem.MolFromSmiles(glcnac_smiles)
    sugar_units, _ = sugar_utils.get_sugar_units(mol)
    name = process_sugar.identify_monosaccharide_v2(mol, sugar_units[0]['ring_atoms'])
    # Logic might return GlcNAc directly from library match
    assert "GlcNAc" in name

def test_sulfate():
    # Sulfated Galactose (at C3?)
    # Gal: OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O
    # C3 is [C@@H](O). Replace O with OS(=O)(=O)O
    # SMILES: OC[C@H]1O[C@@H](O)[C@@H](OS(=O)(=O)O)[C@@H](O)[C@@H]1O
    gal_sulf_smiles = "OC[C@H]1O[C@@H](O)[C@@H](OS(=O)(=O)O)[C@@H](O)[C@@H]1O"
    mol = Chem.MolFromSmiles(gal_sulf_smiles)
    sugar_units, _ = sugar_utils.get_sugar_units(mol)
    
    name = process_sugar.identify_monosaccharide_v2(mol, sugar_units[0]['ring_atoms'])
    assert "Gal" in name
    assert "S" in name

if __name__ == "__main__":
    try:
        print("Testing Extended Sugars...")
        test_identify_extended_sugars()
        print("PASS")
    except Exception as e:
        print(f"FAIL: {e}")
        import traceback
        traceback.print_exc()
        
    try:
        print("Testing Modifications...")
        test_modifications()
        print("PASS")
    except Exception as e:
        print(f"FAIL: {e}")
        traceback.print_exc()
        
    try:
        print("Testing Sulfate...")
        test_sulfate()
        print("PASS")
    except Exception as e:
        print(f"FAIL: {e}")
        traceback.print_exc()
