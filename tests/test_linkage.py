import sys, os
from rdkit import Chem

import traceback
sys.path.append(os.path.abspath(os.path.dirname(os.path.dirname(__file__))))
from lib import sugar_utils
from lib import sugar_sequence

# [TEST DATA ONLY]
TEST_CASES = {
    # === 1-4 连接对照组 ===
    "Maltose (Glc-a1,4-Glc)": "O[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)[C@@H](CO)O1",
    "Cellobiose (Glc-b1,4-Glc)": "O[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)[C@@H](CO)O1",
    "Lactose (Gal-b1,4-Glc)": "O[C@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]2[C@H](O)[C@@H](O)[C@@H](O)[C@@H](CO)O2)[C@@H](CO)O1",
    
    # === 1-6 连接对照组 ===
    "Isomaltose (Glc-a1,6-Glc)": "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO[C@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)O1",
    "Gentiobiose (Glc-b1,6-Glc)": "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO[C@@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)O1",
    "Melibiose (Gal-a1,6-Glc)": "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO[C@H]2[C@H](O)[C@@H](O)[C@@H](O)[C@@H](CO)O2)O1",
    
    # === 1-2 非还原连接对照组 ===
    "Sucrose (Glc-a1,b2-Fru)": "OC[C@@]1(O[C@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO)O2)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    
    # === 复杂多糖对照组 ===
    "Stachyose (Gal-a1,6-Gal-a1,6-Glc-a1,b2-Fru)": "OC[C@@]1(O[C@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](CO[C@H]3[C@H](O)[C@@H](O)[C@@H](O)[C@@H](CO[C@H]4[C@H](O)[C@@H](O)[C@@H](O)[C@@H](CO)O4)O3)O2)[C@@H](O)[C@H](O)[C@@H](CO)O1"
}

for name, smi in TEST_CASES.items():
    print(f"\n{'='*60}")
    print(f"Testing: {name}")
    print(f"SMILES: {smi}")
    try:
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            print("  ERROR: Could not parse SMILES")
            continue
            
        units, atom_map = sugar_utils.get_sugar_units(mol)
        print(f"  Sugar units found: {len(units)}")
        for u in units:
            print(f"    - id={u['id']}, name={u['name']}, anomer={u['anomeric_config']}, mods={u['modifications']}")
            
        linkages = sugar_utils.find_glycosidic_linkages(mol, units)
        print(f"  Linkages found: {len(linkages)}")
        for l in linkages:
            print(f"    - Donor sugar {l['sugar_donor']} -> Acceptor sugar {l['sugar_acceptor']}, link: {l['linkage']}")
        
        seq, mods = sugar_sequence.generate_refined_sequence(mol)
        print(f"  Sugar_Sequence:        {seq}")
        print(f"  Sugar_Functional_Group: {mods}")
    except Exception as e:
        print(f"  CRASHED: {e}")
        traceback.print_exc()
