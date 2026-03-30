"""Diagnose last 4 failures in detail"""
import sys, os, json
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from rdkit import Chem
from lib.glycan_topology import find_mapped_sugar_units
from lib.monosaccharide_identifier import identify_monosaccharide_v10

with open("data/benchmark_200.json", "r", encoding="utf-8") as f:
    entries = json.load(f)

targets = {10: "Vitexin", 23: "Thymidine", 24: "Cytarabine", 49: "Amphotericin B"}

for entry in entries:
    if entry["id"] not in targets:
        continue
    idx = entry["id"]
    smi = entry["smiles"]
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print("[%d] PARSE FAIL" % idx)
        continue

    print("=== [%d] %s ===" % (idx, entry["name"]))
    print("  Expected: %s" % entry.get("expected_sugars", []))
    
    units = find_mapped_sugar_units(mol)
    print("  Sugar units found: %d" % len(units))
    for u in units:
        ra = u.get("ring_atoms", [])
        name = u.get("name", "?")
        mods = u.get("modifications", [])
        print("    name=%s ring_size=%d mods=%s" % (name, len(ra), mods))
        
        # Also run v10 directly
        result = identify_monosaccharide_v10(mol, ra)
        if isinstance(result, tuple):
            v10Name, v10Anom = result
        else:
            v10Name = result
            v10Anom = "?"
        print("    v10 result: %s (%s)" % (v10Name, v10Anom))
    print()
