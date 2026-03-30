"""Fix Saikosaponin A SMILES in benchmark_200.json"""
import json
import os

benchPath = os.path.join(os.path.dirname(__file__), "..", "data", "benchmark_200.json")
with open(benchPath, "r", encoding="utf-8") as f:
    entries = json.load(f)

# Saikosaponin A CID 167928 — correct isomeric SMILES from PubChem
CORRECT_SMI = "C[C@@H]1[C@@H]([C@@H]([C@H]([C@@H](O1)O[C@H]2CC[C@]3([C@H]([C@]2(C)CO)CC[C@@]4([C@@H]3C=C[C@@]56[C@]4(C[C@@H]([C@@]7([C@H]5CC(CC7)(C)C)CO6)O)C)C)C)O)O[C@H]8[C@@H]([C@H]([C@@H]([C@H](O8)CO)O)O)O)O"

for entry in entries:
    if entry["id"] == 13:
        print("Old: %s..." % entry["smiles"][:60])
        entry["smiles"] = CORRECT_SMI
        print("New: %s..." % CORRECT_SMI[:60])
        break

with open(benchPath, "w", encoding="utf-8") as f:
    json.dump(entries, f, indent=2, ensure_ascii=False)
print("Fixed Saikosaponin A SMILES (CID 167928)")
