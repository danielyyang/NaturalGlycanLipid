"""
Splitter 失败尸检脚本 (Splitter Failure Autopsy)
诊断 Convallatoxin, Diosgenin glucoside, Oleandrin 的糖环检测失败
"""
import sys
import os
import json

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from rdkit import Chem
from lib.glycan_topology import find_mapped_sugar_units

# Load benchmark
with open("data/benchmark_200.json", "r", encoding="utf-8") as f:
    entries = json.load(f)

targets = {5: "Convallatoxin", 37: "Diosgenin glucoside", 39: "Oleandrin"}

for entry in entries:
    if entry["id"] not in targets:
        continue
    idx = entry["id"]
    name = entry["name"]
    smi = entry["smiles"]
    mol = Chem.MolFromSmiles(smi)

    print(f"{'='*60}")
    print(f"  [{idx}] {name}")
    print(f"{'='*60}")

    if mol is None:
        print(f"  SMILES PARSE FAILED!")
        print(f"  SMILES: {smi}")
        continue

    print(f"  Atoms: {mol.GetNumHeavyAtoms()}")
    print(f"  SMILES (first 100): {smi[:100]}...")

    # === 1. Find ALL rings ===
    ri = mol.GetRingInfo()
    allRings = ri.AtomRings()
    print(f"\n  Total rings: {len(allRings)}")

    sugarLike = []
    for rIdx, ring in enumerate(allRings):
        size = len(ring)
        if size in (5, 6):
            oCount = sum(1 for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 8)
            cCount = sum(1 for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 6)
            nCount = sum(1 for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 7)
            
            # Check exocyclic OH count
            ringSet = set(ring)
            exoOH = 0
            for atomIdx in ring:
                atom = mol.GetAtomWithIdx(atomIdx)
                if atom.GetAtomicNum() != 6:
                    continue
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() not in ringSet and nbr.GetAtomicNum() == 8:
                        exoOH += 1
            
            marker = " <-- SUGAR-LIKE" if oCount == 1 else ""
            print(f"  Ring#{rIdx}: {size}-membered  O={oCount} C={cCount} N={nCount}  exoOH={exoOH}{marker}")
            if oCount == 1:
                sugarLike.append((rIdx, ring))

    print(f"\n  Sugar-like rings (5/6-membered, O=1): {len(sugarLike)}")

    # === 2. Run find_mapped_sugar_units ===
    print(f"\n  --- find_mapped_sugar_units output ---")
    try:
        units = find_mapped_sugar_units(mol)
        print(f"  Result: {len(units)} sugar units found")
        for u in units:
            ra = u.get("ring_atoms", [])
            print(f"    ring_atoms={ra}")
    except Exception as e:
        print(f"  ERROR: {e}")

    # === 3. Check if SRU detects sugar at the ring level ===
    if sugarLike and not units:
        print(f"\n  *** DIAGNOSIS: {len(sugarLike)} sugar-like rings exist but splitter found 0! ***")
        print(f"  *** Likely cause: heuristic filter rejected the ring(s) ***")
        
        # Check each sugar-like ring individually
        for rIdx, ring in sugarLike:
            ringSet = set(ring)
            
            # Count exocyclic features
            exoO = 0
            exoC = 0
            exoN = 0
            for atomIdx in ring:
                atom = mol.GetAtomWithIdx(atomIdx)
                if atom.GetAtomicNum() != 6:
                    continue
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() not in ringSet:
                        if nbr.GetAtomicNum() == 8:
                            exoO += 1
                        elif nbr.GetAtomicNum() == 6:
                            exoC += 1
                        elif nbr.GetAtomicNum() == 7:
                            exoN += 1
            
            # Check aromatic
            isAromatic = any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)
            
            print(f"\n  Ring#{rIdx} deep analysis:")
            print(f"    Atoms: {[mol.GetAtomWithIdx(i).GetSymbol() for i in ring]}")
            print(f"    ExoO={exoO} ExoC={exoC} ExoN={exoN}")
            print(f"    Aromatic: {isAromatic}")
            print(f"    Indices: {list(ring)}")

    print()
