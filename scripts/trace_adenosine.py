"""
Adenosine 呋喃糖鉴定完整追踪 (Full Furanose Identification Trace)
================================================================
打印 Phase 5 中每一步: 环发现 → C1 确定 → CIP 提取 → 参考库匹配
"""
import sys, os, json
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from rdkit import Chem
from lib.glycan_topology import find_mapped_sugar_units
from lib.monosaccharide_identifier import identify_monosaccharide_v10
from lib.cip_exo_engine import (
    walkSugarRing, extractSugarFingerprint,
    matchSugarFingerprint, classifyExoSubstituent,
    getReferenceFingerprintDb
)

# Load Adenosine from benchmark
with open("data/benchmark_200.json", "r", encoding="utf-8") as f:
    entries = json.load(f)

adenosine = next(e for e in entries if e["id"] == 21)
smi = adenosine["smiles"]
mol = Chem.MolFromSmiles(smi)

print("=" * 70)
print("  ADENOSINE FURANOSE IDENTIFICATION TRACE")
print("  腺苷呋喃糖鉴定完整追踪")
print("=" * 70)
print()
print("  SMILES: %s" % smi)
print("  Heavy atoms: %d" % mol.GetNumHeavyAtoms())
print()

# === Step 1: Ring Discovery ===
print("-" * 70)
print("  Step 1: Ring Discovery (环发现)")
print("-" * 70)
ri = mol.GetRingInfo()
for rIdx, ring in enumerate(ri.AtomRings()):
    size = len(ring)
    oCount = sum(1 for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 8)
    nCount = sum(1 for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 7)
    symbols = [mol.GetAtomWithIdx(i).GetSymbol() for i in ring]
    print("  Ring#%d: %d-membered [%s] O=%d N=%d" % (rIdx, size, ",".join(symbols), oCount, nCount))
print()

# === Step 2: Sugar Unit Detection ===
print("-" * 70)
print("  Step 2: Sugar Unit Detection (find_mapped_sugar_units)")
print("-" * 70)
units = find_mapped_sugar_units(mol)
print("  Found %d sugar unit(s)" % len(units))
for u in units:
    ra = u.get("ring_atoms", [])
    name = u.get("name", "?")
    anom = u.get("anomeric_config", "?")
    mods = u.get("modifications", [])
    print("  -> name=%s anomer=%s ring_atoms=%s mods=%s" % (name, anom, ra, mods))

    # === Step 3: Ring Walk ===
    print()
    print("-" * 70)
    print("  Step 3: Ring Walk (环行走)")
    print("-" * 70)
    path = walkSugarRing(mol, ra)
    if path is None:
        print("  walkSugarRing returned None!")
        continue
    print("  Walk path (C1->C%d): %s" % (len(path), path))
    for pIdx, atomIdx in enumerate(path):
        atom = mol.GetAtomWithIdx(atomIdx)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        cip = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else "?"
        tag = atom.GetChiralTag()
        tagStr = "CW" if str(tag) == "CHI_TETRAHEDRAL_CW" else ("CCW" if str(tag) == "CHI_TETRAHEDRAL_CCW" else "?")

        # Exo neighbors
        ringSet = set(ra)
        exoNbrs = [n for n in atom.GetNeighbors() if n.GetIdx() not in ringSet]
        exoSyms = [n.GetSymbol() for n in exoNbrs]

        print("  C%d (atom#%d): CIP=%s ChiralTag=%s ExoNeighbors=%s" % (
            pIdx + 1, atomIdx, cip, tagStr, exoSyms))

    # === Step 4: Exo Classification ===
    print()
    print("-" * 70)
    print("  Step 4: Exo Classification (环外取代基)")
    print("-" * 70)
    ringSet = set(ra)
    if len(path) > 1:
        exoC2 = classifyExoSubstituent(mol, path[1], ringSet)
        print("  ExoC2 (at C2): %s" % exoC2)
    if len(path) > 3:
        exoC4 = classifyExoSubstituent(mol, path[3], ringSet)
        print("  ExoC4 (at C4, furanose terminal): %s" % exoC4)
    elif len(path) > 2:
        exoC3 = classifyExoSubstituent(mol, path[2], ringSet)
        print("  ExoC3 (at C3): %s" % exoC3)

    # === Step 5: Fingerprint Extraction ===
    print()
    print("-" * 70)
    print("  Step 5: Composite Fingerprint (联合指纹)")
    print("-" * 70)
    fp = extractSugarFingerprint(mol, ra)
    if fp is not None:
        print("  CIP:  C2=%s C3=%s C4=%s C5=%s" % (fp.cipC2, fp.cipC3, fp.cipC4, fp.cipC5))
        print("  Tag:  C2=%s C3=%s C4=%s C5=%s" % (fp.tagC2, fp.tagC3, fp.tagC4, fp.tagC5))
        print("  Exo:  C2=%s C5=%s" % (fp.exoC2, fp.exoC5))
        print("  Ring: %d" % fp.ringSize)
    else:
        print("  Fingerprint extraction returned None!")

    # === Step 6: Reference Matching ===
    print()
    print("-" * 70)
    print("  Step 6: Reference Matching (参考库匹配)")
    print("-" * 70)
    if fp is not None:
        matches = matchSugarFingerprint(fp)
        if matches:
            print("  Top 5 matches:")
            for rank, m in enumerate(matches[:5]):
                print("    #%d: %s(%s) score=%d cipCheck=%d exoC5=%s exoC2=%s" % (
                    rank + 1, m[0], m[1], m[2], m[3], m[4], m[5]))
        else:
            print("  No matches found!")

    # === Step 7: v10 Final Result ===
    print()
    print("-" * 70)
    print("  Step 7: identify_monosaccharide_v10 Final Result")
    print("-" * 70)
    result = identify_monosaccharide_v10(mol, ra)
    if isinstance(result, tuple):
        print("  Result: %s (anomer=%s)" % (result[0], result[1]))
    else:
        print("  Result: %s" % result)

print()
print("=" * 70)
print("  TRACE COMPLETE")
print("=" * 70)
