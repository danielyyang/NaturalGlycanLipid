"""
Glycan Modification Scanner — Test + Pipeline Integration
Tests SMARTS scanning on known modified sugars, then runs on full dataset
"""
import os, sys, time
import pandas as pd
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from rdkit import Chem
from lib.modification_scanner import scanGlycanModifications, scanAndFormat
from lib.glycan_topology import find_mapped_sugar_units
from lib.bond_cleavage_engine import cleaveWithConservation


def header(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")


def testKnownModifications():
    """Test with known modified sugar SMILES [TEST DATA ONLY]"""
    header("Test: Known Modified Sugars")

    testCases = [
        ("GlcNAc (N-Acetylglucosamine)",
         "OC[C@H]1OC(O)[C@H](NC(=O)C)[C@@H](O)[C@@H]1O",
         ["NAc"]),

        ("6-O-Sulfate-GlcNAc (Heparan sulfate unit)",
         "O=S(=O)(O)OC[C@H]1OC(O)[C@H](NC(=O)C)[C@@H](O)[C@@H]1O",
         ["NAc", "Sulfate"]),

        ("2,3-di-O-Acetyl-Glucose",
         "OC[C@H]1OC(O)[C@H](OC(=O)C)[C@@H](OC(=O)C)[C@@H]1O",
         ["O-Ac"]),

        ("3-O-Methyl-Glucose",
         "OC[C@H]1OC(O)[C@H](OC)[C@@H](O)[C@@H]1O",
         ["O-Me"]),

        ("Glucuronic acid (GlcA)",
         "O=C(O)[C@@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
         ["COOH"]),

        ("Glucose (no modifications)",
         "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
         []),

        ("Phospho-mannose (Man-6-P)",
         "O=P(O)(O)OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O",
         ["Phosphate"]),
    ]

    allPassed = True
    for name, smiles, expectedMods in testCases:
        mods = scanGlycanModifications(smiles)
        formatted = scanAndFormat(smiles)
        detectedNames = list(mods.keys())

        # Check expected mods are present
        missing = [m for m in expectedMods if m not in detectedNames]
        status = "PASS" if not missing else "FAIL"
        if missing:
            allPassed = False

        print(f"  {name}")
        print(f"    SMILES: {smiles[:55]}...")
        print(f"    Detected: {formatted if formatted else '(none)'}")
        print(f"    Expected: {expectedMods}")
        print(f"    [{status}] {'Missing: ' + str(missing) if missing else ''}")
        print()

    return allPassed


def testRutinGlycan():
    """Test with real Rutin glycan from Phase 2 cleavage"""
    header("Test: Rutin Glycan (from Phase 2 cleavage)")
    RUTIN = "C[C@@H]1OC(OC[C@H]2OC(Oc3c(-c4ccc(O)c(O)c4)oc4cc(O)cc(O)c4c3=O)C(O)C(O)C2O)C(O)C(O)C1O"
    mol = Chem.MolFromSmiles(RUTIN)
    units = find_mapped_sugar_units(mol)
    glycan, aglycon, meta = cleaveWithConservation(mol, units)

    print(f"  Glycan SMILES: {glycan[:60]}...")
    mods = scanGlycanModifications(glycan)
    formatted = scanAndFormat(glycan)
    print(f"  Modifications: {formatted if formatted else '(none - as expected for underivatized Rutin)'}")

    # Rutin has natural Rha (deoxysugar) but no synthetic modifications
    # O-Me might match if SMARTS hits the rutinoside linkage oxygen
    print(f"  Detail: {mods}")
    return True


def testFullDataset():
    """Run on full dataset and report stats"""
    header("Full Dataset Modification Scan (94,242 rows)")

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    dataPath = os.path.join(baseDir, "reports", "Coconut_Sugar_Check.csv")

    df = pd.read_csv(dataPath, low_memory=False, dtype=str, encoding="utf-8-sig")
    print(f"  Total rows: {len(df)}")

    # Use canonical_smiles (full molecule) since we don't have pre-split Glycan column
    # For actual pipeline, Phase 2 would split first
    smilesCol = "canonical_smiles"
    validMask = df[smilesCol].notna() & (~df[smilesCol].isin(["", "nan", "NULL"]))
    validDf = df[validMask].copy()
    print(f"  Valid SMILES: {len(validDf)}")

    t0 = time.time()
    modResults = []
    modStats = {}

    for smiles in validDf[smilesCol]:
        mods = scanGlycanModifications(str(smiles))
        formatted = scanAndFormat(str(smiles))
        modResults.append(formatted)
        for modName, count in mods.items():
            modStats[modName] = modStats.get(modName, 0) + count

    elapsed = time.time() - t0
    validDf["Glycan_Modifications"] = modResults

    hasModification = sum(1 for r in modResults if r)
    print(f"  Scan time: {elapsed:.1f}s ({elapsed/len(validDf)*1000:.2f}ms/compound)")
    print(f"  With modifications: {hasModification} / {len(validDf)} ({hasModification/len(validDf)*100:.1f}%)")

    print(f"\n  [Modification Distribution — Total Matches]")
    for modName, totalCount in sorted(modStats.items(), key=lambda x: -x[1]):
        # Count how many compounds have this mod
        compoundCount = sum(1 for r in modResults if modName in r)
        print(f"    {modName:15} {totalCount:6} matches across {compoundCount:5} compounds")

    # Top 10 examples
    withMods = validDf[validDf["Glycan_Modifications"] != ""]
    if not withMods.empty:
        print(f"\n  [Sample Compounds with Modifications (top 10)]")
        for _, row in withMods.head(10).iterrows():
            smiles = str(row[smilesCol])[:40]
            mods = row["Glycan_Modifications"]
            print(f"    {mods:35}  SMILES: {smiles}...")

    return hasModification


if __name__ == "__main__":
    print("=" * 70)
    print("  Glycan Modification Scanner — Comprehensive Test Suite")
    print("=" * 70)

    unitPassed = testKnownModifications()
    testRutinGlycan()
    modCount = testFullDataset()

    print(f"\n{'='*70}")
    print(f"  Unit tests: {'ALL PASSED' if unitPassed else 'SOME FAILED'}")
    print(f"  Full dataset: {modCount} compounds with modifications")
    print(f"{'='*70}")
