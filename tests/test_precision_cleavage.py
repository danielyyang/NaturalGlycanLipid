"""
End-to-end test: Precision Cleavage -> Phase 5 -> Sugar ID
Verifies carbon conservation and correct sugar identification
"""
import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from rdkit import Chem
from lib.glycan_topology import find_mapped_sugar_units
from lib.bond_cleavage_engine import cleaveWithConservation, findGlycosidicBonds, findAnomericCarbons
from lib.feature_extractor import processPhase5Row


def header(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")


def testRutin():
    header("Test 1: Rutin (Glc + Rha + Quercetin)")
    # Rutin SMILES: quercetin-3-O-[alpha-L-Rha-(1->6)-beta-D-Glc]
    RUTIN = "C[C@@H]1OC(OC[C@H]2OC(Oc3c(-c4ccc(O)c(O)c4)oc4cc(O)cc(O)c4c3=O)C(O)C(O)C2O)C(O)C(O)C1O"
    mol = Chem.MolFromSmiles(RUTIN)
    units = find_mapped_sugar_units(mol)

    print(f"  Sugar units: {len(units)}")
    for i, u in enumerate(units):
        anomeric = findAnomericCarbons(mol, u["ring_atoms"])
        print(f"    Unit {i}: {u.get('name','?')}, ring={u['ring_atoms']}, anomeric_C={anomeric}")

    # Precision cleavage
    glycan, aglycon, meta = cleaveWithConservation(mol, units)
    print(f"\n  Glycan:  {glycan}")
    print(f"  Aglycon: {aglycon}")
    print(f"  Bonds cut: {meta['bonds_cut']}")
    print(f"  Carbon: {meta['carbon_original']} = {meta['carbon_glycan']} + {meta['carbon_aglycon']}")
    print(f"  Conserved: {meta['carbon_conserved']}")

    assert meta["carbon_conserved"], f"CARBON LOST! {meta['carbon_original']} != {meta['carbon_glycan']} + {meta['carbon_aglycon']}"

    # Phase 5 re-test with REAL cleavage products
    result = processPhase5Row(glycan, aglycon)
    print(f"\n  [Phase 5 on REAL cleavage products]")
    print(f"    Sugar Sequence:         {result['sugar_sequence']}")
    print(f"    Sugar Functional Group: {result['sugar_functional_group']}")
    print(f"    Murcko Scaffold:        {result['murcko_scaffold']}")
    print(f"    Ring Count:             {result['aglycon_ring_count']}")

    # Key assertions
    seq = result["sugar_sequence"]
    # Should NOT contain "Pen" or "Xyl" (those are pentoses, Rutin has only hexoses)
    assert "Xyl" not in seq, f"FAIL: Xyl detected in Rutin! C6 was lost. seq={seq}"
    assert "Pen" not in seq, f"FAIL: Pentose detected in Rutin! C6 was lost. seq={seq}"
    print("\n  >> Carbon conservation: PASSED")
    print("  >> No pentose (Xyl/Pen) detected: PASSED")


def testGinsenosideRg1():
    header("Test 2: Ginsenoside Rg1 (2x Glc + Dammarane)")
    RG1 = "C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@]2(O)C[C@@H](C)O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4O)C[C@H](O)[C@@H]5[C@@]3(C)CC[C@H](C5(C)C)O[C@@H]6O[C@H](CO)[C@@H](O)[C@H](O)[C@H]6O"
    mol = Chem.MolFromSmiles(RG1)
    units = find_mapped_sugar_units(mol)

    print(f"  Sugar units: {len(units)}")
    glycan, aglycon, meta = cleaveWithConservation(mol, units)
    print(f"  Glycan:  {glycan[:60]}...")
    print(f"  Aglycon: {aglycon[:60]}...")
    print(f"  Carbon: {meta['carbon_original']} = {meta['carbon_glycan']} + {meta['carbon_aglycon']}")
    print(f"  Conserved: {meta['carbon_conserved']}")

    assert meta["carbon_conserved"], f"CARBON LOST!"

    result = processPhase5Row(glycan, aglycon)
    print(f"\n  [Phase 5]")
    print(f"    Sugar Sequence:  {result['sugar_sequence']}")
    print(f"    Murcko Scaffold: {result['murcko_scaffold'][:50]}")
    print("  >> PASSED")


def testSingleGlycosideFlavonoide():
    header("Test 3: Puerarin (single C-glycoside, daidzein-8-C-glucoside)")
    PUERARIN = "OC[C@H]1OC([C@@H](O)[C@@H](O)[C@@H]1O)c1c(O)cc2oc(-c3ccc(O)cc3)c(=O)c2c1O"
    mol = Chem.MolFromSmiles(PUERARIN)
    if mol is None:
        print("  Skipped (invalid SMILES)")
        return

    units = find_mapped_sugar_units(mol)
    print(f"  Sugar units: {len(units)}")

    glycan, aglycon, meta = cleaveWithConservation(mol, units)
    print(f"  Glycan:  {glycan}")
    print(f"  Aglycon: {aglycon}")
    print(f"  Carbon: {meta['carbon_original']} = {meta['carbon_glycan']} + {meta['carbon_aglycon']}")
    print(f"  Conserved: {meta['carbon_conserved']}")

    assert meta["carbon_conserved"], f"CARBON LOST!"
    print("  >> PASSED")


def testBatchSummary():
    header("Summary")
    # More test molecules
    testMols = {
        "Rutin": "C[C@@H]1OC(OC[C@H]2OC(Oc3c(-c4ccc(O)c(O)c4)oc4cc(O)cc(O)c4c3=O)C(O)C(O)C2O)C(O)C(O)C1O",
        "Salicin": "OC[C@H]1OC(Oc2ccccc2CO)[C@H](O)[C@@H](O)[C@@H]1O",
        "Amygdalin": "N#C[C@H](OC1OC(CO[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@H]1O)c1ccccc1",
    }
    allPassed = True
    for name, smi in testMols.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"  {name}: invalid SMILES")
            continue
        units = find_mapped_sugar_units(mol)
        glycan, aglycon, meta = cleaveWithConservation(mol, units)
        status = "OK" if meta["carbon_conserved"] else "FAIL"
        if not meta["carbon_conserved"]:
            allPassed = False
        print(f"  {name:15} C={meta['carbon_original']:2} -> G={meta['carbon_glycan']:2} + A={meta['carbon_aglycon']:2}  [{status}]")

    assert allPassed, "Some molecules failed carbon conservation!"
    print("\n  >> All carbon conservation checks PASSED")


if __name__ == "__main__":
    print("=" * 70)
    print("  Precision Cleavage + Phase 5 End-to-End Test")
    print("=" * 70)

    testRutin()
    testGinsenosideRg1()
    testSingleGlycosideFlavonoide()
    testBatchSummary()

    print(f"\n{'='*70}")
    print("  ALL TESTS PASSED")
    print(f"{'='*70}")
