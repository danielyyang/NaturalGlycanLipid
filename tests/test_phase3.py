"""
Phase 3 Test: Nucleotide + Peptide detection with anti-false-positive validation
[TEST DATA ONLY] — All SMILES are from public databases
"""
import sys, os, time
sys.path.append("D:/Glycan_Database")
import pandas as pd
from lib.secondary_fragment_scanner import (
    detectNucleotideSugar, detectPeptideOrAminoAcid,
    scanSecondaryFragments, batchScanSecondary
)


def header(title):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")


def testCase(name, smiles, expectNuc, expectPep):
    """Run and validate a single test case"""
    r = scanSecondaryFragments(smiles)
    nucOk = r["Has_Nucleotide"] == expectNuc
    pepOk = r["Has_Peptide"] == expectPep
    status = "PASS" if (nucOk and pepOk) else "FAIL"
    print(f"\n  [{status}] {name}")
    print(f"    SMILES: {smiles[:70]}{'...' if len(smiles)>70 else ''}")
    print(f"    Has_Nucleotide={r['Has_Nucleotide']} (expect={expectNuc}) {'✓' if nucOk else '✗ MISMATCH'}")
    if r["Nucleotide_Detail"]:
        print(f"      Detail: {r['Nucleotide_Detail']}")
    print(f"    Has_Peptide={r['Has_Peptide']} (expect={expectPep}) {'✓' if pepOk else '✗ MISMATCH'}")
    if r["Peptide_Detail"]:
        print(f"      Detail: {r['Peptide_Detail']}")
    return nucOk and pepOk


def main():
    header("Phase 3: Nucleotide + Peptide Scanner — Test Suite")
    allPassed = True

    # ===== 核苷酸糖测试 (Nucleotide Sugar Tests) =====
    header("Nucleotide Sugar Detection")

    # 1. UDP-Glucose (尿苷二磷酸葡萄糖)
    #    来源: PubChem CID 8629 [TEST DATA ONLY]
    UDP_GLC = "OC[C@H]1OC(OP(=O)(O)OP(=O)(O)OC[C@H]2OC([C@@H](O)[C@@H]2O)n2ccc(=O)[nH]c2=O)[C@H](O)[C@@H](O)[C@@H]1O"
    allPassed &= testCase("UDP-Glucose", UDP_GLC, expectNuc=True, expectPep=False)

    # 2. GDP-Mannose (鸟苷二磷酸甘露糖)
    GDP_MAN = "OC[C@H]1OC(OP(=O)(O)OP(=O)(O)OC[C@H]2OC([C@@H](O)[C@@H]2O)n2cnc3c(=O)[nH]c(N)nc32)[C@@H](O)[C@@H](O)[C@@H]1O"
    allPassed &= testCase("GDP-Mannose", GDP_MAN, expectNuc=True, expectPep=False)

    # 3. CMP-Neu5Ac (胞苷单磷酸唾液酸)
    CMP_NEUAC = "CC(=O)N[C@@H]1[C@@H](O)CC(OP(=O)(O)OC[C@H]2OC([C@@H](O)[C@@H]2O)n2ccc(N)nc2=O)(OC1[C@@H](O)[C@@H](O)CO)C(=O)O"
    allPassed &= testCase("CMP-Neu5Ac", CMP_NEUAC, expectNuc=True, expectPep=False)

    # ===== 糖肽 / 氨基酸糖苷测试 (Peptide Tests) =====
    header("Peptide / Amino Acid Detection")

    # 4. Vancomycin (万古霉素) — 含多个肽键的糖肽抗生素
    #    简化 SMILES [TEST DATA ONLY]
    VANCOMYCIN_FRAG = "CC(O)C(NC(=O)C(CC(=O)N)NC(=O)C(NC(=O)C1CC=CC=C1)C(O)C2=CC=CC=C2)C(=O)O"
    allPassed &= testCase("Vancomycin (fragment)", VANCOMYCIN_FRAG, expectNuc=False, expectPep=True)

    # 5. Serine-O-Glucoside (丝氨酸-O-葡萄糖苷) — 氨基酸糖苷
    SER_GLC = "N[C@@H](CO[C@@H]1OC(CO)[C@@H](O)[C@H](O)[C@H]1O)C(=O)O"
    allPassed &= testCase("Serine-O-Glucoside", SER_GLC, expectNuc=False, expectPep=True)

    # ===== 防误判测试 (Anti-False-Positive Tests) =====
    header("Anti-False-Positive Controls")

    # 6. GlcNAc (N-乙酰氨基葡萄糖) — 含 NAc 但不应报 peptide
    GLCNAC = "OC[C@H]1OC(O)[C@H](NC(=O)C)[C@@H](O)[C@@H]1O"
    allPassed &= testCase("GlcNAc (NAc control)", GLCNAC, expectNuc=False, expectPep=False)

    # 7. Rutin (芦丁) — 黄酮糖苷, 无核苷酸也无肽键
    RUTIN = "C[C@@H]1OC(OC[C@H]2OC(Oc3c(-c4ccc(O)c(O)c4)oc4cc(O)cc(O)c4c3=O)C(O)C(O)C2O)C(O)C(O)C1O"
    allPassed &= testCase("Rutin (negative control)", RUTIN, expectNuc=False, expectPep=False)

    # 8. 普通苷 — 不含碱基的磷酸化糖不应报核苷酸
    PHOS_SUGAR = "O=P(O)(O)OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"
    allPassed &= testCase("Man-6-P (phosphate only, no base)", PHOS_SUGAR, expectNuc=False, expectPep=False)

    # 9. 腺嘌呤核苷 (Adenosine) — 碱基 + 核糖 但无磷酸
    ADENOSINE = "Nc1ncnc2c1ncn2[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"
    allPassed &= testCase("Adenosine (base only, no phosphate)", ADENOSINE, expectNuc=False, expectPep=False)

    # ===== 全量数据测试 (Full Dataset Test) =====
    header("Full Dataset Scan (94,242 rows)")
    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    dataPath = os.path.join(baseDir, "reports", "Coconut_Sugar_Check.csv")
    df = pd.read_csv(dataPath, low_memory=False, dtype=str, encoding="utf-8-sig")
    print(f"  Total rows: {len(df)}")

    t0 = time.time()
    df = batchScanSecondary(df, smilesCol="canonical_smiles")
    elapsed = time.time() - t0

    nucCount = df["Has_Nucleotide"].sum()
    pepCount = df["Has_Peptide"].sum()
    print(f"\n  Scan time: {elapsed:.1f}s ({elapsed/len(df)*1000:.2f}ms/compound)")
    print(f"  Nucleotide Sugars: {nucCount} ({nucCount/len(df)*100:.2f}%)")
    print(f"  Peptide/Amino Acid: {pepCount} ({pepCount/len(df)*100:.2f}%)")

    # Show examples
    nucs = df[df["Has_Nucleotide"]]
    if not nucs.empty:
        print(f"\n  [Sample Nucleotide Sugars (top 5)]")
        for _, row in nucs.head(5).iterrows():
            print(f"    Detail: {row['Nucleotide_Detail']}")
            print(f"    SMILES: {str(row['canonical_smiles'])[:60]}...")

    peps = df[df["Has_Peptide"]]
    if not peps.empty:
        print(f"\n  [Sample Peptide/Amino Acid Compounds (top 5)]")
        for _, row in peps.head(5).iterrows():
            print(f"    Detail: {row['Peptide_Detail']}")
            print(f"    SMILES: {str(row['canonical_smiles'])[:60]}...")

    # Final summary
    header("FINAL RESULT")
    print(f"  Unit tests: {'ALL 9 PASSED ✓' if allPassed else 'SOME FAILED ✗'}")
    print(f"  Nucleotide Sugars in COCONUT: {nucCount}")
    print(f"  Peptide/AA Glycosides in COCONUT: {pepCount}")
    print(f"  Total time: {elapsed:.1f}s")


if __name__ == "__main__":
    main()
