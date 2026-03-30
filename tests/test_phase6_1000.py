"""
Phase 6 Test: 1000-sample classification with distribution stats
Uses canonical_smiles directly (full molecule SMILES).
"""
import os, sys, time
import pandas as pd
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from rdkit import Chem
from lib.chemical_classifier import (
    classifyAglycon,
    buildReferenceFingerprints,
)


def main():
    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    dataPath = os.path.join(baseDir, "reports", "Coconut_Sugar_Check.csv")

    print("=" * 70)
    print("  Phase 6: Intelligent Chemical Classification — 1000-sample Test")
    print("=" * 70)

    df = pd.read_csv(dataPath, low_memory=False, dtype=str, encoding="utf-8-sig")
    print(f"  Total rows: {len(df)}")

    smilesCol = "canonical_smiles"
    classCol = "np_classifier_superclass"

    validSmiles = df[smilesCol].notna() & (~df[smilesCol].isin(["", "nan", "NULL"]))
    missingClass = df[classCol].isna() | df[classCol].isin(["", "nan", "NULL"])

    missingPool = df[validSmiles & missingClass]
    existingPool = df[validSmiles & ~missingClass]

    nMissing = min(500, len(missingPool))
    nExisting = min(500, len(existingPool))

    sampleMissing = missingPool.sample(n=nMissing, random_state=42) if nMissing > 0 else pd.DataFrame()
    sampleExisting = existingPool.sample(n=nExisting, random_state=42)
    sampleDf = pd.concat([sampleMissing, sampleExisting])

    print(f"  Sample: {len(sampleDf)} ({nMissing} missing + {nExisting} with class)")

    # Step 1: Build Tanimoto reference
    print(f"\n[Step 1] Building Tanimoto reference from {nExisting} classified compounds...")
    t0 = time.time()
    refFps, refLabels = buildReferenceFingerprints(sampleExisting, smilesCol, classCol)
    print(f"  Reference FPs: {len(refFps)} in {time.time()-t0:.2f}s")

    # Step 2: Classify
    print(f"\n[Step 2] Running 3-tier classification on {len(sampleDf)} compounds...")
    t1 = time.time()

    results = []
    for _, row in sampleDf.iterrows():
        smiles = str(row.get(smilesCol, ""))
        existingClass = str(row.get(classCol, ""))
        r = classifyAglycon(smiles, refFps, refLabels, tanimotoThreshold=0.85)
        r["existing_class"] = existingClass if existingClass not in ("", "nan", "NULL", "None") else ""
        r["smiles_short"] = smiles[:50]
        results.append(r)

    classifyTime = time.time() - t1
    resultDf = pd.DataFrame(results)
    print(f"  Done in {classifyTime:.2f}s ({classifyTime/len(sampleDf)*1000:.1f}ms/compound)")

    # Results
    print(f"\n{'='*70}")
    print("  RESULTS: Classification Distribution")
    print(f"{'='*70}")

    classDist = resultDf["classification"].value_counts()
    print(f"\n  [All Classifications]")
    for cls, count in classDist.items():
        pct = count / len(resultDf) * 100
        bar = "#" * int(pct / 2)
        print(f"    {cls:50} {count:4} ({pct:5.1f}%) {bar}")

    methodDist = resultDf["classification_method"].value_counts()
    print(f"\n  [Method Distribution]")
    for method, count in methodDist.items():
        pct = count / len(resultDf) * 100
        print(f"    {method:20} {count:5} ({pct:5.1f}%)")

    # Glycolipid capture
    glycolipidRows = resultDf[resultDf["glycolipid_flag"] != ""]
    print(f"\n  [Glycolipid Capture]")
    print(f"    Total flagged: {len(glycolipidRows)} / {len(resultDf)} ({len(glycolipidRows)/len(resultDf)*100:.1f}%)")
    if not glycolipidRows.empty:
        glycoSubtypes = glycolipidRows["glycolipid_flag"].value_counts()
        for subtype, count in glycoSubtypes.items():
            print(f"      {subtype:45} {count:4}")
        print(f"\n    [Sample Glycolipid Captures]")
        for _, r in glycolipidRows.head(5).iterrows():
            print(f"      {r['classification']:35} flag={r['glycolipid_flag']}")
            print(f"        SMILES: {r['smiles_short']}...")

    # Tanimoto accuracy on known compounds
    withExisting = resultDf[resultDf["existing_class"] != ""]
    tanimotoHits = withExisting[withExisting["classification_method"] == "Tanimoto"]
    if not tanimotoHits.empty:
        matches = sum(
            1 for _, r in tanimotoHits.iterrows()
            if r["classification"].split(" (Tanimoto")[0] == r["existing_class"]
        )
        print(f"\n  [Tanimoto Rescue Accuracy]")
        print(f"    Correct: {matches}/{len(tanimotoHits)} ({matches/len(tanimotoHits)*100:.1f}%)")

    # For missing-class compounds: what did we classify?
    if nMissing > 0:
        missingResults = resultDf[resultDf["existing_class"] == ""]
        classified = missingResults[missingResults["classification"] != "Unclassified"]
        print(f"\n  [Missing-class Recovery]")
        print(f"    Classified: {len(classified)}/{len(missingResults)} ({len(classified)/max(len(missingResults),1)*100:.1f}%)")
        if not classified.empty:
            recoveredDist = classified["classification"].value_counts().head(10)
            for cls, count in recoveredDist.items():
                print(f"      {cls:50} {count:4}")

    print(f"\n{'='*70}")
    print(f"  Total time: {time.time()-t0:.1f}s")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
