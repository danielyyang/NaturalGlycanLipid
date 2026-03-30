"""Run LOTUS metadata taxonomy fill with full taxonomy columns."""
import os, sys, time, pandas as pd
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from lib.taxonomy_lotus_matcher import loadLotusDump, fillTaxonomyFromLotus

baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
lotusPath = os.path.join(baseDir, "data", "230106_frozen_metadata.csv.gz")
targetPath = os.path.join(baseDir, "reports", "Coconut_Sugar_Check.csv")
outputPath = os.path.join(baseDir, "reports", "Coconut_Sugar_FullTax.csv")

t0 = time.time()
print("=" * 70)
print("[1] Loading LOTUS metadata index (with full taxonomy)...")
lotusIndex = loadLotusDump(lotusPath)
loadTime = time.time() - t0
print(f"    Index: {len(lotusIndex)} Block-1 entries in {loadTime:.1f}s")
print(f"    Columns: {list(lotusIndex.columns)}")

print("[2] Loading COCONUT target...")
targetDf = pd.read_csv(targetPath, low_memory=False, dtype=str, encoding="utf-8-sig")
print(f"    Rows: {len(targetDf)}")

# Count missing before
def countMissing(df, col):
    if col not in df.columns: return len(df)
    return (df[col].isna() | df[col].astype(str).str.strip().isin(["", "nan", "NULL", "Not Result"])).sum()

bOrg = countMissing(targetDf, "organisms")
bFam = countMissing(targetDf, "Family")
bNpc = countMissing(targetDf, "np_classifier_superclass")

print(f"    Missing: organisms={bOrg}, Family={bFam}, np_classifier={bNpc}")

print("[3] Running merge...")
t1 = time.time()
resultDf, imputed = fillTaxonomyFromLotus(targetDf, lotusIndex)
mergeTime = time.time() - t1

aOrg = countMissing(resultDf, "organisms")
aFam = countMissing(resultDf, "Family")
aNpc = countMissing(resultDf, "np_classifier_superclass")

orgFilled = len([c for c in imputed if c[1] == "organisms"])
famFilled = len([c for c in imputed if c[1] == "Family"])
npcFilled = len([c for c in imputed if c[1] == "np_classifier_superclass"])

resultDf.to_csv(outputPath, index=False, encoding="utf-8-sig")

print()
print("=" * 70)
print("RESULTS")
print("=" * 70)
print(f"  {'':30} {'Before':>10} {'Filled':>10} {'After':>10} {'Rate':>8}")
print(f"  {'organisms':30} {bOrg:10} {orgFilled:10} {aOrg:10} {orgFilled/max(bOrg,1)*100:7.1f}%")
print(f"  {'Family':30} {bFam:10} {famFilled:10} {aFam:10} {famFilled/max(bFam,1)*100:7.1f}%")
print(f"  {'np_classifier_superclass':30} {bNpc:10} {npcFilled:10} {aNpc:10} {npcFilled/max(bNpc,1)*100:7.1f}%")
print(f"  ---")
print(f"  Index load:  {loadTime:.1f}s | Merge: {mergeTime:.3f}s | Total: {time.time()-t0:.1f}s")
print(f"  Output: {outputPath}")

# Show extra taxonomy columns
extraCols = [c for c in resultDf.columns if c.startswith("taxonomy_")]
if extraCols:
    print(f"\n  Extra taxonomy columns added: {extraCols}")

# Sample filled rows
filledIdx = [c[0] for c in imputed if c[1] == "organisms"]
if filledIdx:
    showCols = ["standard_inchi_key", "name", "organisms", "Family"]
    showCols += [c for c in extraCols if c in resultDf.columns]
    sample = resultDf.loc[resultDf.index.isin(filledIdx[:5]), [c for c in showCols if c in resultDf.columns]]
    print("\n[Sample]")
    for _, r in sample.iterrows():
        print(f"  {str(r.get('name',''))[:30]:32} org={str(r.get('organisms',''))[:40]}")
        print(f"  {'':32} fam={str(r.get('Family',''))[:25]}  kingdom={str(r.get('taxonomy_kingdom',''))[:15]}  genus={str(r.get('taxonomy_genus',''))[:20]}")
