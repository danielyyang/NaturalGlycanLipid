"""Quick audit: check modification tags and bond types in Saponin CSV."""
import pandas as pd
import json
import re
from collections import Counter
from pathlib import Path

df = pd.read_csv(r"d:\Glycan_Database\reports\GlycoNP_Saponin_DB_v13.csv", low_memory=False)
df = df[pd.to_numeric(df["Total_Sugar_Count"], errors="coerce") > 0]

# --- 1. All modification tags ---
allMods = []
for mods in df["Glycan_Modifications"].dropna():
    allMods.extend(re.findall(r"\*([A-Za-z\d\-]+)", str(mods)))
print("=== Top 30 Modification Tags ===")
for tag, cnt in Counter(allMods).most_common(30):
    print(f"  {tag:20s}  {cnt:,}")

# --- 2. All bond types ---
allBonds = []
for raw in df["Glycan-Aglycone_Bond_Detail"].dropna():
    try:
        bonds = json.loads(str(raw))
        for b in bonds:
            allBonds.append(b.get("bond", ""))
    except Exception:
        pass
print(f"\n=== All Bond Types (N={len(allBonds):,}) ===")
for bond, cnt in Counter(allBonds).most_common():
    print(f"  {bond:20s}  {cnt:,}")

# --- 3. Scaffold uniqueness ---
scaffolds = df["Murcko_Scaffold"].dropna().value_counts()
print(f"\n=== Scaffold Stats ===")
print(f"  Unique scaffolds: {len(scaffolds):,}")
print(f"  Top 15:")
for scaf, cnt in scaffolds.head(15).items():
    print(f"    count={cnt:,}  SMILES={scaf[:50]}...")
