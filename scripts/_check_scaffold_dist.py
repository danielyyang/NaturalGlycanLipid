import pandas as pd
df = pd.read_csv("reports/GlycoNP_Saponin_DB_v13.csv", low_memory=False)
df = df[df["Total_Sugar_Count"] > 0]
has_scaffold = df[df["Murcko_Scaffold"].notna() & (df["Murcko_Scaffold"] != "")]
total = len(has_scaffold)
no_scaffold = len(df) - total
print(f"Total molecules: {len(df):,}")
print(f"  With Murcko_Scaffold: {total:,}")
print(f"  Without scaffold: {no_scaffold:,}")
print(f"\nUnique scaffolds: {has_scaffold['Murcko_Scaffold'].nunique()}")
top = has_scaffold["Murcko_Scaffold"].value_counts()
cumPct = 0
for i, (smi, cnt) in enumerate(top.head(20).items(), 1):
    pct = cnt / total * 100
    cumPct += pct
    print(f"  #{i}: N={cnt:5,} ({pct:5.1f}%) cumul={cumPct:5.1f}%  SMILES={smi[:50]}")
print(f"\nTop-8 cumulative: {top.head(8).sum()/total*100:.1f}%")
print(f"Top-20 cumulative: {top.head(20).sum()/total*100:.1f}%")
