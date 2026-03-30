import pandas as pd
df = pd.read_csv("reports/GlycoNP_Saponin_DB_v13.csv", low_memory=False)
df = df[df["Total_Sugar_Count"] > 0]
df = df[df["Murcko_Scaffold"].notna() & (df["Murcko_Scaffold"] != "")]
top8 = df["Murcko_Scaffold"].value_counts().head(8)
for i, (smi, cnt) in enumerate(top8.items(), 1):
    npClass = df[df["Murcko_Scaffold"] == smi]["Detailed_NP_Class"].dropna().value_counts().head(2)
    sample = df[df["Murcko_Scaffold"] == smi]["name"].dropna().head(2).tolist()
    print(f"#{i}: N={cnt}  SMILES={smi}")
    print(f"     NP Class: {dict(npClass)}")
    print(f"     Samples: {sample}")
    print()
