"""
查找 Scaffold #6 和 #7 对应的具体天然产物名称。
Find specific named natural products for Scaffolds #6 and #7.
"""
import pandas as pd

df = pd.read_csv("reports/GlycoNP_Saponin_DB_v13.csv", low_memory=False)
df = df[df["Total_Sugar_Count"] > 0]

targets = {
    6: "C1CCC2C(C1)CCC1C2CCC23CCC4(CCCCC42)CCC13",
    7: "CC1CCC2CCC3C4CCC5CCCCC5C4CCC123",
}

for idx, smi in targets.items():
    print(f"\n{'='*70}")
    print(f"  SCAFFOLD #{idx}: {smi}")
    print(f"{'='*70}")
    
    subset = df[df["Murcko_Scaffold"] == smi]
    
    # 找有名字的化合物 (Filter compounds with actual names, not IUPAC)
    named = subset[subset["name"].notna()].copy()
    # 排除过长的 IUPAC 名 (Exclude long IUPAC names)
    named = named[named["name"].str.len() < 40]
    
    print(f"\n  === 已知化合物名 (N={len(named)}) ===")
    nameList = named["name"].value_counts().head(20)
    for name, cnt in nameList.items():
        print(f"    {name} (×{cnt})")
    
    # NP Class 详细分布
    print(f"\n  === NP Class 分布 ===")
    for cls, cnt in subset["Detailed_NP_Class"].value_counts().head(5).items():
        print(f"    {cls}: {cnt}")
    
    # 来源科属 (Source families)
    print(f"\n  === 来源植物科 Top-5 ===")
    families = subset["LOTUS_family"].dropna().value_counts().head(5)
    for fam, cnt in families.items():
        print(f"    {fam}: {cnt}")
    
    # 来源生物 (Source organisms)
    print(f"\n  === 来源生物 Top-5 ===")
    orgs = subset["organisms"].dropna().str.split("|").explode().str.strip().value_counts().head(5)
    for org, cnt in orgs.items():
        print(f"    {org}: {cnt}")
