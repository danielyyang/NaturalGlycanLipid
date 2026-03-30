"""Quick check of V13 Pruned columns and Super_Scaffold_Class distribution."""
import pandas as pd, os
REPORT_DIR = os.path.join(os.path.dirname(__file__), "..", "reports")
df = pd.read_csv(os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched_v13_Pruned.csv"), nrows=5, low_memory=False)
print("=== COLUMNS ===")
for i, c in enumerate(df.columns):
    print(f"  {i:3d}. {c}")
print(f"\n=== Total: {len(df.columns)} columns ===")

# Full load for class distribution
df2 = pd.read_csv(os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched_v13_Pruned.csv"), low_memory=False, usecols=lambda c: c in ['Super_Scaffold_Class','Organism_Type','np_classifier_superclass','Detailed_NP_Class'])
print(f"\n=== Super_Scaffold_Class distribution (top 20) ===")
if 'Super_Scaffold_Class' in df2.columns:
    print(df2['Super_Scaffold_Class'].value_counts().head(20).to_string())
print(f"\n=== Organism_Type distribution ===")
if 'Organism_Type' in df2.columns:
    print(df2['Organism_Type'].value_counts().to_string())
print(f"\n=== np_classifier_superclass distribution (top 15) ===")
if 'np_classifier_superclass' in df2.columns:
    print(df2['np_classifier_superclass'].value_counts().head(15).to_string())
