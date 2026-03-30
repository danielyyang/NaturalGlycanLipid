"""
Quick preview of Saponin CSV structure for chart planning.
快速预览皂苷 CSV 数据结构，用于图表规划。
"""
import pandas as pd

CSV_PATH = r"d:\Glycan_Database\reports\GlycoNP_Saponin_DB_v13.csv"
df = pd.read_csv(CSV_PATH, nrows=5, low_memory=False)

print(f"Total columns: {len(df.columns)}")
print(f"\n=== ALL COLUMNS ===")
for i, col in enumerate(df.columns):
    print(f"  [{i:3d}] {col}")

print(f"\n=== SAMPLE VALUES (first row) ===")
for col in df.columns:
    val = str(df[col].iloc[0])[:80]
    print(f"  {col:40s} -> {val}")
