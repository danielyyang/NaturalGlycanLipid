
import pandas as pd
import os

file_path = r"d:\Lipid Database\data\processed\Coconut_Sugars_Unique.xlsx"
output_info = "sheet_info.txt"
if os.path.exists(file_path):
    try:
        xl = pd.ExcelFile(file_path)
        with open(output_info, "w") as f:
            f.write(f"Sheet names: {xl.sheet_names}\n")
            df = xl.parse(xl.sheet_names[0], nrows=1)
            f.write(f"Columns (first sheet): {df.columns.tolist()}\n")
        print(f"Info written to {output_info}")
    except Exception as e:
        print(f"Error reading excel: {e}")
else:
    print("File not found.")
