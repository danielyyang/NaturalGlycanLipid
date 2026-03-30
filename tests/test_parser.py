import pandas as pd
from rdkit import Chem
import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from lib import glycan_topology

df = pd.read_csv("data/Coconut.csv", nrows=100)
df_sugar = df[df['contains_sugar'].astype(str).str.lower() == 'true'].copy()

failures = 0
successes = 0

for idx, row in df_sugar.iterrows():
    smi = row['canonical_smiles']
    mol = Chem.MolFromSmiles(smi)
    if not mol: continue
    
    aglycan, glycan, is_fp = sugar_utils.strict_cleavage(mol)
    if is_fp:
        failures += 1
        print(f"FAILED (False Positive): {smi}")
    else:
        successes += 1
        print(f"SUCCESS: {smi} | Glycan: {glycan}")

print(f"Total Success: {successes}, Total Failures (FP): {failures}")
