"""
深度诊断脚本 v2：追溯 #6, #7, #8 骨架的原始分子来源。
不依赖 lib/ 目录，直接用 RDKit 内置 Murcko 函数。
"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor, AllChem, rdmolops
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol, MakeScaffoldGeneric
from pathlib import Path

rdDepictor.SetPreferCoordGen(True)
OUT = Path("reports/saponin_figures/_debug_scaffolds")
OUT.mkdir(parents=True, exist_ok=True)

def makeTopologyScaffold(smi: str) -> str:
    """重现 getTopologyScaffoldSmiles 的逻辑 (无需 lib 依赖)。
    Reproduce the topology scaffold extraction logic inline.
    Steps: Murcko → replace heteroatoms with C → saturate non-aromatic bonds → strip chirality
    """
    mol = Chem.MolFromSmiles(smi)
    if not mol: return "NULL"
    try:
        scaffold = GetScaffoldForMol(mol)
        rw = Chem.RWMol(scaffold)
        # 杂原子 → C (Heteroatom → Carbon)
        for atom in rw.GetAtoms():
            if atom.GetAtomicNum() != 6:
                atom.SetAtomicNum(6)
                atom.SetFormalCharge(0)
                atom.SetNumExplicitHs(0)
        # 非芳香键饱和 (Saturate non-aromatic bonds)
        for bond in rw.GetBonds():
            if bond.GetBondType() != Chem.BondType.AROMATIC:
                bond.SetBondType(Chem.BondType.SINGLE)
        Chem.SanitizeMol(rw)
        # 去手性 (Strip chirality)
        rdmolops.RemoveStereochemistry(rw)
        return Chem.MolToSmiles(rw)
    except Exception as e:
        return f"ERROR: {e}"

# 加载数据 (Load data)
df = pd.read_csv("reports/GlycoNP_Saponin_DB_v13.csv", low_memory=False)
df = df[df["Total_Sugar_Count"] > 0]
totalMols = len(df[df["Murcko_Scaffold"].notna() & (df["Murcko_Scaffold"] != "")])

targets = {
    6: "C1CCC2C(C1)CCC1C2CCC23CCC4(CCCCC42)CCC13",
    7: "CC1CCC2CCC3C4CCC5CCCCC5C4CCC123",
    8: "C1CCC2C(C1)CCC1C2CCC2C3CCCC3CCC21",
}

for scaffoldIdx, targetSmi in targets.items():
    print(f"\n{'='*70}")
    print(f"  SCAFFOLD #{scaffoldIdx}: {targetSmi}")
    print(f"{'='*70}")
    
    subset = df[df["Murcko_Scaffold"] == targetSmi]
    pct = len(subset) / totalMols * 100
    print(f"  Total: {len(subset)} molecules ({pct:.1f}% of all)")
    
    # NP Class 分布
    npDist = subset["Detailed_NP_Class"].value_counts().head(3)
    for cls, cnt in npDist.items():
        print(f"    {cls}: {cnt}")
    
    # 取 3 个样本 (Take 3 samples)
    samples = subset[subset["Aglycon_SMILES"].notna()].head(3)
    
    for rowIdx, (_, row) in enumerate(samples.iterrows()):
        aglyconSmi = row["Aglycon_SMILES"]
        storedScaffold = row["Murcko_Scaffold"]
        name = str(row.get("name", "N/A"))[:60]
        
        print(f"\n  --- Sample {rowIdx+1}: {name} ---")
        print(f"  Aglycon_SMILES: {aglyconSmi[:90]}")
        
        aglyconMol = Chem.MolFromSmiles(aglyconSmi)
        if not aglyconMol:
            print(f"  [ERROR] Cannot parse Aglycon_SMILES"); continue
        
        nAtoms = aglyconMol.GetNumHeavyAtoms()
        nRings = aglyconMol.GetRingInfo().NumRings()
        numC = sum(1 for a in aglyconMol.GetAtoms() if a.GetAtomicNum() == 6)
        numO = sum(1 for a in aglyconMol.GetAtoms() if a.GetAtomicNum() == 8)
        oRatio = numO / max(1, numC)
        print(f"  Original: atoms={nAtoms}, rings={nRings}, C={numC}, O={numO}, O/C={oRatio:.2f}")
        
        # Step A: RDKit Murcko (保留杂原子)
        try:
            stdMurcko = GetScaffoldForMol(aglyconMol)
            stdSmi = Chem.MolToSmiles(stdMurcko)
            print(f"  Std Murcko:     rings={stdMurcko.GetRingInfo().NumRings()}, atoms={stdMurcko.GetNumHeavyAtoms()}, SMILES={stdSmi[:80]}")
        except Exception as e:
            print(f"  Std Murcko: ERROR {e}"); stdMurcko = None; stdSmi = None
        
        # Step B: Generic Murcko (杂原子→C, 但不饱和)
        try:
            genericMurcko = MakeScaffoldGeneric(stdMurcko)
            genericSmi = Chem.MolToSmiles(genericMurcko)
            print(f"  Generic Murcko: rings={genericMurcko.GetRingInfo().NumRings()}, atoms={genericMurcko.GetNumHeavyAtoms()}, SMILES={genericSmi[:80]}")
        except Exception as e:
            print(f"  Generic Murcko: ERROR {e}"); genericMurcko = None; genericSmi = None
        
        # Step C: 我们的 Topo (重新计算)
        recomputed = makeTopologyScaffold(aglyconSmi)
        print(f"  Topo Recomputed: {recomputed[:80]}")
        
        match = "✅" if storedScaffold == recomputed else "❌ MISMATCH"
        print(f"  Stored vs Recomputed: {match}")
        if storedScaffold != recomputed:
            print(f"    DB stored:  {storedScaffold[:80]}")
            print(f"    Recomputed: {recomputed[:80]}")
        
        # 渲染 4 步对比图
        mols_to_draw = []
        legends = []
        
        rdDepictor.Compute2DCoords(aglyconMol)
        mols_to_draw.append(aglyconMol)
        legends.append(f"Original\n(rings={nRings}, O/C={oRatio:.2f})")
        
        if stdMurcko:
            rdDepictor.Compute2DCoords(stdMurcko)
            mols_to_draw.append(stdMurcko)
            legends.append(f"Std Murcko\n(rings={stdMurcko.GetRingInfo().NumRings()})")
        
        if genericMurcko:
            rdDepictor.Compute2DCoords(genericMurcko)
            mols_to_draw.append(genericMurcko)
            legends.append(f"Generic Murcko\n(rings={genericMurcko.GetRingInfo().NumRings()})")
        
        if recomputed not in ("NULL", "") and not recomputed.startswith("ERROR"):
            topoM = Chem.MolFromSmiles(recomputed)
            if topoM:
                rdDepictor.Compute2DCoords(topoM)
                mols_to_draw.append(topoM)
                legends.append(f"Topo Scaffold\n(rings={topoM.GetRingInfo().NumRings()})")
        
        img = Draw.MolsToGridImage(mols_to_draw, molsPerRow=4,
                                    subImgSize=(350, 300), legends=legends)
        fname = f"scaffold{scaffoldIdx}_sample{rowIdx+1}.png"
        img.save(str(OUT / fname))
        print(f"  -> Saved: {fname}")

print(f"\nAll diagnostic images saved to: {OUT}")
