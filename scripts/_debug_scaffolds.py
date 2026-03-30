"""
诊断脚本：逐个渲染 Top-8 Murcko Scaffold 并保存为独立 PNG 文件。
Diagnostic: render each scaffold individually for visual inspection.
"""
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor
from pathlib import Path

OUT = Path("reports/saponin_figures/_debug_scaffolds")
OUT.mkdir(parents=True, exist_ok=True)

scaffolds = [
    ("1_Oleanane_5ring",    "C1CCC2C(C1)CCC1C2CCC2C3CCCCC3CCC21"),
    ("2_Dammarane_4ring",   "C1CCC2C(C1)CCC1C3CCCC3CCC21"),
    ("3_Cardenolide",       "CC1CCC(C2CCC3C2CCC2C4CCCCC4CCC23)C1"),
    ("4_Spirostane",        "C1CCC2(CC1)CC1CC3C(CCC4C5CCCCC5CCC43)C1C2"),
    ("5_Furostane",         "C1CCC2C(C1)CCC1C2CCC2C3CCCC3CC21"),
    ("6_BridgedOleanane",   "C1CCC2C(C1)CCC1C2CCC23CCC4(CCCCC42)CCC13"),
    ("7_Lanostane",         "CC1CCC2CCC3C4CCC5CCCCC5C4CCC123"),
    ("8_Lupane",            "C1CCC2C(C1)CCC1C2CCC2C3CCCC3CCC21"),
]

rdDepictor.SetPreferCoordGen(True)

for label, smi in scaffolds:
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        print(f"  [FAIL] {label}: invalid SMILES")
        continue
    nR = mol.GetRingInfo().NumRings()
    nA = mol.GetNumHeavyAtoms()
    
    # 度数为1的原子有多少个 (Count degree-1 atoms)
    deg1 = [a.GetIdx() for a in mol.GetAtoms() if a.GetDegree() == 1]
    
    print(f"{label}: rings={nR}, atoms={nA}, degree1_atoms={len(deg1)}, SMILES={smi}")
    
    # 渲染原始骨架 (Render original scaffold)
    rdDepictor.Compute2DCoords(mol)
    img = Draw.MolToImage(mol, size=(400, 350))
    img.save(str(OUT / f"{label}_original.png"))
    
    # 如果有 degree-1 原子，也渲染修剪后版本 (Render pruned version)
    if deg1:
        rw = Chem.RWMol(mol)
        while True:
            to_rm = [a.GetIdx() for a in rw.GetAtoms() if a.GetDegree() == 1]
            if not to_rm:
                break
            for idx in sorted(to_rm, reverse=True):
                rw.RemoveAtom(idx)
        try:
            Chem.SanitizeMol(rw)
            rdDepictor.Compute2DCoords(rw)
            img2 = Draw.MolToImage(rw, size=(400, 350))
            img2.save(str(OUT / f"{label}_pruned.png"))
            nR2 = rw.GetRingInfo().NumRings()
            print(f"  -> Pruned: rings={nR2}, atoms={rw.GetNumHeavyAtoms()}")
        except Exception as e:
            print(f"  -> Pruned failed: {e}")

print(f"\nAll images saved to: {OUT}")
