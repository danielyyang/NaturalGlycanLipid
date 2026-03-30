import os
import sys

sys.stdout.reconfigure(encoding='utf-8')

# 将 lib 添加至相对搜寻路径
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from rdkit import Chem
from lib.glycan_topology import find_mapped_sugar_units

# [PubChem Verified SMILES Dictionary] 
test_cases = {
    # 1. N-乙酰氨基糖 (测试 NAc 修饰提取与 Hex 拒止)
    "Alpha-D-GlcNAc": "CC(=O)N[C@H]1[C@@H](O)O[C@H](CO)[C@@H](O)[C@@H]1O",

    # 2. 脱氧糖 (测试 L-构型保持与 6-deoxy 自动识别，不被虚拟引擎误伤)
    "Alpha-L-Rhamnose": "C[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",

    # 3. 糖醛酸 (测试羧基酸性糖的 A 修饰与基底合并)
    "Beta-D-Glucuronic Acid": "O=C(O)[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",

    # 4. 卡那霉素 A 核心多氨基多糖骨架 (Kanamycin A, 极其变态的抗生素测试)
    # 包含游离氨基 (-NH2) 和脱氧链霉胺，极易导致传统引擎全面崩溃成 Hex
    "Kanamycin A": "C1[C@H]([C@@H]([C@H]([C@@H]([C@H]1N)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CN)O)O)O)O)O[C@@H]3[C@@H]([C@H]([C@@H]([C@H](O3)CO)O)N)O)N",
    
    # 5. 带有环内双键的假糖 (不饱和糖环测试 - 应被拒止)
    "Unsaturated Pyranoside (Should be rejected)": "CO[C@@H]1O[C@H](CO)[C@@H](O)C=C1"
}

print("🧬 --- 开启特种修饰糖降维打击测试 ---")
for name, smiles in test_cases.items():
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print(f"❌ 解析失败: {name}")
        continue
        
    print(f"\n🧪 测试分子: {name}")
    units = find_mapped_sugar_units(mol)
    
    if not units:
        print("  ⚠️ 未能识别出任何糖元。")
    else:
        for u in units:
            mods_str = ", ".join(u['modifications']) if u['modifications'] else "None"
            print(f"  👉 识别结果: {u['name']} (异头: {u['anomeric_config']}, 修饰: {mods_str})")
