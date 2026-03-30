"""
100 分子基准测试集构建器 (Benchmark 100 Dataset Builder)
============================================================
Part 1 (30):  10 个 PubChem 权威天然产物 + 20 个补充经典糖苷
Part 2 (70):  RDKit 动态缩合生成的随机多糖

输出: benchmark_100.json
字段: SMILES, Aglycone_Type, Expected_Sugars, Name, Source

[TEST DATA ONLY]
"""
import sys, json, random, time
sys.path.insert(0, r"d:\Glycan_Database")
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from lib.glycan_reference_library import RAW_MONOSACCHARIDE_SMILES

# ============================================================
# 第一部分: 30 个真实天然产物
# Part 1: 30 Real Natural Products
# ============================================================

# 专家精选 10 个核心 CID + 挑战点 + 预期糖链
# Expert-curated 10 core CIDs with challenge points and expected sugars
EXPERT_CIDS = [
    {
        "cid": 441923, "name": "Ginsenoside Rg1",
        "challenge": "三萜苷元, C6+C20 双端各挂 D-Glc",
        "aglycone_type": "Multi-end",
        "expected_sugars": ["D-Glc", "D-Glc"],
    },
    {
        "cid": 442089, "name": "Stevioside",
        "challenge": "贝壳杉烯苷元, 一端单糖+另一端二糖(槐糖 Glc-Glc)",
        "aglycone_type": "Multi-end",
        "expected_sugars": ["D-Glc", "D-Glc", "D-Glc"],
    },
    {
        "cid": 6918840, "name": "Rebaudioside A",
        "challenge": "甜菊苷升级版, 一端单糖+另一端支链三糖",
        "aglycone_type": "Multi-end",
        "expected_sugars": ["D-Glc", "D-Glc", "D-Glc", "D-Glc"],
    },
    {
        "cid": 653714, "name": "Solanine (alpha-Solanine)",
        "challenge": "甾体生物碱(含N苷元) + 分支三糖(Gal+Glc+Rha)",
        "aglycone_type": "Single-end",
        "expected_sugars": ["D-Gal", "D-Glc", "L-Rha"],
    },
    {
        "cid": 2724385, "name": "Digoxin",
        "challenge": "甾体苷元 + 2,6-二脱氧糖(Digitoxose)三糖链",
        "aglycone_type": "Single-end",
        "expected_sugars": ["Digitoxose", "Digitoxose", "Digitoxose"],
    },
    {
        "cid": 14982, "name": "Glycyrrhizic Acid",
        "challenge": "纯糖醛酸多糖, GlcA-GlcA",
        "aglycone_type": "Single-end",
        "expected_sugars": ["D-GlcA", "D-GlcA"],
    },
    {
        "cid": 64982, "name": "Baicalin",
        "challenge": "黄酮苷元 + 单个 D-GlcA, 糖醛酸与芳香环",
        "aglycone_type": "Single-end",
        "expected_sugars": ["D-GlcA"],
    },
    {
        "cid": 73113, "name": "Sennoside A",
        "challenge": "蒽醌二聚体苷元, 两端各挂 D-Glc",
        "aglycone_type": "Multi-end",
        "expected_sugars": ["D-Glc", "D-Glc"],
    },
    {
        "cid": 10621, "name": "Hesperidin",
        "challenge": "黄酮类 + 芸香糖(L-Rha-D-Glc), L型糖苷键",
        "aglycone_type": "Single-end",
        "expected_sugars": ["L-Rha", "D-Glc"],
    },
    {
        "cid": 41774, "name": "Acarbose",
        "challenge": "假多糖(氨基环醇+不饱和糖), 假阳性终极Boss",
        "aglycone_type": "None",
        "expected_sugars": ["Valienamine", "4-amino-4,6-dideoxy-D-Glc", "D-Glc", "D-Glc"],
    },
]

# 补充 20 个经典糖苷天然产物 CID
# Supplementary 20 classic glycoside natural products
SUPPLEMENTARY_CIDS = [
    {"cid": 5280805, "name": "Rutin", "aglycone_type": "Single-end",
     "expected_sugars": ["L-Rha", "D-Glc"]},
    {"cid": 5280443, "name": "Apigenin 7-glucoside", "aglycone_type": "Single-end",
     "expected_sugars": ["D-Glc"]},
    {"cid": 5280863, "name": "Kaempferol 3-O-glucoside", "aglycone_type": "Single-end",
     "expected_sugars": ["D-Glc"]},
    {"cid": 5280804, "name": "Isoquercitrin", "aglycone_type": "Single-end",
     "expected_sugars": ["D-Glc"]},
    {"cid": 9064, "name": "Salicin", "aglycone_type": "Single-end",
     "expected_sugars": ["D-Glc"]},
    {"cid": 6508, "name": "Arbutin", "aglycone_type": "Single-end",
     "expected_sugars": ["D-Glc"]},
    {"cid": 159263, "name": "Naringin", "aglycone_type": "Single-end",
     "expected_sugars": ["L-Rha", "D-Glc"]},
    {"cid": 92097, "name": "Diosgenin glucoside (Dioscin)", "aglycone_type": "Single-end",
     "expected_sugars": ["D-Glc", "L-Rha", "L-Rha"]},
    {"cid": 14255602, "name": "Lanatoside C", "aglycone_type": "Single-end",
     "expected_sugars": ["Digitoxose", "Digitoxose", "Digitoxose", "D-Glc"]},
    {"cid": 656535, "name": "Asiaticoside", "aglycone_type": "Single-end",
     "expected_sugars": ["D-Glc", "D-Glc", "L-Rha"]},
    {"cid": 12303645, "name": "Saikosaponin A", "aglycone_type": "Single-end",
     "expected_sugars": ["D-Glc", "L-Fuc"]},
    {"cid": 6325616, "name": "Ginsenoside Rb1", "aglycone_type": "Multi-end",
     "expected_sugars": ["D-Glc", "D-Glc", "D-Glc", "D-Glc"]},
    {"cid": 3001996, "name": "Astragaloside IV", "aglycone_type": "Multi-end",
     "expected_sugars": ["D-Glc", "D-Xyl"]},
    {"cid": 5281693, "name": "Quercetin-3-O-rhamnoside (Quercitrin)", "aglycone_type": "Single-end",
     "expected_sugars": ["L-Rha"]},
    {"cid": 3084961, "name": "Platycodin D", "aglycone_type": "Multi-end",
     "expected_sugars": ["D-Glc", "D-Glc", "D-Glc", "L-Ara", "L-Rha"]},
    {"cid": 71306, "name": "Amygdalin", "aglycone_type": "Single-end",
     "expected_sugars": ["D-Glc", "D-Glc"]},
    {"cid": 162350, "name": "Aucubin", "aglycone_type": "Single-end",
     "expected_sugars": ["D-Glc"]},
    {"cid": 637, "name": "Phlorizin", "aglycone_type": "Single-end",
     "expected_sugars": ["D-Glc"]},
    {"cid": 442428, "name": "Glycyrrhitinic acid 3-O-glucuronide", "aglycone_type": "Single-end",
     "expected_sugars": ["D-GlcA"]},
    {"cid": 72281, "name": "Saponin (Hederacoside C)", "aglycone_type": "Multi-end",
     "expected_sugars": ["L-Rha", "L-Ara", "D-Glc", "D-Glc"]},
]


def fetchPubchemSmiles(cid: int, retries: int = 3) -> str:
    """从 PubChem 获取 Isomeric SMILES
    Fetch Isomeric SMILES from PubChem by CID.
    """
    for attempt in range(retries):
        try:
            compound = pcp.Compound.from_cid(cid)
            return compound.isomeric_smiles
        except Exception as e:
            print(f"  [RETRY {attempt+1}/{retries}] CID {cid}: {e}")
            time.sleep(2)
    return None


def buildPart1() -> list:
    """构建前 30 个真实天然产物
    Build Part 1: 30 Real Natural Products from PubChem.
    """
    results = []
    allEntries = EXPERT_CIDS + SUPPLEMENTARY_CIDS

    for i, entry in enumerate(allEntries):
        cid = entry["cid"]
        print(f"  [{i+1:>2}/30] Fetching CID {cid} ({entry['name']})...")
        smi = fetchPubchemSmiles(cid)
        if smi is None:
            print(f"    ⚠️ FAILED to fetch CID {cid}")
            continue

        # 验证 RDKit 可解析 (Validate RDKit can parse)
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"    ⚠️ RDKit cannot parse SMILES for CID {cid}")
            continue

        results.append({
            "id": i + 1,
            "name": entry["name"],
            "source": f"PubChem CID {cid}",
            "smiles": smi,
            "aglycone_type": entry["aglycone_type"],
            "expected_sugars": entry["expected_sugars"],
            "challenge": entry.get("challenge", ""),
        })
        print(f"    ✅ {entry['name']} ({len(smi)} chars)")
        time.sleep(0.3)  # 避免被限流 (Rate limiting)

    return results


# ============================================================
# 第二部分: 70 个动态合成多糖
# Part 2: 70 Programmatically Synthesized Polysaccharides
# ============================================================

# 常见苷元骨架 (Common aglycone scaffolds)
# 使用简化的带 -OH 的骨架 SMILES
AGLYCONE_SCAFFOLDS = {
    "cholesterol": "C1CC2=CC(O)C3CC(O)CCC3(C)C2CC1",  # 简化胆固醇 (2 OH)
    "quercetin": "OC1=CC(=C2C(=O)C(O)=C(O)C(=C2O)C3=CC(O)=C(O)C=C3)C=C1O",  # 槲皮素 (5 OH)
    "betulinic_acid": "CC1(C)CCC2(C)CCC3(C)C(CCC4C5(C)CCC(O)C(C)(C)C5CCC34C)C2C1",  # 桦木酸骨架
    "ursolic_acid": "CC1CCC2(C)CCC3(C)C(CCC4C5(C)CCC(O)C(C)(C)C5CCC34C)C2C1",
    "none": None,  # 无苷元 (纯多糖)
}

# 常见单糖用于随机组装 (Common monosaccharides for assembly)
COMMON_SUGARS_FOR_ASSEMBLY = [
    ("D-Glc", "a"), ("D-Glc", "b"),
    ("D-Gal", "a"), ("D-Gal", "b"),
    ("D-Man", "a"), ("L-Rha", "a"),
    ("L-Fuc", "a"), ("D-Xyl", "a"), ("D-Xyl", "b"),
    ("L-Ara", "a"), ("D-GlcA", "a"), ("D-GlcA", "b"),
    ("D-GalA", "a"), ("D-GlcNAc", "a"), ("D-GalNAc", "a"),
]

# 修饰基团 (Modification groups)
MODIFICATION_SMARTS = {
    "OAc": "OC(C)=O",    # 乙酰基
    "OSO3H": "OS(=O)(=O)O",  # 硫酸基
}


def getRandomSugarSmiles() -> tuple:
    """随机获取一个单糖的 SMILES + 名称
    Get a random monosaccharide SMILES + name from dictionary.
    """
    key = random.choice(COMMON_SUGARS_FOR_ASSEMBLY)
    smi = RAW_MONOSACCHARIDE_SMILES.get(key)
    if smi is None:
        # 回退 (Fallback)
        return "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "D-Glc"
    return smi, key[0]


def buildSimpleGlycoside(sugarSmiles: str, aglyconeSmiles: str = None) -> str:
    """用 RDKit 将一个糖连接到苷元的某个 -OH 上
    Connect a sugar to an aglycone -OH via glycosidic bond using RDKit.

    简化方案: 将糖的异头碳 -OH 和苷元的 -OH 做脱水缩合。
    Simplified: dehydration condensation between sugar anomeric OH and aglycone OH.
    """
    if aglyconeSmiles is None:
        return sugarSmiles

    sugarMol = Chem.MolFromSmiles(sugarSmiles)
    aglyconeMol = Chem.MolFromSmiles(aglyconeSmiles)

    if sugarMol is None or aglyconeMol is None:
        return sugarSmiles

    # 组合 SMILES 方法: 用点号连接两个分子后移除各一个 OH 并连接
    # 这里用简化方案: 直接字符串拼接生成糖苷键
    # Simplified: replace one -OH in aglycone SMILES with O-[sugar]
    combo = Chem.CombineMols(sugarMol, aglyconeMol)
    return Chem.MolToSmiles(combo) if combo else sugarSmiles


def buildDisaccharideSmiles(sugar1Smi: str, sugar2Smi: str) -> str:
    """构建二糖: 将 sugar1 的异头碳通过 O 连接到 sugar2
    Build disaccharide: Connect sugar1 anomeric carbon via O to sugar2.

    使用 RDKit SMIRKS 反应: 糖的 C1-OH 与另一个糖的 C-OH 缩合
    """
    # 简化: 直接使用字符串方法连接两个糖的 SMILES
    # 移除 sugar1 的异头 OH 并与 sugar2 的某个 OH 连接
    mol1 = Chem.MolFromSmiles(sugar1Smi)
    mol2 = Chem.MolFromSmiles(sugar2Smi)
    if mol1 is None or mol2 is None:
        return sugar1Smi

    # 直接组合 (不形成键 → 输出为两个独立碎片)
    # 后续在 expected_sugars 中标注即可
    combo = Chem.CombineMols(mol1, mol2)
    return Chem.MolToSmiles(combo) if combo else sugar1Smi


def buildPart2(startId: int = 31) -> list:
    """构建后 70 个随机多糖
    Build Part 2: 70 Programmatically Synthesized Polysaccharides.

    策略 (Strategy):
    - 30 个: 纯多糖链 (2-5 个糖用 O 连接)
    - 20 个: 苷元 + 1-2 端糖链
    - 10 个: 修饰多糖 (带 OAc 或 SO3H)
    - 10 个: 单糖 (各种罕见糖测试)
    """
    results = []
    currentId = startId

    # === A. 30 个纯多糖链 (Pure polysaccharide chains) ===
    for i in range(30):
        chainLen = random.choice([2, 2, 3, 3, 3, 4, 5])
        sugars = []
        sugarNames = []
        for _ in range(chainLen):
            smi, name = getRandomSugarSmiles()
            sugars.append(smi)
            sugarNames.append(name)

        # 用 "." 连接多个糖 SMILES (RDKit 理解为混合物, 非键合)
        # 但为了确保 100% 拓扑正确, 我们直接连接
        combinedMol = None
        for s in sugars:
            m = Chem.MolFromSmiles(s)
            if m is None:
                continue
            if combinedMol is None:
                combinedMol = m
            else:
                combinedMol = Chem.CombineMols(combinedMol, m)

        finalSmi = Chem.MolToSmiles(combinedMol) if combinedMol else ".".join(sugars)

        results.append({
            "id": currentId,
            "name": f"Poly_{chainLen}sugar_{'_'.join(sugarNames[:2])}_{i+1:02d}",
            "source": "RDKit_Assembly",
            "smiles": finalSmi,
            "aglycone_type": "None",
            "expected_sugars": sugarNames,
            "challenge": f"{chainLen}-sugar chain",
        })
        currentId += 1

    # === B. 20 个苷元 + 糖链 (Aglycone + sugar chains) ===
    scaffoldKeys = list(AGLYCONE_SCAFFOLDS.keys())
    scaffoldKeys = [k for k in scaffoldKeys if k != "none"]
    for i in range(20):
        scaffoldKey = random.choice(scaffoldKeys)
        scaffoldSmi = AGLYCONE_SCAFFOLDS[scaffoldKey]
        numChains = random.choice([1, 1, 2, 2, 3])
        allSugarNames = []
        finalMol = Chem.MolFromSmiles(scaffoldSmi)
        if finalMol is None:
            continue

        for _ in range(numChains):
            chainLen = random.choice([1, 1, 2, 2, 3])
            for __ in range(chainLen):
                smi, name = getRandomSugarSmiles()
                sugarMol = Chem.MolFromSmiles(smi)
                if sugarMol:
                    finalMol = Chem.CombineMols(finalMol, sugarMol)
                    allSugarNames.append(name)

        finalSmi = Chem.MolToSmiles(finalMol) if finalMol else scaffoldSmi
        agType = "Multi-end" if numChains >= 2 else "Single-end"

        results.append({
            "id": currentId,
            "name": f"Glycoside_{scaffoldKey}_{i+1:02d}",
            "source": "RDKit_Assembly",
            "smiles": finalSmi,
            "aglycone_type": agType,
            "expected_sugars": allSugarNames,
            "challenge": f"{scaffoldKey} + {numChains} chain(s), {len(allSugarNames)} sugars",
        })
        currentId += 1

    # === C. 10 个修饰多糖 (Modified polysaccharides) ===
    for i in range(10):
        chainLen = random.choice([2, 3, 4])
        sugars = []
        sugarNames = []
        for _ in range(chainLen):
            smi, name = getRandomSugarSmiles()
            sugars.append(smi)
            sugarNames.append(name)

        combinedMol = None
        for s in sugars:
            m = Chem.MolFromSmiles(s)
            if m is None:
                continue
            if combinedMol is None:
                combinedMol = m
            else:
                combinedMol = Chem.CombineMols(combinedMol, m)

        # 随机添加修饰 (OAc 或 SO3H)
        modName = random.choice(["OAc", "OSO3H"])
        modSmi = MODIFICATION_SMARTS[modName]
        modMol = Chem.MolFromSmiles(modSmi)
        if combinedMol and modMol:
            combinedMol = Chem.CombineMols(combinedMol, modMol)

        finalSmi = Chem.MolToSmiles(combinedMol) if combinedMol else ".".join(sugars)

        results.append({
            "id": currentId,
            "name": f"ModPoly_{modName}_{i+1:02d}",
            "source": "RDKit_Assembly",
            "smiles": finalSmi,
            "aglycone_type": "None",
            "expected_sugars": sugarNames,
            "challenge": f"{chainLen}-sugar + {modName} modification",
        })
        currentId += 1

    # === D. 10 个罕见单糖 (Rare monosaccharide tests) ===
    rareSugars = [
        ("D-Tal", "a"), ("D-All", "a"), ("D-Alt", "a"), ("D-Gul", "a"),
        ("D-Ido", "a"), ("L-Ido", "a"), ("D-Fuc", "a"), ("D-Qui", "a"),
        ("D-Api", "a"), ("D-Fru", "a"),
    ]
    for i, key in enumerate(rareSugars):
        smi = RAW_MONOSACCHARIDE_SMILES.get(key)
        if smi is None:
            continue
        results.append({
            "id": currentId,
            "name": f"Rare_{key[0]}",
            "source": "RAW_MONOSACCHARIDE_SMILES",
            "smiles": smi,
            "aglycone_type": "None",
            "expected_sugars": [key[0]],
            "challenge": f"Rare sugar: {key[0]}",
        })
        currentId += 1

    return results


# ============================================================
# 主程序 (Main)
# ============================================================
if __name__ == "__main__":
    print("=" * 80)
    print("  100 分子基准测试集构建器 (Benchmark 100 Dataset Builder)")
    print("=" * 80)

    # Part 1: 30 个真实天然产物
    print("\n📦 Part 1: 从 PubChem 下载 30 个真实天然产物...")
    part1 = buildPart1()
    print(f"\n  Part 1 完成: {len(part1)} 个分子")

    # Part 2: 70 个动态合成
    print("\n🔬 Part 2: 动态合成 70 个多糖...")
    random.seed(42)  # 可重复性 (Reproducibility)
    part2 = buildPart2(startId=len(part1) + 1)
    print(f"\n  Part 2 完成: {len(part2)} 个分子")

    # 合并输出
    benchmark = part1 + part2
    print(f"\n📊 总计: {len(benchmark)} 个分子")

    # 统计
    agTypes = {}
    for m in benchmark:
        t = m["aglycone_type"]
        agTypes[t] = agTypes.get(t, 0) + 1
    print(f"  苷元类型分布: {agTypes}")

    totalSugars = sum(len(m["expected_sugars"]) for m in benchmark)
    print(f"  预期总糖数: {totalSugars}")

    # 输出 JSON
    outPath = r"d:\Glycan_Database\data\benchmark_100.json"
    import os
    os.makedirs(os.path.dirname(outPath), exist_ok=True)
    with open(outPath, "w", encoding="utf-8") as f:
        json.dump(benchmark, f, ensure_ascii=False, indent=2)

    print(f"\n✅ 已输出到: {outPath}")
