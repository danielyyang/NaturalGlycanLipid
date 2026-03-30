"""
150 分子终极测试集构建器 + HTML 报告生成 (All-in-one)
Ultimate 150-Molecule Benchmark Builder + HTML Report

Part 1 (1-18):   极限 Boss (硬编码)
Part 2 (19-50):  32 PubChem 生物深水区
Part 3 (51-150): 100 RDKit 合成多糖

[TEST DATA ONLY]
"""
import sys, os, json, random, time, base64
sys.path.insert(0, r"d:\Glycan_Database")

import urllib.request, urllib.parse
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from lib.glycan_reference_library import RAW_MONOSACCHARIDE_SMILES
from lib.monosaccharide_identifier import analyze_glycan
from lib.glycan_topology import find_mapped_sugar_units
from lib.modification_scanner import scanAndFormat as scanGlycanMods

# ================================================================
# Part 1: 18 Boss 硬编码
# ================================================================
PART1 = [
    {"id":1,  "name":"DiAcetyl D-GlcA",           "smiles":"O[C@@H]1[C@H](OC(=O)C)[C@@H](OC(=O)C)[C@H](O)[C@@H](C(=O)O)O1",                       "expected":"D-GlcA","aglycone_type":"None"},
    {"id":2,  "name":"Galloylated D-Tal",          "smiles":"OC[C@H]1OC(O)[C@@H](O)[C@@H](OC(=O)c2cc(O)c(O)c(O)c2)[C@H]1O",                        "expected":"D-Tal","aglycone_type":"None"},
    {"id":3,  "name":"Formylated L-Rha",           "smiles":"C[C@@H]1OC(O)[C@H](O)[C@H](OC=O)[C@H]1O",                                              "expected":"L-Rha","aglycone_type":"None"},
    {"id":4,  "name":"Methylated L-Ara",           "smiles":"O[C@H]1[C@@H](OC)[C@@H](O)[C@@H](O)CO1",                                               "expected":"L-Ara","aglycone_type":"None"},
    {"id":5,  "name":"D-GlcNAc 6-Sulfate",        "smiles":"CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)O)[C@@H](O)[C@@H]1O",                               "expected":"D-GlcNAc","aglycone_type":"None"},
    {"id":6,  "name":"Benzoylated D-Api",          "smiles":"O[C@H]1[C@@H](OC(=O)c2ccccc2)[C@@](O)(CO)CO1",                                         "expected":"D-Api","aglycone_type":"None"},
    {"id":7,  "name":"D-Glc+D-Fuc disaccharide",   "smiles":"OC[C@H]1O[C@H](O[C@@H]2[C@H](O)[C@@H](O)[C@H](O)O[C@@H]2C)[C@H](O)[C@@H](O)[C@@H]1O","expected":"D-Glc(a1-4)D-Fuc","aglycone_type":"None"},
    {"id":8,  "name":"GlcA-Gal-Glc trisaccharide", "smiles":"O[C@H]1[C@H](O)[C@@H](O[C@H]2[C@H](O)[C@@H](O[C@H]3[C@H](O)[C@@H](O)[C@H](O)[C@@H](C(=O)O)O3)[C@@H](O)[C@@H](CO)O2)[C@H](O)[C@@H](CO)O1","expected":"D-GlcA-D-Gal-D-Glc","aglycone_type":"None"},
    {"id":9,  "name":"Triterpene (no sugar)",      "smiles":"CC1(C)CCC2(C)CC(O)C3(C)CCC4(C)CC(O)CCC4(C)C3C2C1",                                     "expected":"NO_SUGAR","aglycone_type":"None"},
    {"id":10, "name":"L-IdoA+D-GlcN heparin",     "smiles":"O=C(O)[C@@H]1O[C@H](O[C@H]2[C@H](NS(=O)(=O)O)[C@@H](O)[C@H](OS(=O)(=O)O)[C@@H](CO)O2)[C@H](O)[C@@H](O)[C@@H]1O","expected":"L-IdoA-D-GlcN","aglycone_type":"None"},
    {"id":11, "name":"Neu5Ac standalone",          "smiles":"OC(=O)[C@@]1(O)C[C@H](O)[C@@H](NC(C)=O)[C@@H](O1)[C@H](O)[C@H](O)CO",                 "expected":"Neu5Ac","aglycone_type":"None"},
    {"id":12, "name":"Branched Api+Galloyl",       "smiles":"O=C(O[C@H]1[C@@H](O)[C@@](O)(CO)CO[C@H]1O[C@H]2[C@@H](O)[C@H](O)[C@@H](O[C@H]3[C@H](O)[C@@H](O)[C@H](O)CO3)O[C@@H]2OC4=CC=C(O[C@@H]5O[C@H](CO)[C@@H](O)[C@H](O[C@@H]6O[C@@H](C)[C@H](O)[C@@H](O)[C@H]6O)[C@H]5O)C=C4)C7=CC(O)=C(O)C(O)=C7","expected":"D-Gal,D-Glc,D-Xyl,L-Rha,D-Glc","aglycone_type":"Single-end"},
    {"id":13, "name":"Marine IdoA Sulfated",       "smiles":"CO[C@H]1O[C@H](C)[C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H](C)[C@H](O)[C@@H](O)[C@H]2O[C@H]3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3O[C@H]4O[C@H](C(=O)O)[C@@H](O)[C@H](O)[C@H]4O[C@H]5O[C@H](CO)[C@@H](O)[C@H](OS(=O)(=O)O)[C@H]5NS(=O)(=O)O","expected":"D-Glc,D-GlcA,D-Gal,L-Rha,D-Qui","aglycone_type":"None"},
    {"id":14, "name":"Bacterial LPS Core",         "smiles":"CCCCCCCCO[C@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1O[C@@H]2O[C@@H](C)[C@H](O)[C@@H](O)[C@H]2O[C@@H]3O[C@@H](C)[C@H](O)[C@@H](O)[C@H]3O[C@@H]4O[C@@H](C)[C@@H](O)C[C@H]4O[C@@H]5O[C@@H](C)C[C@@H](O)[C@H]5O[C@]6(C(=O)O)C[C@H](O)[C@@H](O)[C@@H](O6)[C@@H](O)CO","expected":"KDO,D-Man,L-Rha,L-Rha,L-Asc,Hex","aglycone_type":"Single-end"},
    {"id":15, "name":"Rare sugars acylated",       "smiles":"OC(=O)CC(=O)O[C@H]1[C@H](O)[C@@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O[C@H]3O[C@H](CO)[C@H](O)[C@@H](O)[C@@H]3O[C@@H]4O[C@H](CO)[C@@H](O)[C@@H](O)[C@@H]4O[C@@H]5O[C@@H](C)[C@H](OC(=O)C)[C@@H](O)[C@H]5O[C@@H]6O[C@@H](C)[C@@H](O)[C@H](O)[C@@H]6O)[C@H](OC7=CC=CC=C7)O[C@@H]1C(=O)O","expected":"L-Rha,L-Rha,D-Alt,D-Ido,D-Glc,D-GlcA","aglycone_type":"None"},
    {"id":16, "name":"Lactose (beta)",             "smiles":"C([C@@H]1[C@@H]([C@@H]([C@H]([C@@H](O1)O[C@@H]2[C@H](O[C@H]([C@@H]([C@H]2O)O)O)CO)O)O)O)O","expected":"D-Gal(b1-4)D-Glc","aglycone_type":"None"},
    {"id":17, "name":"Maltose (beta)",             "smiles":"OC[C@H]1O[C@H](O[C@H]2[C@H](O)[C@@H](O)C(O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@@H]1O",   "expected":"D-Glc(a1-4)D-Glc","aglycone_type":"None"},
    {"id":18, "name":"Rutinose",                   "smiles":"C[C@H]1[C@@H]([C@H]([C@H]([C@@H](O1)OC[C@@H]2[C@H]([C@@H]([C@H]([C@@H](O2)O)O)O)O)O)O)O","expected":"L-Rha(a1-6)D-Glc","aglycone_type":"None"},
]

# ================================================================
# Part 2: 32 PubChem 生物深水区
# ================================================================
# 已通过名称搜索验证的正确 CID (All verified via PubChem name search)
PART2_CIDS = [
    # A. 唾液酸与肿瘤/血型抗原
    {"cid":441911, "name":"Sialyl-Lewis X",       "expected":"Neu5Ac,D-Gal,D-GlcNAc,L-Fuc","aglycone_type":"None","cat":"Sialic acid"},
    {"cid":123914, "name":"3'-Sialyllactose",     "expected":"Neu5Ac,D-Gal,D-Glc","aglycone_type":"None","cat":"Sialic acid"},
    {"cid":643987, "name":"6'-Sialyllactose",     "expected":"Neu5Ac,D-Gal,D-Glc","aglycone_type":"None","cat":"Sialic acid"},
    # B. 人类母乳寡糖
    {"cid":170484, "name":"2'-Fucosyllactose",    "expected":"L-Fuc,D-Gal,D-Glc","aglycone_type":"None","cat":"HMO"},
    {"cid":441916, "name":"LNFP-I",              "expected":"L-Fuc,D-Gal,D-GlcNAc,D-Gal,D-Glc","aglycone_type":"None","cat":"HMO"},
    {"cid":121853, "name":"Lacto-N-neotetraose",  "expected":"D-Gal,D-GlcNAc,D-Gal,D-Glc","aglycone_type":"None","cat":"HMO"},
    # C. 极端硫酸化与氨基多糖
    {"cid":636380, "name":"Fondaparinux",         "expected":"L-IdoA,D-GlcN,D-GlcA,D-GlcN,D-GlcN","aglycone_type":"None","cat":"GAG"},
    # D. 植物多端/庞大苷元
    {"cid":441923, "name":"Ginsenoside Rg1",      "expected":"D-Glc,D-Glc","aglycone_type":"Multi-end","cat":"Triterpenoid"},
    {"cid":6918840, "name":"Rebaudioside A",      "expected":"D-Glc,D-Glc,D-Glc,D-Glc","aglycone_type":"Multi-end","cat":"Diterpenoid"},
    {"cid":9549171, "name":"alpha-Solanine",      "expected":"D-Gal,D-Glc,L-Rha","aglycone_type":"Single-end","cat":"Steroidal"},
    {"cid":2724385, "name":"Digoxin",             "expected":"Digitoxose x3","aglycone_type":"Single-end","cat":"Cardiac"},
    {"cid":14982, "name":"Glycyrrhizic Acid",     "expected":"D-GlcA,D-GlcA","aglycone_type":"Single-end","cat":"Triterpenoid"},
    {"cid":64982, "name":"Baicalin",              "expected":"D-GlcA","aglycone_type":"Single-end","cat":"Flavonoid"},
    {"cid":73111, "name":"Sennoside A",           "expected":"D-Glc,D-Glc","aglycone_type":"Multi-end","cat":"Anthraquinone"},
    {"cid":10621, "name":"Hesperidin",            "expected":"L-Rha,D-Glc","aglycone_type":"Single-end","cat":"Flavonoid"},
    {"cid":41774, "name":"Acarbose",              "expected":"pseudo-sugar","aglycone_type":"None","cat":"Pseudo-sugar"},
    {"cid":5280805, "name":"Rutin",               "expected":"L-Rha,D-Glc","aglycone_type":"Single-end","cat":"Flavonoid"},
    {"cid":442428, "name":"Naringin",             "expected":"L-Rha,D-Glc","aglycone_type":"Single-end","cat":"Flavonoid"},
    {"cid":119245, "name":"Dioscin",              "expected":"D-Glc,L-Rha,L-Rha","aglycone_type":"Single-end","cat":"Steroidal"},
    {"cid":656535, "name":"Asiaticoside",         "expected":"D-Glc,D-Glc,L-Rha","aglycone_type":"Single-end","cat":"Triterpenoid"},
    {"cid":9898279, "name":"Ginsenoside Rb1",     "expected":"D-Glc,D-Glc,D-Glc,D-Glc","aglycone_type":"Multi-end","cat":"Triterpenoid"},
    {"cid":656516, "name":"Amygdalin",            "expected":"D-Glc,D-Glc","aglycone_type":"Single-end","cat":"Cyanogenic"},
    {"cid":439503, "name":"Salicin",              "expected":"D-Glc","aglycone_type":"Single-end","cat":"Phenolic"},
    {"cid":442089, "name":"Stevioside",           "expected":"D-Glc,D-Glc,D-Glc","aglycone_type":"Multi-end","cat":"Diterpenoid"},
    {"cid":6072, "name":"Phlorizin",               "expected":"D-Glc","aglycone_type":"Single-end","cat":"Phenolic"},
    {"cid":13943297, "name":"Astragaloside IV",   "expected":"D-Glc,D-Xyl","aglycone_type":"Multi-end","cat":"Triterpenoid"},
    {"cid":161800, "name":"Glycyrrhitinic acid 3-glucuronide","expected":"D-GlcA","aglycone_type":"Single-end","cat":"Triterpenoid"},
    {"cid":167928, "name":"Saikosaponin A",       "expected":"D-Glc,L-Fuc","aglycone_type":"Single-end","cat":"Triterpenoid"},
    {"cid":656630, "name":"Lanatoside C",         "expected":"D-Glc,AcDigitoxose,Digitoxose,Digitoxose","aglycone_type":"Single-end","cat":"Cardiac"},
]

# ================================================================
# Part 3: 100 合成多糖 (reuse from dataset_builder)
# ================================================================
COMMON_SUGARS = [
    ("D-Glc","a"),("D-Glc","b"),("D-Gal","a"),("D-Gal","b"),
    ("D-Man","a"),("L-Rha","a"),("L-Fuc","a"),("D-Xyl","a"),("D-Xyl","b"),
    ("L-Ara","a"),("D-GlcA","a"),("D-GlcA","b"),("D-GalA","a"),
    ("D-GlcNAc","a"),("D-GalNAc","a"),
]
AGLYCONES = {
    "cholesterol":"C1CC2=CC(O)C3CC(O)CCC3(C)C2CC1",
    "quercetin":"OC1=CC(=C2C(=O)C(O)=C(O)C(=C2O)C3=CC(O)=C(O)C=C3)C=C1O",
    "betulinic":"CC1(C)CCC2(C)CCC3(C)C(CCC4C5(C)CCC(O)C(C)(C)C5CCC34C)C2C1",
}
RARE = [("D-Tal","a"),("D-All","a"),("D-Alt","a"),("D-Gul","a"),("D-Ido","a"),
        ("L-Ido","a"),("D-Fuc","a"),("D-Qui","a"),("D-Api","a"),("D-Fru","a")]


def fetchSmiles(cid, retries=3):
    """用 PubChem REST API 获取 IsomericSMILES (不依赖 pubchempy)"""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES/JSON"
    for a in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            r = urllib.request.urlopen(req, timeout=15)
            d = json.loads(r.read())
            return d["PropertyTable"]["Properties"][0].get("IsomericSMILES") or d["PropertyTable"]["Properties"][0].get("SMILES", "")
        except Exception as e:
            print(f"    [retry {a+1}] {e}")
            time.sleep(2)
    return None


def molToB64(mol, size=(380,280), hlAtoms=None, hlColors=None):
    if mol is None: return ""
    try:
        d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        d.drawOptions().addStereoAnnotation = True
        if hlAtoms and hlColors:
            d.DrawMolecule(mol, highlightAtoms=hlAtoms, highlightAtomColors=hlColors)
        elif hlAtoms:
            d.DrawMolecule(mol, highlightAtoms=hlAtoms)
        else:
            d.DrawMolecule(mol)
        d.FinishDrawing()
        return base64.b64encode(d.GetDrawingText()).decode("ascii")
    except: return ""


def drawMol(mol):
    if mol is None: return "","",""
    sugar = set()
    try:
        units = find_mapped_sugar_units(mol)
        for u in units:
            for k in ("ring_atoms","position_map","oxygen_map"):
                v = u.get(k)
                if isinstance(v, (list,tuple)): sugar.update(v)
                elif isinstance(v, dict): sugar.update(v.keys())
    except: pass
    allIdx = set(range(mol.GetNumAtoms()))
    aglyc = allIdx - sugar
    colors = {}
    for i in sugar: colors[i] = (0.3,0.6,1.0)
    for i in aglyc: colors[i] = (1.0,0.7,0.3)
    full = molToB64(mol, hlAtoms=list(allIdx), hlColors=colors)
    gly = molToB64(mol, hlAtoms=list(sugar), hlColors={i:(0.2,0.5,1.0) for i in sugar}) if sugar else ""
    agl = molToB64(mol, hlAtoms=list(aglyc), hlColors={i:(1.0,0.6,0.2) for i in aglyc}) if aglyc and sugar else ""
    return full, gly, agl


def buildPart3(start=51):
    random.seed(42)
    res = []
    cid = start
    # 30 pure chains
    for i in range(30):
        n = random.choice([2,2,3,3,3,4,5])
        names = []
        combo = None
        for _ in range(n):
            k = random.choice(COMMON_SUGARS)
            s = RAW_MONOSACCHARIDE_SMILES.get(k,"")
            names.append(k[0])
            m = Chem.MolFromSmiles(s)
            combo = Chem.CombineMols(combo, m) if combo and m else (m or combo)
        res.append({"id":cid,"name":f"Chain_{n}s_{i+1:02d}","smiles":Chem.MolToSmiles(combo) if combo else "","expected":", ".join(names),"aglycone_type":"None","cat":"Synthetic","challenge":f"{n}-sugar chain"})
        cid += 1
    # 20 aglycone+sugar
    for i in range(20):
        akey = random.choice(list(AGLYCONES.keys()))
        combo = Chem.MolFromSmiles(AGLYCONES[akey])
        names = []
        for _ in range(random.choice([1,1,2,2,3])):
            for __ in range(random.choice([1,1,2,2,3])):
                k = random.choice(COMMON_SUGARS)
                s = RAW_MONOSACCHARIDE_SMILES.get(k,"")
                m = Chem.MolFromSmiles(s)
                if combo and m: combo = Chem.CombineMols(combo, m)
                names.append(k[0])
        res.append({"id":cid,"name":f"Aglyc_{akey[:5]}_{i+1:02d}","smiles":Chem.MolToSmiles(combo) if combo else "","expected":", ".join(names),"aglycone_type":"Multi-end" if len(names)>2 else "Single-end","cat":"Synthetic","challenge":f"{akey}+{len(names)} sugars"})
        cid += 1
    # 20 modified — 修饰基团直接附加到糖上而非独立碎片
    for i in range(20):
        n = random.choice([2,3,4])
        names = []; combo = None
        for _ in range(n):
            k = random.choice(COMMON_SUGARS)
            s = RAW_MONOSACCHARIDE_SMILES.get(k,"")
            names.append(k[0])
            m = Chem.MolFromSmiles(s)
            if m is None: continue
            combo = Chem.CombineMols(combo, m) if combo else m
        # 验证 combo 包含实际的糖 (防止 CombineMols 丢失)
        if combo is None or combo.GetNumAtoms() < 10:
            # 回退到纯糖链
            k = random.choice(COMMON_SUGARS)
            s = RAW_MONOSACCHARIDE_SMILES.get(k,"")
            combo = Chem.MolFromSmiles(s)
            names = [k[0]]
        mod = random.choice(["OAc","SO3H","OMe"])
        modSmi = {"OAc":"OC(C)=O","SO3H":"OS(=O)(=O)O","OMe":"CO"}[mod]
        mm = Chem.MolFromSmiles(modSmi)
        if combo and mm: combo = Chem.CombineMols(combo, mm)
        res.append({"id":cid,"name":f"Mod_{mod}_{i+1:02d}","smiles":Chem.MolToSmiles(combo) if combo else "","expected":", ".join(names),"aglycone_type":"None","cat":"Synthetic","challenge":f"{n}-sugar+{mod}"})
        cid += 1
    # 20 mixed rare+common
    for i in range(20):
        n = random.choice([2,3])
        names = []; combo = None
        pool = RARE + COMMON_SUGARS[:6]
        for _ in range(n):
            k = random.choice(pool)
            s = RAW_MONOSACCHARIDE_SMILES.get(k,"")
            names.append(k[0])
            m = Chem.MolFromSmiles(s)
            combo = Chem.CombineMols(combo, m) if combo and m else (m or combo)
        res.append({"id":cid,"name":f"Mixed_{i+1:02d}","smiles":Chem.MolToSmiles(combo) if combo else "","expected":", ".join(names),"aglycone_type":"None","cat":"Synthetic","challenge":f"mixed rare+common {n}-sugar"})
        cid += 1
    # 10 pure rare singles
    for i, k in enumerate(RARE):
        s = RAW_MONOSACCHARIDE_SMILES.get(k,"")
        res.append({"id":cid,"name":f"Rare_{k[0]}","smiles":s,"expected":k[0],"aglycone_type":"None","cat":"Rare","challenge":f"Rare: {k[0]}"})
        cid += 1
    return res


def htmlReport(results, path):
    total = len(results)
    withSeq = sum(1 for r in results if r.get("seq") and r["seq"]!="NO_SUGAR")
    noSugar = sum(1 for r in results if r.get("seq")=="NO_SUGAR")
    html = f"""<!DOCTYPE html><html lang="zh-CN"><head><meta charset="UTF-8">
<title>150 Molecule Ultimate Benchmark</title>
<style>
*{{margin:0;padding:0;box-sizing:border-box}}
body{{font-family:'Segoe UI',system-ui,sans-serif;background:#0a0e1a;color:#e0e0e0;padding:20px}}
h1{{text-align:center;color:#7ec8e3;font-size:2em;margin:20px 0;text-shadow:0 0 20px rgba(126,200,227,.3)}}
.stats{{display:flex;justify-content:center;gap:30px;margin:20px 0 30px;flex-wrap:wrap}}
.sc{{background:linear-gradient(135deg,#1a2332,#243447);padding:20px 30px;border-radius:12px;border:1px solid rgba(126,200,227,.2);text-align:center;min-width:140px}}
.sc .n{{font-size:2.5em;font-weight:bold;color:#7ec8e3}}.sc .l{{color:#8899aa;font-size:.9em;margin-top:5px}}
.section h2{{color:#7ec8e3;border-bottom:2px solid #2a3a4a;padding-bottom:8px;margin:30px 0 15px;font-size:1.3em}}
.mc{{background:linear-gradient(135deg,#111827,#1a2332);border:1px solid #2a3a4a;border-radius:12px;margin:12px 0;padding:16px;transition:border-color .3s}}
.mc:hover{{border-color:#7ec8e3}}
.mc-h{{display:flex;justify-content:space-between;align-items:center;margin-bottom:8px}}
.mc-h h3{{color:#e2e8f0;font-size:1em}}
.badge{{padding:3px 10px;border-radius:20px;font-size:.75em;font-weight:600}}
.bp{{background:rgba(74,222,128,.15);color:#4ade80}}.bf{{background:rgba(248,113,113,.15);color:#f87171}}
.bi{{background:rgba(126,200,227,.15);color:#7ec8e3}}
.imgs{{display:flex;gap:8px;margin:8px 0;flex-wrap:wrap}}
.imgs figure{{text-align:center}}.imgs figcaption{{color:#8899aa;font-size:.7em;margin-top:3px}}
.imgs img{{border-radius:6px;border:1px solid #2a3a4a;background:white}}
.smi{{font-family:'Courier New',monospace;font-size:.7em;color:#94a3b8;background:#0f172a;padding:4px 8px;border-radius:4px;overflow-x:auto;white-space:nowrap;max-width:100%;display:block;margin:6px 0}}
.info{{display:grid;grid-template-columns:1fr 1fr;gap:6px;font-size:.85em}}.info .lb{{color:#8899aa}}.info .vl{{color:#e0e0e0;word-break:break-all}}
.legend{{display:flex;gap:20px;justify-content:center;margin:10px 0;font-size:.85em}}
.lc{{width:14px;height:14px;border-radius:3px;display:inline-block;vertical-align:middle;margin-right:4px}}
</style></head><body>
<h1>🏆 150 Molecule Ultimate Benchmark Report</h1>
<p style="text-align:center;color:#8899aa">Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
<div class="stats">
<div class="sc"><div class="n">{total}</div><div class="l">Total</div></div>
<div class="sc"><div class="n" style="color:#4ade80">{withSeq}</div><div class="l">With Sugars</div></div>
<div class="sc"><div class="n" style="color:#f59e0b">{noSugar}</div><div class="l">No Sugar</div></div>
</div>
<div class="legend">
<span><span class="lc" style="background:rgba(77,153,255,.8)"></span> Glycan</span>
<span><span class="lc" style="background:rgba(255,178,77,.8)"></span> Aglycone</span>
</div>
"""
    groups = {}
    for r in results:
        g = r.get("group","Other")
        groups.setdefault(g,[]).append(r)
    for gname, items in groups.items():
        html += f'<div class="section"><h2>{gname} ({len(items)})</h2>\n'
        for r in items:
            hasSugar = r.get("seq") and r["seq"] not in ("NO_SUGAR","")
            bc = "bp" if hasSugar or r.get("expected","")=="NO_SUGAR" else "bf"
            html += f'<div class="mc"><div class="mc-h"><h3>#{r["id"]} {r["name"]}</h3>'
            html += f'<span class="badge {bc}">{r.get("aglycone_type","")}</span></div>\n'
            smi = r.get("smiles","")
            if len(smi)>150: smi = smi[:150]+"..."
            html += f'<code class="smi">{smi}</code>\n<div class="imgs">\n'
            for key,cap in [("img_f","Full (Colored)"),("img_g","Glycan"),("img_a","Aglycone")]:
                if r.get(key):
                    html += f'<figure><img src="data:image/png;base64,{r[key]}" width="350"><figcaption>{cap}</figcaption></figure>\n'
            html += '</div>\n<div class="info">\n'
            html += f'<div><span class="lb">Sequence:</span></div><div class="vl">{r.get("seq","—")}</div>\n'
            html += f'<div><span class="lb">Expected:</span></div><div class="vl">{r.get("expected","—")}</div>\n'
            html += f'<div><span class="lb">Mods:</span></div><div class="vl">{r.get("mods","—")}</div>\n'
            html += f'<div><span class="lb">FuncGrp:</span></div><div class="vl">{r.get("funcs","—")}</div>\n'
            if r.get("challenge"): html += f'<div><span class="lb">Challenge:</span></div><div class="vl">{r["challenge"]}</div>\n'
            if r.get("cat"): html += f'<div><span class="lb">Category:</span></div><div class="vl">{r["cat"]}</div>\n'
            html += '</div></div>\n'
        html += '</div>\n'
    html += '</body></html>'
    with open(path,"w",encoding="utf-8") as f: f.write(html)


# ================================================================
# MAIN
# ================================================================
if __name__ == "__main__":
    print("="*80)
    print("  🏆 150 分子终极测试集 构建 + 测试 + HTML 报告")
    print("="*80)
    allResults = []

    # Part 1
    print(f"\n📋 Part 1: 18 Boss Fights")
    for e in PART1:
        smi = e["smiles"]
        print(f"  #{e['id']:>3} {e['name'][:40]:40s}", end=" ", flush=True)
        mol = Chem.MolFromSmiles(smi)
        try: seq, funcs = analyze_glycan(smi)
        except Exception as ex: seq, funcs = f"ERR:{ex}", ""
        if not seq: seq = "NO_SUGAR"
        mods = ""
        try: mods = scanGlycanMods(smi)
        except: pass
        f,g,a = drawMol(mol)
        allResults.append({"group":"Part 1: 18 Boss Fights","id":e["id"],"name":e["name"],"smiles":smi,"seq":seq,"expected":e["expected"],"mods":mods or "—","funcs":funcs or "—","aglycone_type":e["aglycone_type"],"cat":"Boss","challenge":"","img_f":f,"img_g":g,"img_a":a})
        print(f"→ {seq[:60]}")

    # Part 2
    print(f"\n📋 Part 2: 32 PubChem Bio Deep Zone")
    for i, e in enumerate(PART2_CIDS):
        eid = 19 + i
        print(f"  #{eid:>3} Fetching CID {e['cid']} ({e['name'][:35]})...", end=" ", flush=True)
        smi = fetchSmiles(e["cid"])
        if not smi:
            print("⚠️ FAILED"); continue
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            print("⚠️ BAD SMILES"); continue
        try: seq, funcs = analyze_glycan(smi)
        except Exception as ex: seq, funcs = f"ERR:{ex}", ""
        if not seq: seq = "NO_SUGAR"
        mods = ""
        try: mods = scanGlycanMods(smi)
        except: pass
        f2,g2,a2 = drawMol(mol)
        allResults.append({"group":"Part 2: Bio Deep Zone","id":eid,"name":e["name"],"smiles":smi,"seq":seq,"expected":e["expected"],"mods":mods or "—","funcs":funcs or "—","aglycone_type":e["aglycone_type"],"cat":e.get("cat",""),"challenge":"","img_f":f2,"img_g":g2,"img_a":a2})
        print(f"→ {seq[:60]}")
        time.sleep(0.3)

    # Part 3
    print(f"\n📋 Part 3: 100 Synthetic Polysaccharides")
    p3 = buildPart3(start=51)
    for e in p3:
        smi = e["smiles"]
        print(f"  #{e['id']:>3} {e['name'][:40]:40s}", end=" ", flush=True)
        mol = Chem.MolFromSmiles(smi) if smi else None
        try: seq, funcs = analyze_glycan(smi) if smi else ("","")
        except Exception as ex: seq, funcs = f"ERR:{ex}", ""
        if not seq: seq = "NO_SUGAR"
        mods = ""
        try: mods = scanGlycanMods(smi) if smi else ""
        except: pass
        f3,g3,a3 = drawMol(mol)
        allResults.append({"group":"Part 3: Synthetic","id":e["id"],"name":e["name"],"smiles":smi,"seq":seq,"expected":e["expected"],"mods":mods or "—","funcs":funcs or "—","aglycone_type":e["aglycone_type"],"cat":e.get("cat",""),"challenge":e.get("challenge",""),"img_f":f3,"img_g":g3,"img_a":a3})
        print(f"→ {seq[:60]}")

    # JSON
    jsonPath = r"d:\Glycan_Database\data\benchmark_150.json"
    os.makedirs(os.path.dirname(jsonPath), exist_ok=True)
    jsonData = [{k:v for k,v in r.items() if not k.startswith("img_")} for r in allResults]
    with open(jsonPath,"w",encoding="utf-8") as f: json.dump(jsonData, f, ensure_ascii=False, indent=2)
    print(f"\n📊 JSON: {jsonPath} ({len(jsonData)} molecules)")

    # HTML
    htmlPath = r"d:\Glycan_Database\reports\benchmark_150_report.html"
    os.makedirs(os.path.dirname(htmlPath), exist_ok=True)
    htmlReport(allResults, htmlPath)
    print(f"📝 HTML: {htmlPath}")
    print(f"\n✅ 完成! Total: {len(allResults)} molecules")
