"""
PPT Slide Demo: Sugar Identification Pipeline Visualization
PPT 演示: 糖识别管线可视化

用知名天然产物 (Rutin, Hesperidin, Digoxin, Baicalin, Salicin, Phlorizin)
走完整管线: SMILES → 糖环检测 → 切分 → 单糖鉴定 → 序列输出
生成带着色分子图 + 流程图的 HTML 文件

[TEST DATA ONLY]
"""
import sys, os, base64, time
sys.path.insert(0, r"d:\Glycan_Database")

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol

from lib.glycan_topology import find_mapped_sugar_units
from lib.monosaccharide_identifier import analyze_glycan
from lib.bond_cleavage_engine import cleaveWithConservation
from lib.modification_scanner import scanAndFormat as scanGlycanMods


# ================================================================
# 知名天然产物测试集 (Well-known Natural Products)
# 全部可以在 PubChem 搜索到
# ================================================================
SHOWCASE_MOLECULES = [
    {
        "name": "Rutin (芦丁)",
        "pubchem_cid": "CID 5280805",
        "smiles": "O[C@@H]1[C@H](O)[C@@H](OC[C@@H]2OC(Oc3c(oc4cc(O)cc(O)c4c3=O)-c3ccc(O)c(O)c3)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](C)O1",
        "category": "Flavonoid Glycoside (黄酮苷)",
        "description": "Common flavonoid in buckwheat, citrus fruits. Contains rutinose (L-Rha + D-Glc) linked to quercetin.",
        "description_cn": "荞麦、柑橘中常见的黄酮苷。含芦丁糖 (L-Rha + D-Glc) 连接槲皮素苷元。",
    },
    {
        "name": "Hesperidin (橙皮苷)",
        "pubchem_cid": "CID 10621",
        "smiles": "O[C@@H]1[C@@H](O)[C@H](OC[C@@H]2OC(Oc3cc(O)c4C(=O)C[C@H](Oc5ccc(OC)c(O)c5)c4c3)[C@H](O)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](C)O1",
        "category": "Flavanone Glycoside (黄烷酮苷)",
        "description": "Major flavanone in citrus peel. Contains neohesperidose (L-Rha + D-Glc) linked to hesperetin.",
        "description_cn": "柑橘皮中的主要黄烷酮苷。含新橙皮糖 (L-Rha + D-Glc) 连接橙皮素。",
    },
    {
        "name": "Digoxin (地高辛)",
        "pubchem_cid": "CID 2724385",
        "smiles": "C[C@@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O[C@@H]3[C@@H](O)[C@H](O[C@H]4CC[C@@]5(C)[C@@H](CC[C@@H]5[C@H]5CC[C@]6(C)[C@@H](CC=7COC(=O)C=7)CC[C@@H]6[C@H]5O)C4)O[C@@H]3C)O[C@@H]2C)[C@@H](O)[C@H]1O",
        "category": "Cardiac Glycoside (强心苷)",
        "description": "Life-saving cardiac drug from foxglove. Contains 3 digitoxose sugars linked to a steroid aglycone.",
        "description_cn": "来自毛地黄的救命强心药。含3个洋地黄毒糖连接甾体苷元。",
    },
    {
        "name": "Baicalin (黄芩苷)",
        "pubchem_cid": "CID 64982",
        "smiles": "O[C@@H]1[C@@H](O)[C@H](Oc2cc3oc(-c4ccccc4)cc(=O)c3c(O)c2O)O[C@@H](C(=O)O)[C@@H]1O",
        "category": "Flavone Glucuronide (黄酮葡萄糖醛酸苷)",
        "description": "Major compound from Scutellaria baicalensis (黄芩). Contains D-GlcA linked to baicalein.",
        "description_cn": "黄芩的主要活性成分。含 D-葡萄糖醛酸 (D-GlcA) 连接黄芩素。",
    },
    {
        "name": "Salicin (水杨苷)",
        "pubchem_cid": "CID 439503",
        "smiles": "OC[C@H]1OC(Oc2ccccc2CO)[C@H](O)[C@@H](O)[C@@H]1O",
        "category": "Phenolic Glycoside (酚苷)",
        "description": "Aspirin precursor from willow bark. Contains a single D-Glc linked to salicyl alcohol.",
        "description_cn": "来自柳树皮的阿司匹林前体。含单个 D-Glc 连接水杨醇。",
    },
    {
        "name": "Phlorizin (根皮苷)",
        "pubchem_cid": "CID 6072",
        "smiles": "OC[C@H]1OC(Oc2cc(O)cc(O)c2C(=O)CCc2ccc(O)cc2)[C@H](O)[C@@H](O)[C@@H]1O",
        "category": "Dihydrochalcone Glycoside (二氢查尔酮苷)",
        "description": "SGLT2 inhibitor lead from apple tree. A single D-Glc linked to phloretin.",
        "description_cn": "苹果树来源的 SGLT2 抑制剂先导化合物。单个 D-Glc 连接根皮素。",
    },
]


# ================================================================
# RDKit 分子绘图工具 (Drawing Utilities)
# ================================================================

def drawMoleculeColored(mol, sugarAtomIdxs, size=(500, 350)):
    """生成带糖/苷元着色的分子 PNG → base64
    Generate colored molecule PNG with sugar (blue) / aglycon (orange) → base64
    """
    if mol is None:
        return ""
    try:
        AllChem.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        opts = drawer.drawOptions()
        opts.addStereoAnnotation = True

        allIdx = list(range(mol.GetNumAtoms()))
        colors = {}
        # 糖原子: 蓝色 (Sugar: blue)
        for i in sugarAtomIdxs:
            colors[i] = (0.25, 0.55, 0.95)
        # 苷元原子: 橙色 (Aglycon: orange)
        for i in allIdx:
            if i not in sugarAtomIdxs:
                colors[i] = (1.0, 0.65, 0.25)

        drawer.DrawMolecule(mol, highlightAtoms=allIdx,
                            highlightAtomColors=colors)
        drawer.FinishDrawing()
        return base64.b64encode(drawer.GetDrawingText()).decode("ascii")
    except Exception as e:
        print(f"  [WARN] Draw error: {e}")
        return ""


def drawPlainMol(mol, size=(350, 250)):
    """纯分子图 (Plain molecule image)"""
    if mol is None:
        return ""
    try:
        AllChem.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        opts = drawer.drawOptions()
        opts.addStereoAnnotation = True
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return base64.b64encode(drawer.GetDrawingText()).decode("ascii")
    except Exception:
        return ""


def findSugarAtoms(mol):
    """找到所有糖环原子索引 (Find all sugar ring atom indices)"""
    sugarAtoms = set()
    try:
        units = find_mapped_sugar_units(mol)
        for u in units:
            for key in ("ring_atoms", "position_map", "oxygen_map"):
                val = u.get(key)
                if isinstance(val, (list, tuple)):
                    sugarAtoms.update(val)
                elif isinstance(val, dict):
                    sugarAtoms.update(val.keys())
    except Exception:
        pass
    return sugarAtoms


# ================================================================
# 运行管线 (Run Pipeline on Each Molecule)
# ================================================================

def processMolecule(entry):
    """对单个分子运行完整管线 (Run full pipeline on one molecule)"""
    smi = entry["smiles"]
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None

    result = {**entry}

    # Step 1: 糖环检测 (Sugar ring detection)
    sugarAtoms = findSugarAtoms(mol)
    result["n_sugar_atoms"] = len(sugarAtoms)
    result["n_total_atoms"] = mol.GetNumAtoms()

    # Step 2: 着色分子图 (Colored molecule image)
    result["img_colored"] = drawMoleculeColored(mol, sugarAtoms, (520, 360))

    # Step 3: 糖苷键切分 (Glycosidic bond cleavage)
    try:
        sugarUnits = find_mapped_sugar_units(mol)
        glycanSmi, aglyconSmi, meta = cleaveWithConservation(mol, sugarUnits)
        result["glycan_smiles"] = glycanSmi if glycanSmi and glycanSmi != "NULL" else "—"
        result["aglycon_smiles"] = aglyconSmi if aglyconSmi and aglyconSmi != "NULL" else "—"
        result["cleavage_detail"] = f"bonds_cut={meta.get('bonds_cut',0)}, C_conserved={meta.get('carbon_conserved','?')}"
    except Exception as e:
        result["glycan_smiles"] = f"Error: {e}"
        result["aglycon_smiles"] = "—"
        result["cleavage_detail"] = "—"
        glycanSmi = None
        aglyconSmi = None

    # Step 3b: 片段分子图 (Fragment images)
    glycanMol = Chem.MolFromSmiles(glycanSmi) if glycanSmi else None
    aglyconMol = Chem.MolFromSmiles(aglyconSmi) if aglyconSmi else None
    result["img_glycan"] = drawPlainMol(glycanMol, (300, 220)) if glycanMol else ""
    result["img_aglycon"] = drawPlainMol(aglyconMol, (300, 220)) if aglyconMol else ""

    # Step 4: 单糖鉴定 + 序列生成 (Monosaccharide identification)
    try:
        seq, mods = analyze_glycan(smi)
        result["sugar_sequence"] = seq if seq else "NO_SUGAR"
        result["sugar_mods"] = mods if mods else "—"
    except Exception as e:
        result["sugar_sequence"] = f"Error: {e}"
        result["sugar_mods"] = "—"

    # Step 5: 修饰扫描 (Modification scanning)
    try:
        modDetail = scanGlycanMods(smi)
        result["mod_detail"] = modDetail if modDetail else "—"
    except Exception:
        result["mod_detail"] = "—"

    # Step 6: Murcko 骨架 (Murcko scaffold)
    if aglyconMol:
        try:
            scaffold = GetScaffoldForMol(aglyconMol)
            result["murcko"] = Chem.MolToSmiles(scaffold) if scaffold else "—"
        except Exception:
            result["murcko"] = "—"
    else:
        result["murcko"] = "—"

    return result


# ================================================================
# HTML 报告生成 (HTML Report Generator)
# ================================================================

def generateHtml(results, outputPath):
    """生成漂亮的 HTML 演示报告 (Generate beautiful HTML demo report)"""

    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>GlycoNP Pipeline — Sugar Identification Demo</title>
<style>
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{
    font-family: 'Segoe UI', system-ui, -apple-system, sans-serif;
    background: linear-gradient(135deg, #0a0e1a 0%, #101829 100%);
    color: #e0e6f0;
    padding: 30px;
    min-height: 100vh;
}}
h1 {{
    text-align: center;
    font-size: 2.2em;
    color: #7ec8e3;
    margin: 20px 0 5px;
    text-shadow: 0 0 30px rgba(126,200,227,0.25);
}}
.subtitle {{
    text-align: center;
    color: #8899aa;
    font-size: 0.95em;
    margin-bottom: 30px;
}}

/* ---- Pipeline Flowchart ---- */
.pipeline-box {{
    background: linear-gradient(135deg, #111827, #1a2332);
    border: 1px solid #2a3a5a;
    border-radius: 16px;
    padding: 30px;
    margin: 25px auto;
    max-width: 1100px;
}}
.pipeline-box h2 {{
    color: #7ec8e3;
    font-size: 1.3em;
    margin-bottom: 18px;
    text-align: center;
}}
.flow {{
    display: flex;
    align-items: center;
    justify-content: center;
    gap: 6px;
    flex-wrap: wrap;
    margin: 10px 0;
}}
.flow-step {{
    background: linear-gradient(135deg, #1e293b, #263548);
    border: 1px solid #3a4a6a;
    border-radius: 10px;
    padding: 12px 16px;
    text-align: center;
    min-width: 130px;
    transition: transform 0.2s, border-color 0.3s;
}}
.flow-step:hover {{ border-color: #7ec8e3; transform: translateY(-2px); }}
.flow-step .step-num {{
    font-size: 0.7em;
    color: #7ec8e3;
    font-weight: 700;
    letter-spacing: 1px;
    margin-bottom: 4px;
}}
.flow-step .step-title {{
    font-size: 0.85em;
    color: #e2e8f0;
    font-weight: 600;
}}
.flow-step .step-detail {{
    font-size: 0.7em;
    color: #8899aa;
    margin-top: 3px;
}}
.flow-arrow {{
    color: #4a6a8a;
    font-size: 1.5em;
    font-weight: bold;
}}

/* ---- Legend ---- */
.legend {{
    display: flex;
    gap: 25px;
    justify-content: center;
    margin: 15px 0;
    font-size: 0.9em;
}}
.legend-item {{ display: flex; align-items: center; gap: 6px; }}
.legend-color {{
    width: 18px;
    height: 18px;
    border-radius: 4px;
    border: 1px solid rgba(255,255,255,0.15);
}}

/* ---- Molecule Card ---- */
.mol-card {{
    background: linear-gradient(135deg, #111827, #1a2332);
    border: 1px solid #2a3a5a;
    border-radius: 16px;
    margin: 25px auto;
    max-width: 1100px;
    padding: 24px;
    transition: border-color 0.3s;
}}
.mol-card:hover {{ border-color: #7ec8e3; }}
.mol-header {{
    display: flex;
    justify-content: space-between;
    align-items: center;
    margin-bottom: 14px;
    flex-wrap: wrap;
    gap: 10px;
}}
.mol-header h3 {{
    font-size: 1.2em;
    color: #e2e8f0;
}}
.badge {{
    padding: 5px 14px;
    border-radius: 20px;
    font-size: 0.78em;
    font-weight: 600;
    background: rgba(126,200,227,0.12);
    color: #7ec8e3;
    border: 1px solid rgba(126,200,227,0.25);
}}
.badge-cid {{
    background: rgba(168,85,247,0.12);
    color: #a855f7;
    border-color: rgba(168,85,247,0.25);
}}
.mol-desc {{
    color: #94a3b8;
    font-size: 0.85em;
    margin-bottom: 14px;
    line-height: 1.5;
}}

/* ---- Images Row ---- */
.img-row {{
    display: flex;
    gap: 12px;
    margin: 14px 0;
    flex-wrap: wrap;
    justify-content: center;
}}
.img-box {{
    text-align: center;
    background: #ffffff;
    border-radius: 10px;
    padding: 8px;
    border: 1px solid #2a3a5a;
}}
.img-box img {{
    border-radius: 6px;
    max-width: 100%;
}}
.img-box figcaption {{
    font-size: 0.72em;
    color: #555;
    margin-top: 4px;
    font-weight: 600;
}}

/* ---- Result Grid ---- */
.result-grid {{
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 10px 20px;
    font-size: 0.88em;
    margin-top: 14px;
}}
.result-grid .label {{
    color: #7ec8e3;
    font-weight: 600;
    font-size: 0.82em;
}}
.result-grid .value {{
    color: #e2e8f0;
    word-break: break-all;
    font-family: 'Cascadia Code', 'Fira Code', monospace;
    font-size: 0.85em;
    background: rgba(15,23,42,0.6);
    padding: 4px 8px;
    border-radius: 5px;
}}
.result-grid .value.highlight {{
    color: #4ade80;
    font-weight: 700;
    font-size: 0.95em;
}}

/* ---- Arrow Steps ---- */
.step-arrow-row {{
    display: flex;
    align-items: center;
    gap: 8px;
    margin: 8px 0;
    flex-wrap: wrap;
}}
.step-pill {{
    display: inline-flex;
    align-items: center;
    gap: 5px;
    background: linear-gradient(135deg, #1e293b, #263548);
    border: 1px solid #3a4a6a;
    border-radius: 8px;
    padding: 6px 12px;
    font-size: 0.78em;
}}
.step-pill .sp-num {{
    color: #7ec8e3;
    font-weight: 700;
}}
.step-pill .sp-text {{
    color: #cbd5e1;
}}
.step-pill-arrow {{
    color: #4a6a8a;
    font-size: 1.2em;
}}

/* ---- Benchmark Section ---- */
.benchmark {{
    background: linear-gradient(135deg, #111827, #1a2332);
    border: 1px solid #2a3a5a;
    border-radius: 16px;
    padding: 24px;
    margin: 25px auto;
    max-width: 1100px;
}}
.benchmark h2 {{
    color: #7ec8e3;
    font-size: 1.3em;
    margin-bottom: 15px;
}}
.bm-grid {{
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 14px;
    margin: 15px 0;
}}
.bm-stat {{
    background: linear-gradient(135deg, #1e293b, #263548);
    border: 1px solid #3a4a6a;
    border-radius: 10px;
    padding: 18px;
    text-align: center;
}}
.bm-stat .bm-num {{
    font-size: 2.2em;
    font-weight: 800;
    color: #4ade80;
}}
.bm-stat .bm-label {{
    font-size: 0.8em;
    color: #8899aa;
    margin-top: 5px;
}}

.footer {{
    text-align: center;
    color: #4a5568;
    font-size: 0.75em;
    margin-top: 40px;
    padding: 15px;
}}
</style>
</head>
<body>

<h1>🔬 GlycoNP Sugar Identification Pipeline</h1>
<p class="subtitle">
    Automated Glycan Recognition from SMILES — Live Demo with Natural Products<br>
    <span style="font-size:0.85em;">Generated: {timestamp}</span>
</p>

<!-- ============ Pipeline Flowchart ============ -->
<div class="pipeline-box">
    <h2>Pipeline Flowchart — From Raw SMILES to Sugar Sequence</h2>
    <div class="flow">
        <div class="flow-step">
            <div class="step-num">STEP 1</div>
            <div class="step-title">Input SMILES</div>
            <div class="step-detail">Raw molecular structure</div>
        </div>
        <div class="flow-arrow">→</div>
        <div class="flow-step">
            <div class="step-num">STEP 2</div>
            <div class="step-title">Sugar Ring Detection</div>
            <div class="step-detail">RDKit ring analysis<br>+ hydroxyl density gate</div>
        </div>
        <div class="flow-arrow">→</div>
        <div class="flow-step">
            <div class="step-num">STEP 3</div>
            <div class="step-title">Glycosidic Bond Cleavage</div>
            <div class="step-detail">Anomeric carbon → C1-O-X<br>Carbon conservation assert</div>
        </div>
        <div class="flow-arrow">→</div>
        <div class="flow-step">
            <div class="step-num">STEP 4</div>
            <div class="step-title">Monosaccharide ID</div>
            <div class="step-detail">120+ SMILES library<br>3-tier matching engine</div>
        </div>
        <div class="flow-arrow">→</div>
        <div class="flow-step">
            <div class="step-num">STEP 5</div>
            <div class="step-title">Sequence & Mods</div>
            <div class="step-detail">IUPAC-like sequence<br>+ modification annotation</div>
        </div>
    </div>

    <div class="legend" style="margin-top:18px;">
        <div class="legend-item">
            <div class="legend-color" style="background:rgba(64,143,242,0.85);"></div>
            <span>Sugar atoms (糖原子)</span>
        </div>
        <div class="legend-item">
            <div class="legend-color" style="background:rgba(255,166,64,0.85);"></div>
            <span>Aglycone atoms (苷元原子)</span>
        </div>
    </div>

    <div style="text-align:center;margin-top:12px;">
        <div style="display:inline-flex;gap:12px;flex-wrap:wrap;justify-content:center;">
            <div class="step-pill"><span class="sp-num">Library</span><span class="sp-text">120+ monosaccharide SMILES with full stereochemistry</span></div>
            <div class="step-pill"><span class="sp-num">Validation</span><span class="sp-text">150-molecule benchmark (PubChem + synthetic)</span></div>
            <div class="step-pill"><span class="sp-num">Output</span><span class="sp-text">e.g. α-D-Glc-(1→6)-α-L-Rha</span></div>
        </div>
    </div>
</div>

<!-- ============ Molecule Demos ============ -->
"""

    for r in results:
        truncSmi = r["smiles"][:120] + "..." if len(r["smiles"]) > 120 else r["smiles"]
        html += f"""
<div class="mol-card">
    <div class="mol-header">
        <h3>🧪 {r['name']}</h3>
        <div style="display:flex;gap:8px;flex-wrap:wrap;">
            <span class="badge">{r['category']}</span>
            <span class="badge badge-cid">{r['pubchem_cid']}</span>
        </div>
    </div>
    <div class="mol-desc">
        {r['description']}<br>
        <span style="color:#6b7a8a;">{r['description_cn']}</span>
    </div>

    <!-- Pipeline for this molecule -->
    <div class="step-arrow-row">
        <div class="step-pill"><span class="sp-num">①</span><span class="sp-text">Input SMILES</span></div>
        <span class="step-pill-arrow">→</span>
        <div class="step-pill"><span class="sp-num">②</span><span class="sp-text">Detect {r.get('n_sugar_atoms', 0)} sugar atoms / {r.get('n_total_atoms', 0)} total</span></div>
        <span class="step-pill-arrow">→</span>
        <div class="step-pill"><span class="sp-num">③</span><span class="sp-text">Cleave glycosidic bond</span></div>
        <span class="step-pill-arrow">→</span>
        <div class="step-pill"><span class="sp-num">④</span><span class="sp-text">Identify sugars</span></div>
        <span class="step-pill-arrow">→</span>
        <div class="step-pill" style="border-color:#4ade80;"><span class="sp-num" style="color:#4ade80;">✓</span><span class="sp-text" style="color:#4ade80;">{r.get('sugar_sequence','—')}</span></div>
    </div>

    <!-- Images -->
    <div class="img-row">
"""
        if r.get("img_colored"):
            html += f"""
        <div class="img-box">
            <img src="data:image/png;base64,{r['img_colored']}" width="500">
            <figcaption>Full Molecule — Blue = Sugar, Orange = Aglycone</figcaption>
        </div>
"""
        if r.get("img_glycan"):
            html += f"""
        <div class="img-box">
            <img src="data:image/png;base64,{r['img_glycan']}" width="280">
            <figcaption>Glycan Fragment (糖部分)</figcaption>
        </div>
"""
        if r.get("img_aglycon"):
            html += f"""
        <div class="img-box">
            <img src="data:image/png;base64,{r['img_aglycon']}" width="280">
            <figcaption>Aglycone Fragment (苷元部分)</figcaption>
        </div>
"""
        html += f"""
    </div>

    <!-- Results -->
    <div class="result-grid">
        <div><span class="label">Sugar Sequence (糖序列):</span></div>
        <div><span class="value highlight">{r.get('sugar_sequence', '—')}</span></div>

        <div><span class="label">Modifications (修饰):</span></div>
        <div><span class="value">{r.get('mod_detail', '—')}</span></div>

        <div><span class="label">Glycan SMILES:</span></div>
        <div><span class="value">{r.get('glycan_smiles', '—')[:80]}{'...' if len(r.get('glycan_smiles','')) > 80 else ''}</span></div>

        <div><span class="label">Murcko Scaffold:</span></div>
        <div><span class="value">{r.get('murcko', '—')}</span></div>
    </div>
</div>
"""

    # Benchmark stats section
    html += f"""
<!-- ============ Benchmark Validation ============ -->
<div class="benchmark">
    <h2>✅ Validation — 150-Molecule Benchmark Test Set</h2>
    <p style="color:#94a3b8;font-size:0.88em;margin-bottom:15px;">
        Our pipeline is validated against a curated benchmark of 150 molecules covering extreme edge cases.
        <br>测试集覆盖 3 个难度梯度: 手工设计极限用例 + PubChem 真实生物分子 + RDKit 合成组合。
    </p>
    <div class="bm-grid">
        <div class="bm-stat">
            <div class="bm-num">18</div>
            <div class="bm-label">Boss Fights<br><span style="color:#6b7a8a;font-size:0.9em;">Hand-crafted edge cases<br>硬编码极限用例</span></div>
        </div>
        <div class="bm-stat">
            <div class="bm-num">32</div>
            <div class="bm-label">PubChem Real<br><span style="color:#6b7a8a;font-size:0.9em;">Sialyl Lewis X, Fondaparinux, Digoxin...<br>唾液酸、肝素、强心苷…</span></div>
        </div>
        <div class="bm-stat">
            <div class="bm-num">100</div>
            <div class="bm-label">Synthetic Combos<br><span style="color:#6b7a8a;font-size:0.9em;">RDKit-generated multi-sugar chains<br>合成多糖链组合</span></div>
        </div>
    </div>
    <div style="display:flex;gap:14px;margin-top:15px;flex-wrap:wrap;justify-content:center;">
        <div class="bm-stat" style="min-width:200px;">
            <div class="bm-num" style="color:#7ec8e3;">120+</div>
            <div class="bm-label">Monosaccharide Library<br>单糖参考库 SMILES 条目</div>
        </div>
        <div class="bm-stat" style="min-width:200px;">
            <div class="bm-num" style="color:#7ec8e3;">3</div>
            <div class="bm-label">Matching Tiers<br>三层退避匹配引擎</div>
        </div>
        <div class="bm-stat" style="min-width:200px;">
            <div class="bm-num" style="color:#7ec8e3;">7</div>
            <div class="bm-label">Gate Controls<br>多重门控验证</div>
        </div>
    </div>
</div>

<div class="footer">
    GlycoNP Pipeline — Sugar Identification Demo | Generated {timestamp}
</div>

</body>
</html>"""

    with open(outputPath, "w", encoding="utf-8") as f:
        f.write(html)


# ================================================================
# MAIN
# ================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("  🔬 GlycoNP Sugar Identification Pipeline Demo")
    print("  Generating HTML with real natural product examples")
    print("=" * 70)

    results = []
    for entry in SHOWCASE_MOLECULES:
        print(f"\n  Processing: {entry['name']}...")
        r = processMolecule(entry)
        if r:
            results.append(r)
            print(f"    → Sequence: {r.get('sugar_sequence', '—')}")
            print(f"    → Sugar atoms: {r.get('n_sugar_atoms', 0)} / {r.get('n_total_atoms', 0)}")
        else:
            print(f"    → FAILED")

    outPath = r"d:\Glycan_Database\reports\slide_sugar_pipeline_demo.html"
    os.makedirs(os.path.dirname(outPath), exist_ok=True)
    generateHtml(results, outPath)
    print(f"\n  ✅ HTML saved: {outPath}")
    print(f"  Total molecules: {len(results)}")
    print("=" * 70)
