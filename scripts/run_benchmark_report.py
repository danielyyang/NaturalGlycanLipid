"""
118 分子综合测试 + HTML 报告生成器
Comprehensive 118-Molecule Test Runner + HTML Report Generator

包含:
  - 18 个原始测试分子 (test_corrected_18)
  - 100 个 benchmark 分子 (benchmark_100.json)
  - RDKit MolToImage 分子图 (Full/Glycan/Aglycan 带着色)
  - 完整 HTML 报告输出

[TEST DATA ONLY]
"""
import sys, os, json, base64, io, time
sys.path.insert(0, r"d:\Glycan_Database")

from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from lib.monosaccharide_identifier import analyze_glycan
from lib.glycan_topology import find_mapped_sugar_units, classify_sugar_parts
from lib.modification_scanner import scanAndFormat as scanGlycanMods

# ============================================================
# 18 个原始测试分子 (Original 18 Test Molecules)
# ============================================================
ORIGINAL_18 = [
    (1,  "DiAcetyl D-GlcA",           "O[C@@H]1[C@H](OC(=O)C)[C@@H](OC(=O)C)[C@H](O)[C@@H](C(=O)O)O1",   "D-GlcA"),
    (2,  "Galloylated D-Tal",          "OC[C@H]1OC(O)[C@@H](O)[C@@H](OC(=O)c2cc(O)c(O)c(O)c2)[C@H]1O",     "D-Tal"),
    (3,  "Formylated L-Rha",           "C[C@@H]1OC(O)[C@H](O)[C@H](OC=O)[C@H]1O",                           "L-Rha"),
    (4,  "Methylated L-Ara",           "O[C@H]1[C@@H](OC)[C@@H](O)[C@@H](O)CO1",                            "L-Ara"),
    (5,  "D-GlcNAc 6-Sulfate",        "CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)O)[C@@H](O)[C@@H]1O",           "D-GlcNAc"),
    (6,  "Benzoylated D-Api",          "O[C@H]1[C@@H](OC(=O)c2ccccc2)[C@@](O)(CO)CO1",                      "D-Api"),
    (7,  "D-Glc-D-Fuc disaccharide",   "OC[C@H]1O[C@H](O[C@@H]2[C@H](O)[C@@H](O)[C@H](O)O[C@@H]2C)[C@H](O)[C@@H](O)[C@@H]1O", "D-Glc-D-Fuc"),
    (8,  "GlcA-Gal-Glc trisaccharide", "O[C@H]1[C@H](O)[C@@H](O[C@H]2[C@H](O)[C@@H](O[C@H]3[C@H](O)[C@@H](O)[C@H](O)[C@@H](C(=O)O)O3)[C@@H](O)[C@@H](CO)O2)[C@H](O)[C@@H](CO)O1", "D-GlcA-D-Gal-D-Glc"),
    (9,  "Triterpene (no sugar)",      "CC1(C)CCC2(C)CC(O)C3(C)CCC4(C)CC(O)CCC4(C)C3C2C1",                  "NO_SUGAR"),
    (10, "L-IdoA+D-GlcN heparin",     "O=C(O)[C@@H]1O[C@H](O[C@H]2[C@H](NS(=O)(=O)O)[C@@H](O)[C@H](OS(=O)(=O)O)[C@@H](CO)O2)[C@H](O)[C@@H](O)[C@@H]1O", "L-IdoA-D-GlcN"),
    (11, "Sialic acid (Neu5Ac)",       "OC(=O)[C@@]1(O)C[C@H](O)[C@@H](NC(C)=O)[C@@H](O1)[C@H](O)[C@H](O)CO", "Neu5Ac"),
    (12, "Branched Api+Galloyl",       "O=C(O[C@H]1[C@@H](O)[C@@](O)(CO)CO[C@H]1O[C@H]2[C@@H](O)[C@H](O)[C@@H](O[C@H]3[C@H](O)[C@@H](O)[C@H](O)CO3)O[C@@H]2OC4=CC=C(O[C@@H]5O[C@H](CO)[C@@H](O)[C@H](O[C@@H]6O[C@@H](C)[C@H](O)[C@@H](O)[C@H]6O)[C@H]5O)C=C4)C7=CC(O)=C(O)C(O)=C7", "multi-sugar"),
    (13, "Marine IdoA Sulfated",       "CO[C@H]1O[C@H](C)[C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H](C)[C@H](O)[C@@H](O)[C@H]2O[C@H]3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3O[C@H]4O[C@H](C(=O)O)[C@@H](O)[C@H](O)[C@H]4O[C@H]5O[C@H](CO)[C@@H](O)[C@H](OS(=O)(=O)O)[C@H]5NS(=O)(=O)O", "5 sugars"),
    (14, "Bacterial LPS Core",         "CCCCCCCCO[C@H]1[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]1O[C@@H]2O[C@@H](C)[C@H](O)[C@@H](O)[C@H]2O[C@@H]3O[C@@H](C)[C@H](O)[C@@H](O)[C@H]3O[C@@H]4O[C@@H](C)[C@@H](O)C[C@H]4O[C@@H]5O[C@@H](C)C[C@@H](O)[C@H]5O[C@]6(C(=O)O)C[C@H](O)[C@@H](O)[C@@H](O6)[C@@H](O)CO", "6 sugars"),
    (15, "Rare sugars acylated",       "OC(=O)CC(=O)O[C@H]1[C@H](O)[C@@H](O[C@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O[C@H]3O[C@H](CO)[C@H](O)[C@@H](O)[C@@H]3O[C@@H]4O[C@H](CO)[C@@H](O)[C@@H](O)[C@@H]4O[C@@H]5O[C@@H](C)[C@H](OC(=O)C)[C@@H](O)[C@H]5O[C@@H]6O[C@@H](C)[C@@H](O)[C@H](O)[C@@H]6O)[C@H](OC7=CC=CC=C7)O[C@@H]1C(=O)O", "6 rare sugars"),
    (16, "Lactose (PubChem 6134)",     "C([C@@H]1[C@@H]([C@@H]([C@H]([C@@H](O1)O[C@@H]2[C@H](O[C@H]([C@@H]([C@H]2O)O)O)CO)O)O)O)O", "D-Gal(b1-4)D-Glc"),
    (17, "Maltose (PubChem 10991489)", "OC[C@H]1O[C@H](O[C@H]2[C@H](O)[C@@H](O)C(O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@@H]1O", "D-Glc(a1-4)D-Glc"),
    (18, "Rutinose (PubChem 441427)",  "C[C@H]1[C@@H]([C@H]([C@H]([C@@H](O1)OC[C@@H]2[C@H]([C@@H]([C@H]([C@@H](O2)O)O)O)O)O)O)O", "L-Rha(a1-6)D-Glc"),
]


def molToBase64Png(mol, size: tuple = (350, 250), highlightAtoms: list = None,
                    highlightColors: dict = None) -> str:
    """将 RDKit Mol 渲染为 base64 PNG
    Render RDKit Mol to base64 PNG string for HTML embedding.
    """
    if mol is None:
        return ""
    try:
        drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        opts = drawer.drawOptions()
        opts.addStereoAnnotation = True

        if highlightAtoms and highlightColors:
            drawer.DrawMolecule(mol, highlightAtoms=highlightAtoms,
                                highlightAtomColors=highlightColors)
        elif highlightAtoms:
            drawer.DrawMolecule(mol, highlightAtoms=highlightAtoms)
        else:
            drawer.DrawMolecule(mol)

        drawer.FinishDrawing()
        png = drawer.GetDrawingText()
        return base64.b64encode(png).decode("ascii")
    except Exception as e:
        return ""


def drawColoredMolecule(mol, smi: str) -> tuple:
    """生成带糖环着色的全分子图 + Glycan/Aglycan 片段图
    Generate full molecule image with sugar ring coloring + Glycan/Aglycan fragment images.

    Returns: (fullMolB64, glycanB64, aglycanB64)
    """
    if mol is None:
        return "", "", ""

    # 找糖环原子 (Find sugar ring atoms)
    sugarAtoms = set()
    try:
        units = find_mapped_sugar_units(mol)
        for u in units:
            for idx in u.get("ring_atoms", []):
                sugarAtoms.add(idx)
            for idx in u.get("position_map", {}).keys():
                sugarAtoms.add(idx)
            for idx in u.get("oxygen_map", {}).keys():
                sugarAtoms.add(idx)
    except Exception:
        pass

    # 着色: 糖原子=蓝色, 非糖原子=默认
    # Coloring: sugar atoms = blue, non-sugar atoms = default
    highlightAtoms = list(sugarAtoms)
    highlightColors = {}
    for idx in sugarAtoms:
        highlightColors[idx] = (0.3, 0.6, 1.0)  # 蓝色 (Blue)

    # 非糖原子标记为苷元颜色 (橙色)
    # Non-sugar atoms marked as aglycone color (orange)
    allAtoms = set(range(mol.GetNumAtoms()))
    aglyconAtoms = allAtoms - sugarAtoms
    for idx in aglyconAtoms:
        if idx not in highlightColors:
            highlightColors[idx] = (1.0, 0.7, 0.3)  # 橙色 (Orange)

    allHighlight = list(sugarAtoms | aglyconAtoms)
    fullB64 = molToBase64Png(mol, size=(400, 300),
                              highlightAtoms=allHighlight,
                              highlightColors=highlightColors)

    # Glycan 片段 (仅糖原子)
    glycanB64 = ""
    if sugarAtoms:
        glycanB64 = molToBase64Png(mol, size=(350, 250),
                                    highlightAtoms=list(sugarAtoms),
                                    highlightColors={idx: (0.2, 0.5, 1.0) for idx in sugarAtoms})

    # Aglycan 片段 (仅非糖原子)
    aglycanB64 = ""
    if aglyconAtoms and sugarAtoms:  # 只有同时有糖和非糖时才画
        aglycanB64 = molToBase64Png(mol, size=(350, 250),
                                     highlightAtoms=list(aglyconAtoms),
                                     highlightColors={idx: (1.0, 0.6, 0.2) for idx in aglyconAtoms})

    return fullB64, glycanB64, aglycanB64


def generateHtmlReport(results: list, outputPath: str):
    """生成完整 HTML 报告
    Generate comprehensive HTML report with molecular images.
    """
    passCount = sum(1 for r in results if r.get("status") == "PASS" or r.get("sequence"))
    totalCount = len(results)

    html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<title>Glycan Analysis Benchmark Report</title>
<style>
  * {{ margin: 0; padding: 0; box-sizing: border-box; }}
  body {{ font-family: 'Segoe UI', system-ui, sans-serif; background: #0a0e1a; color: #e0e0e0; padding: 20px; }}
  h1 {{ text-align: center; color: #7ec8e3; font-size: 2em; margin: 20px 0; text-shadow: 0 0 20px rgba(126,200,227,0.3); }}
  .stats {{ display: flex; justify-content: center; gap: 30px; margin: 20px 0 30px; flex-wrap: wrap; }}
  .stat-card {{ background: linear-gradient(135deg, #1a2332, #243447); padding: 20px 30px; border-radius: 12px;
    border: 1px solid rgba(126,200,227,0.2); text-align: center; min-width: 150px; }}
  .stat-card .num {{ font-size: 2.5em; font-weight: bold; color: #7ec8e3; }}
  .stat-card .label {{ color: #8899aa; font-size: 0.9em; margin-top: 5px; }}
  .section {{ margin: 30px 0; }}
  .section h2 {{ color: #7ec8e3; border-bottom: 2px solid #2a3a4a; padding-bottom: 8px; margin-bottom: 15px; font-size: 1.4em; }}
  .mol-card {{ background: linear-gradient(135deg, #111827, #1a2332); border: 1px solid #2a3a4a;
    border-radius: 12px; margin: 15px 0; padding: 20px; transition: border-color 0.3s; }}
  .mol-card:hover {{ border-color: #7ec8e3; }}
  .mol-card.pass {{ border-left: 4px solid #4ade80; }}
  .mol-card.fail {{ border-left: 4px solid #f87171; }}
  .mol-header {{ display: flex; justify-content: space-between; align-items: center; margin-bottom: 12px; }}
  .mol-header h3 {{ color: #e2e8f0; font-size: 1.1em; }}
  .badge {{ padding: 4px 12px; border-radius: 20px; font-size: 0.8em; font-weight: 600; }}
  .badge-pass {{ background: rgba(74,222,128,0.15); color: #4ade80; }}
  .badge-fail {{ background: rgba(248,113,113,0.15); color: #f87171; }}
  .badge-info {{ background: rgba(126,200,227,0.15); color: #7ec8e3; }}
  .mol-images {{ display: flex; gap: 10px; margin: 12px 0; flex-wrap: wrap; }}
  .mol-images figure {{ text-align: center; }}
  .mol-images figcaption {{ color: #8899aa; font-size: 0.75em; margin-top: 4px; }}
  .mol-images img {{ border-radius: 8px; border: 1px solid #2a3a4a; background: white; }}
  .mol-info {{ display: grid; grid-template-columns: 1fr 1fr; gap: 8px; font-size: 0.9em; }}
  .mol-info .label {{ color: #8899aa; }}
  .mol-info .value {{ color: #e0e0e0; word-break: break-all; }}
  .smiles {{ font-family: 'Courier New', monospace; font-size: 0.75em; color: #94a3b8;
    background: #0f172a; padding: 6px 10px; border-radius: 6px; overflow-x: auto;
    white-space: nowrap; max-width: 100%; display: block; margin: 8px 0; }}
  .legend {{ display: flex; gap: 20px; justify-content: center; margin: 10px 0; }}
  .legend-item {{ display: flex; align-items: center; gap: 6px; font-size: 0.85em; }}
  .legend-color {{ width: 16px; height: 16px; border-radius: 4px; }}
</style>
</head>
<body>
<h1>🔬 Glycan Analysis Benchmark Report</h1>
<p style="text-align:center;color:#8899aa;">Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}</p>

<div class="stats">
  <div class="stat-card"><div class="num">{totalCount}</div><div class="label">Total Molecules</div></div>
  <div class="stat-card"><div class="num" style="color:#4ade80;">{passCount}</div><div class="label">With Sequence</div></div>
</div>

<div class="legend">
  <div class="legend-item"><div class="legend-color" style="background:rgba(77,153,255,0.8);"></div> Glycan (糖)</div>
  <div class="legend-item"><div class="legend-color" style="background:rgba(255,178,77,0.8);"></div> Aglycone (苷元)</div>
</div>
"""

    # 分组: Original 18 + Benchmark 100
    orig = [r for r in results if r.get("group") == "original_18"]
    bench = [r for r in results if r.get("group") == "benchmark_100"]

    for sectionTitle, items in [("Original 18 Test Molecules", orig), ("Benchmark 100", bench)]:
        html += f'<div class="section"><h2>{sectionTitle} ({len(items)} molecules)</h2>\n'
        for r in items:
            status = "pass" if r.get("sequence") and r["sequence"] != "NO_SUGAR" or r.get("expected") == "NO_SUGAR" else "fail"
            badgeClass = "badge-pass" if status == "pass" else "badge-fail"
            badgeText = "✅" if status == "pass" else "⚠️"

            html += f'''<div class="mol-card {status}">
  <div class="mol-header">
    <h3>#{r.get("id","?")} {r.get("name","Unknown")}</h3>
    <span class="badge {badgeClass}">{badgeText} {r.get("aglycone_type","")}</span>
  </div>
  <code class="smiles">{r.get("smiles","")[:200]}{"..." if len(r.get("smiles","")) > 200 else ""}</code>
  <div class="mol-images">
'''
            # 图片
            if r.get("img_full"):
                html += f'    <figure><img src="data:image/png;base64,{r["img_full"]}" width="400"><figcaption>Full Molecule (Colored)</figcaption></figure>\n'
            if r.get("img_glycan"):
                html += f'    <figure><img src="data:image/png;base64,{r["img_glycan"]}" width="350"><figcaption>Glycan (Blue)</figcaption></figure>\n'
            if r.get("img_aglycan"):
                html += f'    <figure><img src="data:image/png;base64,{r["img_aglycan"]}" width="350"><figcaption>Aglycone (Orange)</figcaption></figure>\n'

            html += '  </div>\n  <div class="mol-info">\n'
            html += f'    <div><span class="label">Sequence:</span></div><div class="value">{r.get("sequence","—")}</div>\n'
            html += f'    <div><span class="label">Expected:</span></div><div class="value">{r.get("expected","—")}</div>\n'
            html += f'    <div><span class="label">Modifications:</span></div><div class="value">{r.get("mods","—")}</div>\n'
            html += f'    <div><span class="label">Func Groups:</span></div><div class="value">{r.get("func_groups","—")}</div>\n'
            if r.get("challenge"):
                html += f'    <div><span class="label">Challenge:</span></div><div class="value">{r.get("challenge","")}</div>\n'
            html += '  </div>\n</div>\n'

        html += '</div>\n'

    html += '</body></html>'

    with open(outputPath, "w", encoding="utf-8") as f:
        f.write(html)


# ============================================================
# 主程序 (Main)
# ============================================================
if __name__ == "__main__":
    print("=" * 80)
    print("  118 分子综合测试 + HTML 报告生成器")
    print("=" * 80)

    allResults = []

    # --- Original 18 ---
    print("\n📋 Running Original 18 Test Molecules...")
    for num, name, smi, exp in ORIGINAL_18:
        print(f"  [{num:>2}/18] {name}...", end=" ", flush=True)
        mol = Chem.MolFromSmiles(smi)
        seq, funcs = analyze_glycan(smi)
        mods = scanGlycanMods(smi) if smi else ""
        if not seq:
            seq = "NO_SUGAR"

        imgFull, imgGlycan, imgAglycan = drawColoredMolecule(mol, smi)

        allResults.append({
            "group": "original_18",
            "id": num,
            "name": name,
            "smiles": smi,
            "sequence": seq,
            "expected": exp,
            "mods": mods if mods else "—",
            "func_groups": funcs if funcs else "—",
            "aglycone_type": "None",
            "challenge": "",
            "img_full": imgFull,
            "img_glycan": imgGlycan,
            "img_aglycan": imgAglycan,
        })
        print(f"→ {seq[:50]}")

    # --- Benchmark 100 ---
    print("\n📋 Running Benchmark 100 Molecules...")
    benchPath = r"d:\Glycan_Database\data\benchmark_100.json"
    with open(benchPath, "r", encoding="utf-8") as f:
        benchData = json.load(f)

    for i, entry in enumerate(benchData):
        smi = entry["smiles"]
        name = entry["name"]
        print(f"  [{i+1:>3}/100] {name[:40]:40s}...", end=" ", flush=True)

        mol = Chem.MolFromSmiles(smi)
        try:
            seq, funcs = analyze_glycan(smi)
        except Exception as e:
            seq, funcs = f"ERROR: {e}", ""

        mods = ""
        try:
            mods = scanGlycanMods(smi) if smi else ""
        except Exception:
            pass

        if not seq:
            seq = "NO_SUGAR"

        imgFull, imgGlycan, imgAglycan = drawColoredMolecule(mol, smi)

        allResults.append({
            "group": "benchmark_100",
            "id": entry["id"],
            "name": name,
            "smiles": smi,
            "sequence": seq,
            "expected": ", ".join(entry.get("expected_sugars", [])),
            "mods": mods if mods else "—",
            "func_groups": funcs if funcs else "—",
            "aglycone_type": entry.get("aglycone_type", ""),
            "challenge": entry.get("challenge", ""),
            "img_full": imgFull,
            "img_glycan": imgGlycan,
            "img_aglycan": imgAglycan,
        })
        print(f"→ {seq[:50] if seq else 'NO_SUGAR'}")

    # --- Generate HTML Report ---
    outPath = r"d:\Glycan_Database\reports\benchmark_report.html"
    os.makedirs(os.path.dirname(outPath), exist_ok=True)
    print(f"\n📝 Generating HTML report...")
    generateHtmlReport(allResults, outPath)
    print(f"\n✅ HTML Report: {outPath}")
    print(f"   Total molecules: {len(allResults)}")
