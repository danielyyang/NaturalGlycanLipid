"""
generate_benchmark_debug_reports.py
===================================
生成 Tier A (50 mol) 和 Benchmark 150 的 HTML 调试报告。
Generate Tier A (50 mol) and Benchmark 150 HTML debug reports.

报告格式与 run_v12_full_pipeline.py 的 step4HtmlReport 一致:
  - 三色高亮分子图 (红=糖环, 黄=修饰, 蓝=苷元)
  - Sugar Sequence (引擎输出 vs 期望值)
  - Bond Detail (糖苷键)
  - Modifications
  - 匹配状态 (PASS/FAIL)

[TEST DATA ONLY]
"""
import os, sys, json, time, re, base64
sys.path.insert(0, r"d:\Glycan_Database")

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

from lib.glycan_topology import find_mapped_sugar_units, get_split_smiles
from lib.monosaccharide_identifier import identify_monosaccharide_v10, generate_refined_sequence

# 复用 run_v12_full_pipeline.py 的糖苷键检测
# Reuse glycosidic bond detection from run_v12_full_pipeline.py
from scripts.run_v12_full_pipeline import (
    molToHighlightedBase64Png, molToBase64Png,
    detectAllGlycosidicBonds, _checkPureSugarMolecule,
)

REPORT_DIR = os.path.join(r"d:\Glycan_Database", "reports")


# =====================================================================
# 评估映射 (Eval Synonym Map) — 从 run_benchmark_200_tier_a.py 复用
# =====================================================================
SYNONYM_MAP = {
    "D-Glucose": "D-Glc", "D-Galactose": "D-Gal", "D-Mannose": "D-Man",
    "D-Xylose": "D-Xyl", "L-Arabinose": "L-Ara", "L-Rhamnose": "L-Rha",
    "L-Fucose": "L-Fuc", "D-Ribose": "D-Rib", "D-Arabinose": "D-Ara",
    "D-Glucuronic acid": "D-GlcA", "D-Galacturonic acid": "D-GalA",
    "N-Acetylneuraminic acid": "Neu5Ac", "D-Glucosamine": "D-GlcN",
    "N-Acetylglucosamine": "D-GlcNAc", "D-Digitoxose": "D-Dtx",
    "L-Oleandrose": "L-Ole", "D-Cymarose": "D-Cym", "L-Thevetose": "L-Tve",
    "D-Mycosamine": "D-Myc", "2-Deoxy-D-ribose": "dRib",
    "L-Cladinose": "L-Cla", "D-Desosamine": "D-Des",
}


def normalizeName(name: str) -> str:
    return SYNONYM_MAP.get(name, name)


# =====================================================================
# HTML 报告生成器 (HTML Report Generator)
# =====================================================================
def generateBenchmarkReport(entries: list, title: str, outputPath: str) -> dict:
    """为基准测试集生成完整 HTML 调试报告。
    Generate full HTML debug report for a benchmark test set.

    Args:
        entries: 测试条目列表 (list of test entries)
        title: 报告标题 (report title)
        outputPath: 输出文件路径 (output file path)

    Returns:
        dict: 聚合统计 (aggregate stats)
    """
    t0 = time.time()

    passCount = 0
    failCount = 0
    errorCount = 0
    totalLatency = 0.0

    html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<title>{title}</title>
<style>
  * {{ margin:0; padding:0; box-sizing:border-box; }}
  body {{ font-family:'Segoe UI',system-ui,sans-serif; background:#0a0e1a; color:#e0e0e0; padding:8px; }}
  h1 {{ text-align:center; color:#7ec8e3; font-size:1.5em; margin:8px 0; }}
  .meta {{ text-align:center; color:#6b7280; font-size:0.75em; margin-bottom:6px; }}
  table {{ width:100%; border-collapse:collapse; font-size:0.72em; table-layout:fixed; }}
  th {{ background:#1a2332; color:#7ec8e3; padding:5px 3px; border:1px solid #2a3a4a;
    position:sticky; top:0; z-index:10; text-align:center; font-size:0.85em; }}
  td {{ padding:3px 2px; border:1px solid #1a2332; vertical-align:top;
    word-break:break-word; background:#111827; overflow:hidden; }}
  tr:hover td {{ background:#1a2332; }}
  tr.pass td {{ border-left:3px solid #22c55e; }}
  tr.fail td {{ border-left:3px solid #ef4444; }}
  img {{ border-radius:2px; border:1px solid #2a3a4a; background:white;
    max-width:100%; height:auto; display:block; }}
  .seq {{ font-family:'Courier New',monospace; font-size:0.85em; color:#a0d0f0; }}
  .expected {{ color:#6b7280; font-size:0.8em; }}
  .pass-tag {{ color:#22c55e; font-weight:bold; }}
  .fail-tag {{ color:#ef4444; font-weight:bold; }}
  .generic {{ color:#f87171; }}
  .bond-tag {{ display:inline-block; padding:1px 4px; border-radius:3px; font-size:0.8em; margin:1px 0; }}
  .bond-O {{ background:#2563eb33; color:#60a5fa; border:1px solid #2563eb; }}
  .bond-S {{ background:#d9770633; color:#fbbf24; border:1px solid #d97706; }}
  .bond-N {{ background:#7c3aed33; color:#c084fc; border:1px solid #7c3aed; }}
  .bond-C {{ background:#05966933; color:#34d399; border:1px solid #059669; }}
  .mod-tag {{ display:inline-block; padding:1px 4px; border-radius:3px; font-size:0.8em;
    margin:1px 0; background:#7c3aed22; color:#c084fc; border:1px solid #7c3aed44; }}
  .no-cleavage {{ color:#f87171; font-style:italic; font-size:0.8em; }}
  col.c1 {{ width:6%; }}  col.c2 {{ width:12%; }} col.c3 {{ width:10%; }}
  col.c4 {{ width:10%; }} col.c5 {{ width:12%; }} col.c6 {{ width:8%; }}
  col.c7 {{ width:6%; }}  col.c8 {{ width:6%; }}  col.c9 {{ width:3%; }}
  col.c10 {{ width:3%; }} col.c11 {{ width:6%; }}
</style>
</head>
<body>
<h1>{title}</h1>
<p class="meta">Entries: {len(entries)} | Generated: {time.strftime('%Y-%m-%d %H:%M')}</p>
"""
    headers = [
        "ID / Name", "Full Molecule", "Glycan Part", "Aglycone Part",
        "Sugar Sequence", "Bond Detail", "Modification", "Expected",
        "Status", "Sugars", "Latency",
    ]
    html += '<table>\n<colgroup>'
    for i in range(1, len(headers) + 1):
        html += f'<col class="c{i}">'
    html += '</colgroup>\n<thead><tr>'
    for h in headers:
        html += f'<th>{h}</th>'
    html += '</tr></thead>\n<tbody>\n'

    for entry in entries:
        entryId = entry.get("id", "?")
        name = entry.get("name", "?")
        smi = entry.get("smiles", "")
        expectedSugars = entry.get("expected_sugars", [])

        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None:
            html += f'<tr class="fail"><td><b>{entryId}</b><br><small>{name}</small></td>'
            html += '<td colspan="10"><span class="fail-tag">SMILES PARSE ERROR</span></td></tr>\n'
            errorCount += 1
            continue

        t1 = time.perf_counter()

        # 运行完整鉴定管线 (Run full identification pipeline)
        try:
            newSeq, newMods = generate_refined_sequence(mol)
        except Exception:
            newSeq, newMods = "", ""
        latencyMs = (time.perf_counter() - t1) * 1000
        totalLatency += latencyMs

        # Glycan/Aglycon split
        try:
            aglyconSmi, glycanSmi = get_split_smiles(mol)
        except Exception:
            aglyconSmi, glycanSmi = "", ""

        # Sugar units + bond detection
        try:
            units = find_mapped_sugar_units(mol)
            _, bondJson = detectAllGlycosidicBonds(mol, units)
        except Exception:
            units = []
            bondJson = "[]"

        # 判定 PASS/FAIL (Determine PASS/FAIL)
        # Bug 3 修复: 核心骨架匹配代替数量匹配
        # Bug 3 fix: Core skeleton matching replaces count-only matching
        gotSugars = []
        for u in units:
            ra = u.get("ring_atoms", [])
            try:
                result = identify_monosaccharide_v10(mol, ra)
                sugarName = result[0] if isinstance(result, tuple) else result
                gotSugars.append(sugarName)
            except Exception:
                gotSugars.append("?")

        def extractCore(name: str) -> str:
            """提取核心骨架名 (Extract core skeleton name).
            'D-Glc(NAc,S)' → 'D-Glc', 'L-Rha(CIP:high)' → 'L-Rha'
            """
            core = re.sub(r'\(.*\)$', '', name)  # 去除括号修饰
            core = re.sub(r'\(CIP:.*\)', '', core)  # 去除 CIP 标签
            return core.strip()

        expNorm = [normalizeName(e) for e in expectedSugars]
        gotCount = len(gotSugars)
        expCount = len(expNorm)

        # 核心骨架匹配: 数量一致 + 每个核心名对应匹配
        # Core match: count consistent + each core name matches
        if gotCount == expCount:
            gotCores = sorted([extractCore(g) for g in gotSugars])
            expCores = sorted([extractCore(e) for e in expNorm])
            matched = gotCores == expCores
        else:
            matched = False
        status = "PASS" if matched else "FAIL"
        if matched:
            passCount += 1
        else:
            failCount += 1

        # 三色高亮图 (Three-color highlighted molecule)
        imgFull = molToHighlightedBase64Png(smi, size=(380, 250))
        imgGlycan = molToBase64Png(glycanSmi, size=(280, 200)) if glycanSmi else ""
        imgAglycon = molToBase64Png(aglyconSmi, size=(280, 200)) if aglyconSmi else ""

        # Sequence display
        seqDisplay = newSeq if newSeq else "—"
        for g in ("Hex", "Pen", "dHex", "HexA", "HexN", "PenN"):
            if g in seqDisplay:
                seqDisplay = seqDisplay.replace(g, f'<span class="generic">{g}</span>', 1)

        # Bond Detail
        bondHtml = ""
        try:
            bonds = json.loads(bondJson) if bondJson else []
            for bd in bonds:
                sugar = bd.get("sugar", "?")
                bond = bd.get("bond", "?")
                bondClass = ""
                if "-O-" in bond: bondClass = "bond-O"
                elif "-S-" in bond: bondClass = "bond-S"
                elif "-N-" in bond: bondClass = "bond-N"
                elif "-C-" in bond: bondClass = "bond-C"
                bondHtml += f'<span class="bond-tag {bondClass}">{sugar}→Aglycon: {bond}</span><br>'
        except Exception:
            pass
        if not bondHtml:
            bondHtml = "—"

        # Modifications
        modHtml = ""
        for u in units:
            for m in u.get("modifications", []):
                modHtml += f'<span class="mod-tag">{m}</span> '
        if not modHtml:
            modHtml = "—"

        # Expected display
        expHtml = ", ".join(expectedSugars) if expectedSugars else "—"

        rowClass = "pass" if status == "PASS" else "fail"
        statusHtml = f'<span class="pass-tag">PASS</span>' if status == "PASS" else f'<span class="fail-tag">FAIL</span>'

        imgFullHtml = f'<img src="data:image/png;base64,{imgFull}">' if imgFull else "—"
        imgGlycanHtml = f'<img src="data:image/png;base64,{imgGlycan}">' if imgGlycan else '<span class="no-cleavage">No Cleavage</span>'
        imgAglyconHtml = f'<img src="data:image/png;base64,{imgAglycon}">' if imgAglycon else '<span class="no-cleavage">No Aglycone</span>'

        html += f'<tr class="{rowClass}">'
        html += f'<td><b>{entryId}</b><br><small>{name[:30]}</small></td>'
        html += f'<td>{imgFullHtml}</td>'
        html += f'<td>{imgGlycanHtml}</td>'
        html += f'<td>{imgAglyconHtml}</td>'
        html += f'<td class="seq">{seqDisplay}<br><span class="expected">Got: {", ".join(gotSugars)}</span></td>'
        html += f'<td>{bondHtml}</td>'
        html += f'<td>{modHtml}</td>'
        html += f'<td class="expected">{expHtml}</td>'
        html += f'<td style="text-align:center">{statusHtml}<br>{gotCount}/{expCount}</td>'
        html += f'<td style="text-align:center">{gotCount}</td>'
        html += f'<td style="text-align:center">{latencyMs:.0f}ms</td>'
        html += '</tr>\n'

    total = passCount + failCount + errorCount
    html += '</tbody>\n</table>\n'
    html += f"""
<div style="text-align:center; margin:12px 0; color:#7ec8e3; font-size:0.9em;">
  <b>Summary:</b> {passCount} PASS / {failCount} FAIL / {errorCount} ERROR
  | Total: {total} | Pass Rate: {passCount/total*100:.1f}%
  | Avg Latency: {totalLatency/total:.1f}ms
  | Generated in {time.time()-t0:.0f}s
</div>
</body>
</html>"""

    with open(outputPath, "w", encoding="utf-8") as f:
        f.write(html)

    stats = {
        "total": total, "pass": passCount, "fail": failCount,
        "error": errorCount, "avg_latency_ms": totalLatency / total if total else 0,
        "pass_rate": passCount / total * 100 if total else 0,
    }
    return stats


def main():
    # === Tier A (50 mol) ===
    print("=" * 70)
    print("  Generating Tier A Debug Report (50 mol)")
    print("=" * 70)
    with open(r"d:\Glycan_Database\data\benchmark_200.json", "r", encoding="utf-8") as f:
        tierA = json.load(f)
    outA = os.path.join(REPORT_DIR, "TierA_50_debug_report.html")
    statsA = generateBenchmarkReport(tierA, "GlycoNP Tier A Debug Report (50 mol)", outA)
    print(f"  Tier A: {statsA['pass']}/{statsA['total']} PASS ({statsA['pass_rate']:.1f}%)")
    print(f"  Avg latency: {statsA['avg_latency_ms']:.1f}ms")
    print(f"  Output: {outA}")

    # === Benchmark 150 ===
    print("\n" + "=" * 70)
    print("  Generating Benchmark 150 Debug Report (147 mol)")
    print("=" * 70)
    with open(r"d:\Glycan_Database\data\benchmark_150.json", "r", encoding="utf-8") as f:
        bench150 = json.load(f)
    out150 = os.path.join(REPORT_DIR, "Benchmark150_debug_report.html")
    stats150 = generateBenchmarkReport(bench150, "GlycoNP Benchmark 150 Debug Report", out150)
    print(f"  Bench150: {stats150['pass']}/{stats150['total']} PASS ({stats150['pass_rate']:.1f}%)")
    print(f"  Avg latency: {stats150['avg_latency_ms']:.1f}ms")
    print(f"  Output: {out150}")

    # === Unified Benchmark (226 mol) ===
    unifiedPath = r"d:\Glycan_Database\data\benchmark_unified.json"
    if os.path.exists(unifiedPath):
        print("\n" + "=" * 70)
        print("  Generating Unified Benchmark Debug Report (226 mol)")
        print("=" * 70)
        with open(unifiedPath, "r", encoding="utf-8") as f:
            unified = json.load(f)
        outU = os.path.join(REPORT_DIR, "Unified_benchmark_debug_report.html")
        statsU = generateBenchmarkReport(unified, "GlycoNP Unified Benchmark Debug Report (226 mol)", outU)
        print(f"  Unified: {statsU['pass']}/{statsU['total']} PASS ({statsU['pass_rate']:.1f}%)")
        print(f"  Avg latency: {statsU['avg_latency_ms']:.1f}ms")
        print(f"  Output: {outU}")

    print("\n" + "=" * 70)
    print("  ALL REPORTS GENERATED!")
    print("=" * 70)


if __name__ == "__main__":
    main()
