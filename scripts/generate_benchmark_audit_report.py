"""
Benchmark 150 可视化审计报告生成器 (Visual Audit Report Generator)
================================================================

为每个分子生成:
  1. 2D 分子结构图 (Base64 内联)
  2. expected vs actual 对比 (标红差异)
  3. 每个糖环的详细鉴定结果
  4. 修饰检测结果

Generates a visual HTML report for each benchmark molecule showing:
  1. 2D molecule image (inline Base64)
  2. Expected vs Actual comparison (red-highlighted differences)
  3. Per-ring identification details
  4. Modification detection results

Usage:
  python scripts/generate_benchmark_audit_report.py
"""
import base64
import io
import json
import os
import re
import sys
import traceback
from typing import List, Tuple, Dict, Optional

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from lib.monosaccharide_identifier import (
    identify_monosaccharide_v10,
    generate_refined_sequence,
)
from lib.glycan_topology import find_mapped_sugar_units


def molToBase64Png(mol: Optional[Chem.Mol], size: Tuple[int, int] = (500, 350)) -> str:
    """将 RDKit Mol 转为 Base64 PNG 字符串。
    Convert RDKit Mol to Base64 PNG string for HTML embedding.
    """
    if mol is None:
        return ""
    try:
        AllChem.Compute2DCoords(mol)
        img = Draw.MolToImage(mol, size=size)
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        return base64.b64encode(buf.getvalue()).decode("utf-8")
    except Exception:
        return ""


def analyzeOneMolecule(entry: Dict) -> Dict:
    """对单个 benchmark 分子运行完整诊断。
    Run full diagnostic on a single benchmark molecule.
    """
    result = {
        "id": entry.get("id", "?"),
        "name": entry.get("name", ""),
        "group": entry.get("group", ""),
        "cat": entry.get("cat", ""),
        "smiles": entry.get("smiles", ""),
        "expected_field": entry.get("expected", ""),
        "expected_sugars": entry.get("expected_sugars", []),
        "expected_sequence": entry.get("expected_sequence", ""),
        "seq_field": entry.get("seq", ""),
        "mods_field": entry.get("mods", ""),
        "actual_sequence": "",
        "actual_mods": "",
        "per_ring": [],
        "status": "ERROR",
        "img_b64": "",
        "issues": [],
    }

    smi = result["smiles"]
    if not smi:
        result["status"] = "NO_SMILES"
        return result

    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        result["status"] = "PARSE_FAIL"
        return result

    result["img_b64"] = molToBase64Png(mol)

    # 检查 expected 字段一致性 (Check expected field consistency)
    # expected 字段 vs expected_sugars 是否逻辑一致
    expectedField = result["expected_field"]
    expectedSugars = result["expected_sugars"]
    if expectedField in ("NO_SUGAR", "N/A", ""):
        result["status"] = "NO_SUGAR"
        return result

    # 引擎运行 (Run engine)
    try:
        seq, mods = generate_refined_sequence(mol)
        result["actual_sequence"] = seq or "EMPTY"
        result["actual_mods"] = mods or ""
    except Exception as e:
        result["actual_sequence"] = f"EXCEPTION: {str(e)[:80]}"
        result["status"] = "ENGINE_ERROR"
        return result

    # 逐环鉴定 (Per-ring identification)
    try:
        units = find_mapped_sugar_units(mol)
        for i, u in enumerate(units):
            ringAtoms = u.get("ring_atoms", [])
            ringSize = len(ringAtoms)
            # 使用 v10 (生产级鉴定引擎, 含 CIP/Exo + Rescue)
            # Using v10 (production-level ID engine with CIP/Exo + Rescue)
            identifyResult = identify_monosaccharide_v10(mol, ringAtoms)
            if isinstance(identifyResult, tuple):
                name, anomer = identifyResult
            else:
                name = str(identifyResult)
                anomer = "?"
            result["per_ring"].append({
                "ring_idx": i,
                "ring_size": ringSize,
                "name": name,
                "anomer": anomer,
            })
    except Exception as e:
        result["issues"].append(f"Ring analysis error: {str(e)[:80]}")

    # 对比评估 (Comparison)
    # 使用 expected_sequence 作为 ground truth (严格比对)
    normActual = re.sub(r'\(CIP:[^)]+\)', '', result["actual_sequence"]).strip()
    normExpected = result["expected_sequence"].strip()

    # 提取 sugar tokens (参考两个字段交叉验证)
    tokenPat = r'Neu5Ac|Neu5Gc|KDO|[DL]-[A-Za-z]+|Hex|dHex|Pen|HexA|HexN|HexNAc'
    actualTokens = re.findall(tokenPat, normActual)

    # 一致性检查: expected_sugars vs expected 字段
    # 用 expected 字段 (由专家标注) 作为权威来源
    expectedTokensFromField = []
    if expectedField:
        # expected 字段格式: "L-IdoA,D-GlcN,D-GlcA,D-GlcN,D-GlcN"
        for part in re.split(r'[,;/\s]+', expectedField):
            part = part.strip()
            if part and part not in ("x3", "x2", ""):
                expectedTokensFromField.append(part)

    # 检查 expected_sugars 与 expected 字段是否一致
    if expectedTokensFromField and expectedSugars:
        # 简单检查: 长度和内容
        if len(expectedTokensFromField) != len(expectedSugars):
            result["issues"].append(
                f"LABEL CONFLICT: expected field has {len(expectedTokensFromField)} "
                f"sugars ({expectedField}) but expected_sugars has "
                f"{len(expectedSugars)} ({expectedSugars})"
            )
        else:
            for i, (ef, es) in enumerate(zip(expectedTokensFromField, expectedSugars)):
                if ef != es and not ef.startswith(es):
                    result["issues"].append(
                        f"LABEL CONFLICT at sugar #{i+1}: "
                        f"expected field='{ef}' vs expected_sugars='{es}'"
                    )

    # 序列比对
    if normActual == normExpected:
        result["status"] = "PASS"
    else:
        result["status"] = "DIFF"
        result["issues"].append(
            f"Sequence mismatch: expected='{normExpected}' vs actual='{normActual}'"
        )

    return result


def generateHtmlReport(results: List[Dict], outputPath: str):
    """生成 HTML 审计报告。
    Generate HTML audit report with molecule images and comparison tables.
    """
    passCount = sum(1 for r in results if r["status"] == "PASS")
    failCount = sum(1 for r in results if r["status"] == "DIFF")
    noSugarCount = sum(1 for r in results if r["status"] == "NO_SUGAR")
    errorCount = sum(1 for r in results if r["status"] in ("ERROR", "ENGINE_ERROR", "PARSE_FAIL"))
    issueCount = sum(1 for r in results if r.get("issues"))

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Benchmark 150 Audit Report</title>
<style>
body {{ font-family: 'Segoe UI', Arial, sans-serif; margin: 20px; background: #f5f5f5; }}
h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
.summary {{ 
    display: flex; gap: 15px; margin: 20px 0; flex-wrap: wrap;
}}
.stat-card {{
    padding: 15px 25px; border-radius: 8px; color: white; font-weight: bold;
    min-width: 150px; text-align: center; font-size: 14px;
}}
.stat-pass {{ background: #27ae60; }}
.stat-fail {{ background: #e74c3c; }}
.stat-issue {{ background: #e67e22; }}
.stat-nosug {{ background: #95a5a6; }}
.stat-card .num {{ font-size: 36px; display: block; }}

.molecule {{
    background: white; border-radius: 10px; margin: 15px 0; padding: 20px;
    box-shadow: 0 2px 8px rgba(0,0,0,0.1); page-break-inside: avoid;
}}
.molecule.status-PASS {{ border-left: 6px solid #27ae60; }}
.molecule.status-DIFF {{ border-left: 6px solid #e74c3c; }}
.molecule.status-NO_SUGAR {{ border-left: 6px solid #95a5a6; }}
.molecule.status-ERROR, .molecule.status-ENGINE_ERROR {{ border-left: 6px solid #8e44ad; }}

.mol-header {{
    display: flex; justify-content: space-between; align-items: center;
    margin-bottom: 10px;
}}
.mol-header h3 {{ margin: 0; color: #2c3e50; }}
.badge {{
    padding: 4px 12px; border-radius: 12px; font-size: 12px;
    font-weight: bold; color: white;
}}
.badge-PASS {{ background: #27ae60; }}
.badge-DIFF {{ background: #e74c3c; }}
.badge-NO_SUGAR {{ background: #95a5a6; }}
.badge-ERROR, .badge-ENGINE_ERROR {{ background: #8e44ad; }}

.mol-body {{ display: flex; gap: 20px; flex-wrap: wrap; }}
.mol-img {{ flex: 0 0 auto; }}
.mol-img img {{ border: 1px solid #ddd; border-radius: 5px; }}
.mol-details {{ flex: 1; min-width: 300px; }}

table.detail {{ width: 100%; border-collapse: collapse; font-size: 13px; }}
table.detail th {{ background: #ecf0f1; text-align: left; padding: 6px 10px; }}
table.detail td {{ padding: 6px 10px; border-bottom: 1px solid #eee; }}
table.detail td.label {{ font-weight: bold; width: 180px; color: #555; }}

.issue-box {{
    background: #ffeaa7; border: 1px solid #fdcb6e; border-radius: 5px;
    padding: 8px 12px; margin-top: 10px; font-size: 12px; color: #856404;
}}
.issue-box strong {{ color: #e74c3c; }}

.ring-table {{ width: 100%; border-collapse: collapse; margin-top: 8px; font-size: 12px; }}
.ring-table th {{ background: #dfe6e9; padding: 4px 8px; text-align: left; }}
.ring-table td {{ padding: 4px 8px; border-bottom: 1px solid #eee; }}

.seq-compare {{ font-family: 'Consolas', monospace; font-size: 13px; }}
.seq-match {{ color: #27ae60; }}
.seq-diff {{ color: #e74c3c; font-weight: bold; background: #ffe0e0; padding: 1px 3px; }}

.filter-bar {{ margin: 15px 0; }}
.filter-bar button {{
    padding: 8px 16px; margin-right: 8px; border: 1px solid #bdc3c7;
    border-radius: 5px; cursor: pointer; background: white; font-size: 13px;
}}
.filter-bar button.active {{ background: #3498db; color: white; border-color: #3498db; }}
.filter-bar button:hover {{ background: #ecf0f1; }}
</style>
</head>
<body>
<h1>Benchmark 150 Audit Report</h1>
<p style="color:#666;">Generated for expert review. Each molecule shows: structure image, 
expected vs actual sugar identification, per-ring details, and detected issues.</p>

<div class="summary">
    <div class="stat-card stat-pass"><span class="num">{passCount}</span>PASS</div>
    <div class="stat-card stat-fail"><span class="num">{failCount}</span>DIFF</div>
    <div class="stat-card stat-issue"><span class="num">{issueCount}</span>HAS ISSUES</div>
    <div class="stat-card stat-nosug"><span class="num">{noSugarCount}</span>NO SUGAR</div>
</div>

<div class="filter-bar">
    <button class="active" onclick="filterMols('all')">All ({len(results)})</button>
    <button onclick="filterMols('DIFF')">Failures ({failCount})</button>
    <button onclick="filterMols('issue')">Has Issues ({issueCount})</button>
    <button onclick="filterMols('PASS')">Pass ({passCount})</button>
</div>
<script>
function filterMols(type) {{
    document.querySelectorAll('.filter-bar button').forEach(b => b.classList.remove('active'));
    event.target.classList.add('active');
    document.querySelectorAll('.molecule').forEach(el => {{
        if (type === 'all') {{ el.style.display = ''; return; }}
        if (type === 'issue') {{
            el.style.display = el.dataset.hasIssue === 'true' ? '' : 'none';
            return;
        }}
        el.style.display = el.dataset.status === type ? '' : 'none';
    }});
}}
</script>
"""

    for r in results:
        status = r["status"]
        hasIssue = "true" if r.get("issues") else "false"
        html += f'<div class="molecule status-{status}" data-status="{status}" data-has-issue="{hasIssue}">\n'

        # Header
        html += f'<div class="mol-header">'
        html += f'<h3>#{r["id"]} — {r["name"]} <small style="color:#999">({r["cat"]}, {r["group"]})</small></h3>'
        html += f'<span class="badge badge-{status}">{status}</span>'
        html += '</div>\n'

        if status == "NO_SUGAR":
            html += '<p style="color:#999">No sugar expected / detected.</p></div>\n'
            continue

        # Body
        html += '<div class="mol-body">\n'

        # Image
        if r["img_b64"]:
            html += f'<div class="mol-img"><img src="data:image/png;base64,{r["img_b64"]}" width="450" height="320"></div>\n'

        # Details table
        html += '<div class="mol-details">\n'
        html += '<table class="detail">\n'
        html += f'<tr><td class="label">Expected (expert field)</td><td><b>{r["expected_field"]}</b></td></tr>\n'
        html += f'<tr><td class="label">expected_sugars (JSON)</td><td>{", ".join(r["expected_sugars"])}</td></tr>\n'
        html += f'<tr><td class="label">expected_sequence (JSON)</td><td class="seq-compare">{r["expected_sequence"]}</td></tr>\n'
        html += f'<tr><td class="label">Actual sequence (engine)</td><td class="seq-compare"><b>{r["actual_sequence"]}</b></td></tr>\n'
        html += f'<tr><td class="label">Actual mods (engine)</td><td>{r["actual_mods"]}</td></tr>\n'
        html += f'<tr><td class="label">Benchmark mods</td><td>{r["mods_field"]}</td></tr>\n'
        html += '</table>\n'

        # Per-ring table
        if r["per_ring"]:
            html += '<table class="ring-table">\n'
            html += '<tr><th>Ring #</th><th>Size</th><th>Identified As</th><th>Anomer</th></tr>\n'
            for pr in r["per_ring"]:
                html += f'<tr><td>{pr["ring_idx"]}</td><td>{pr["ring_size"]}</td>'
                html += f'<td><b>{pr["name"]}</b></td><td>{pr["anomer"]}</td></tr>\n'
            html += '</table>\n'

        # Issues
        if r.get("issues"):
            html += '<div class="issue-box">\n'
            for iss in r["issues"]:
                html += f'<strong>Issue:</strong> {iss}<br>\n'
            html += '</div>\n'

        html += '</div>\n'  # mol-details
        html += '</div>\n'  # mol-body
        html += '</div>\n'  # molecule

    html += """
<footer style="margin-top:30px; padding:15px; color:#999; font-size:12px; border-top:1px solid #eee;">
GlycoNP-Pipeline Benchmark Audit Report | Generated for expert review
</footer>
</body></html>"""

    with open(outputPath, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"Report saved: {outputPath}")


def main():
    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    benchPath = os.path.join(baseDir, "data", "benchmark_150.json")
    outputPath = os.path.join(baseDir, "reports", "benchmark_150_test_report.html")

    with open(benchPath, "r", encoding="utf-8") as f:
        benchmarks = json.load(f)

    print(f"Loaded {len(benchmarks)} benchmark entries")
    print(f"Running diagnostics...")

    results = []
    from tqdm import tqdm
    for entry in tqdm(benchmarks, desc="Analyzing", ncols=80):
        try:
            r = analyzeOneMolecule(entry)
        except Exception as e:
            r = {
                "id": entry.get("id", "?"),
                "name": entry.get("name", ""),
                "group": entry.get("group", ""),
                "cat": entry.get("cat", ""),
                "smiles": entry.get("smiles", ""),
                "expected_field": entry.get("expected", ""),
                "expected_sugars": entry.get("expected_sugars", []),
                "expected_sequence": entry.get("expected_sequence", ""),
                "seq_field": entry.get("seq", ""),
                "mods_field": entry.get("mods", ""),
                "actual_sequence": f"CRASH: {str(e)[:60]}",
                "actual_mods": "",
                "per_ring": [],
                "status": "ERROR",
                "img_b64": "",
                "issues": [f"Exception: {traceback.format_exc()[:200]}"],
            }
        results.append(r)

    generateHtmlReport(results, outputPath)

    # Summary
    statusCounts = {}
    issueMols = []
    for r in results:
        statusCounts[r["status"]] = statusCounts.get(r["status"], 0) + 1
        if r.get("issues"):
            issueMols.append((r["id"], r["name"], r["issues"]))

    print(f"\nStatus summary:")
    for s, c in sorted(statusCounts.items()):
        print(f"  {s}: {c}")

    if issueMols:
        print(f"\nMolecules with issues ({len(issueMols)}):")
        for mid, mname, issues in issueMols:
            print(f"  #{mid} {mname}:")
            for iss in issues:
                print(f"    - {iss[:100]}")


if __name__ == "__main__":
    main()
