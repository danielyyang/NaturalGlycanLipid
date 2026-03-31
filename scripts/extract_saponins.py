"""
extract_saponins.py
===================
[EN] GlycoNP Saponin Exporter
     Domain-specific script to extract and visualize Saponin-class molecules 
     from the global natural product dataset. It builds a user-friendly HTML 
     report containing structural badges, interactive 2D renderings, and literature DOIs.

[CN] GlycoNP 皂苷专属提取器与可视化生成管线
     基于在全量天然产物数据库中的挖掘，本脚本专业提取“皂苷(Saponin)”门类（包括
     三萜类与甾体类皂苷）。通过生成带有直观结构徽章 (Badges)、糖断键图以及 2D
     RDKit 渲染图的 HTML 专家级交互报告，为化学工作者提供复核依据。

[TEST DATA ONLY] - GlycoNP Project
"""

import pandas as pd
import os, sys, time, json, re, base64
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from tqdm import tqdm

from lib.pipeline_utils import (
    molToHighlightedBase64Png, molToBase64Png,
    detectAllGlycosidicBonds
)
from lib.glycan_topology import find_mapped_sugar_units

REPORT_DIR = os.path.join(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")), "reports")
SAPONIN_CSV = os.path.join(REPORT_DIR, "GlycoNP_Saponin_DB.csv")


def buildSaponinHtmlReport(df: pd.DataFrame, max_rows: int = 1000) -> str:
    """Generate the GlycoNP-style HTML debug report for Saponin subset.
    Matches the column layout: ID/Name | Full Mol | Glycan | Aglycone |
    Sugar Seq | Bond Detail | Mod | Organism | Classification | Sugars | Chain | DOI | Bioactivity.
    """
    t0 = time.time()

    sampleDf = df.head(max_rows).reset_index(drop=True)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>GlycoNP Saponin Report</title>
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
  img {{ border-radius:2px; border:1px solid #2a3a4a; background:white;
    max-width:100%; height:auto; display:block; margin:auto; }}
  .seq {{ font-family:'Courier New',monospace; font-size:0.8em; color:#a0d0f0; }}
  .generic {{ color:#f87171; }}
  .bond-tag {{ display:inline-block; padding:1px 4px; border-radius:3px; font-size:0.8em;
    margin:1px 0; }}
  .bond-O {{ background:#2563eb33; color:#60a5fa; border:1px solid #2563eb; }}
  .bond-S {{ background:#d9770633; color:#fbbf24; border:1px solid #d97706; }}
  .bond-N {{ background:#7c3aed33; color:#c084fc; border:1px solid #7c3aed; }}
  .bond-C {{ background:#05966933; color:#34d399; border:1px solid #059669; }}
  .mod-tag {{ display:inline-block; padding:1px 4px; border-radius:3px; font-size:0.8em;
    margin:1px 0; background:#7c3aed22; color:#c084fc; border:1px solid #7c3aed44; }}
  .org-Plant   {{ color:#4ade80; font-weight:bold; }}
  .org-Fungi   {{ color:#fb923c; font-weight:bold; }}
  .org-Bacteria {{ color:#60a5fa; font-weight:bold; }}
  .org-Marine  {{ color:#38bdf8; font-weight:bold; }}
  .org-Unknown {{ color:#6b7280; font-weight:bold; }}
  .no-cleavage {{ color:#f87171; font-style:italic; font-size:0.8em; display:block; text-align:center; padding:8px; }}
  .doi {{ font-size:0.75em; color:#93c5fd; word-break:break-all; }}
  .bio {{ font-size:0.75em; color:#86efac; }}
  col.c1 {{ width:6%; }}  col.c2 {{ width:11%; }} col.c3 {{ width:9%; }}
  col.c4 {{ width:9%; }}  col.c5 {{ width:10%; }} col.c6 {{ width:5%; }}
  col.c7 {{ width:5%; }}  col.c8 {{ width:7%; }}  col.c9 {{ width:6%; }}
  col.c10 {{ width:3%; }} col.c11 {{ width:3%; }} col.c12 {{ width:9%; }}
  col.c13 {{ width:10%; }}
</style>
</head>
<body>
<h1>GlycoNP Saponin Report</h1>
<p class="meta">Total Saponins in Database: {len(df):,} | Showing first {min(max_rows, len(df)):,} rows | Generated: {time.strftime('%Y-%m-%d %H:%M')}</p>
<table>
<colgroup>
{"".join(f'<col class="c{i}">' for i in range(1,14))}
</colgroup>
<thead><tr>
  <th>ID / Name</th><th>Full Molecule</th><th>Glycan Part</th><th>Aglycone Part</th>
  <th>Sugar Sequence</th><th>Bond Detail</th><th>Modification</th><th>Organism</th>
  <th>Classification</th><th>Sugars</th><th>Chain</th><th>DOI / References</th>
  <th>Bioactivity</th>
</tr></thead>
<tbody>
"""

    for _, row in tqdm(sampleDf.iterrows(), total=len(sampleDf), desc="Rendering HTML"):
        # ── SMILES ──────────────────────────────────────────────────────────────
        fullSmi    = str(row.get("canonical_smiles", "")) if pd.notna(row.get("canonical_smiles")) else ""
        glycanSmi  = str(row.get("Glycan_SMILES",   "")) if pd.notna(row.get("Glycan_SMILES"))   else ""
        aglyconSmi = str(row.get("Aglycon_SMILES",  "")) if pd.notna(row.get("Aglycon_SMILES"))  else ""

        # ── Molecular images — FULL molecule gets three-color highlights ─────────
        imgFull    = molToHighlightedBase64Png(fullSmi, size=(380, 250))
        imgGlycan  = molToBase64Png(glycanSmi,  size=(280, 200))
        imgAglycon = molToBase64Png(aglyconSmi, size=(280, 200))

        imgFullHtml    = f'<img src="data:image/png;base64,{imgFull}">'    if imgFull    else "—"
        imgGlycanHtml  = (f'<img src="data:image/png;base64,{imgGlycan}">'  if imgGlycan
                          else ('<span class="no-cleavage">No Cleavage</span>'  if not glycanSmi  else '<span class="no-cleavage">Parse Error</span>'))
        imgAglyconHtml = (f'<img src="data:image/png;base64,{imgAglycon}">' if imgAglycon
                          else ('<span class="no-cleavage">No Aglycone</span>' if not aglyconSmi else '<span class="no-cleavage">Parse Error</span>'))

        # ── Sugar Sequence ────────────────────────────────────────────────────
        # v13: column is Consensus_Sugar_Sequence (LLM-validated)
        seq = str(row.get("Consensus_Sugar_Sequence", "")) if pd.notna(row.get("Consensus_Sugar_Sequence")) else ""
        seqHtml = seq
        for g in ("Hex", "Pen", "dHex", "HexA", "Hept", "Oct"):
            if g in seqHtml:
                seqHtml = seqHtml.replace(g, f'<span class="generic">{g}</span>', 1)
        seqHtml = seqHtml or "—"

        # ── Bond Detail — parse Glycan-Aglycone_Bond_Detail JSON ────────────
        # v13: column is "Glycan-Aglycone_Bond_Detail"
        bondDetailRaw = str(row.get("Glycan-Aglycone_Bond_Detail", "[]")) if pd.notna(row.get("Glycan-Aglycone_Bond_Detail")) else "[]"
        bondHtml = ""
        try:
            bonds = json.loads(bondDetailRaw) if bondDetailRaw and bondDetailRaw not in ("nan", "", "[]") else []
            seenSugars: set = set()
            for bd in bonds:
                target = bd.get("target", "")
                if "Aglycon" not in target and "aglycon" not in target.lower():
                    continue
                sugar = bd.get("sugar", "?")
                if sugar in seenSugars:
                    continue
                seenSugars.add(sugar)
                bond = bd.get("bond", "?")
                bondClass = ("bond-O" if "-O-" in bond else
                             "bond-S" if "-S-" in bond else
                             "bond-N" if "-N-" in bond else
                             "bond-C" if "-C-" in bond else "")
                bondHtml += f'<span class="bond-tag {bondClass}">{sugar}→Aglycon: {bond}</span><br>'
        except Exception:
            pass
        if not bondHtml:
            # fallback: use Aglycone_Linkage_Type if present
            rootBond = str(row.get("Aglycone_Linkage_Type", "")) if pd.notna(row.get("Aglycone_Linkage_Type")) else ""
            if rootBond and rootBond != "nan":
                bondClass = ("bond-O" if "-O-" in rootBond else
                             "bond-S" if "-S-" in rootBond else
                             "bond-N" if "-N-" in rootBond else "")
                bondHtml = f'<span class="bond-tag {bondClass}">{rootBond}</span>'
            else:
                bondHtml = "—"

        # ── Modification ──────────────────────────────────────────────────────
        modHtml = ""
        glycanMods = str(row.get("Glycan_Modifications", "")) if pd.notna(row.get("Glycan_Modifications")) else ""
        if glycanMods and glycanMods not in ("nan", "None", "", "NULL"):
            modTokens = re.findall(r'\*([A-Za-z0-9_-]+)', glycanMods)
            seen: set = set()
            for mt in modTokens:
                if mt not in seen:
                    modHtml += f'<span class="mod-tag">{mt}</span> '
                    seen.add(mt)
        if not modHtml:
            modHtml = "—"

        # ── Organism ─────────────────────────────────────────────────────────
        orgType  = str(row.get("Organism_Type", "Unknown")) if pd.notna(row.get("Organism_Type")) else "Unknown"
        organism = (str(row.get("organisms", ""))[:35] if pd.notna(row.get("organisms")) else "")
        family   = (str(row.get("LOTUS_family", ""))[:20] if pd.notna(row.get("LOTUS_family")) else "")
        orgHtml = f'<span class="org-{orgType}"><b>{orgType}</b></span><br>{organism}<br><small>{family}</small>'

        # ── Classification ────────────────────────────────────────────────────
        npSuper      = str(row.get("np_classifier_superclass", "")) if pd.notna(row.get("np_classifier_superclass")) else ""
        npClass      = str(row.get("np_classifier_class", ""))      if pd.notna(row.get("np_classifier_class"))      else ""
        scaffoldCls  = str(row.get("Super_Scaffold_Class", ""))     if pd.notna(row.get("Super_Scaffold_Class"))     else ""
        detailedCls  = str(row.get("Detailed_NP_Class", ""))        if pd.notna(row.get("Detailed_NP_Class"))        else ""
        displayClass = npSuper or scaffoldCls or "Saponins"
        subClass     = detailedCls or npClass
        classHtml = f'<b>{displayClass[:22]}</b><br><small>{subClass[:22]}</small>'

        # ── Counts ───────────────────────────────────────────────────────────
        totalSugar = str(int(row.get("Total_Sugar_Count", 0))) if pd.notna(row.get("Total_Sugar_Count")) else "—"
        maxChain   = str(int(row.get("Max_Chain_Length",  0))) if pd.notna(row.get("Max_Chain_Length"))  else "—"

        # ── DOI / References — real `dois` column (pipe-separated) ────────────
        doisRaw = str(row.get("dois", "")) if pd.notna(row.get("dois")) else ""
        if doisRaw and doisRaw not in ("nan", "None", ""):
            # Split pipe-separated DOIs, truncate list for display
            doi_list = [d.strip() for d in doisRaw.split("|") if d.strip()]
            # Show first 3 as clickable links, then a count badge
            doi_links = []
            for d in doi_list[:3]:
                doi_links.append(f'<a href="https://doi.org/{d}" class="doi" target="_blank">{d[:45]}</a>')
            if len(doi_list) > 3:
                doi_links.append(f'<span style="color:#6b7280;font-size:0.75em">+{len(doi_list)-3} more</span>')
            doiHtml = "<br>".join(doi_links)
        else:
            # Fallback: show PubChem CID if no DOI
            pubchemCid = str(row.get("PubChem_CID", "")) if pd.notna(row.get("PubChem_CID")) else ""
            doiHtml = f'<span class="doi">PubChem: {pubchemCid}</span>' if pubchemCid and pubchemCid != "nan" else "—"

        # ── Bioactivity ───────────────────────────────────────────────────────
        bioRaw      = str(row.get("bioactivity_summary", "")) if pd.notna(row.get("bioactivity_summary")) else ""
        chemblTgts  = str(row.get("ChEMBL_Targets",      "")) if pd.notna(row.get("ChEMBL_Targets"))      else ""
        bioDisplay  = bioRaw
        if chemblTgts and chemblTgts not in ("nan", ""):
            bioDisplay += (" " if bioDisplay else "") + f"[{chemblTgts[:40]}]"
        bioHtml = f'<span class="bio">{bioDisplay[:80]}</span>' if bioDisplay else "—"

        # ── ID / Name ─────────────────────────────────────────────────────────
        identifier = str(row.get("identifier", ""))[:15] if pd.notna(row.get("identifier")) else ""
        nameTxt    = str(row.get("name",       ""))[:35] if pd.notna(row.get("name"))        else ""
        nameHtml   = f"<b>{identifier}</b><br><small>{nameTxt}</small>"

        # ── Row ───────────────────────────────────────────────────────────────
        html += '<tr>'
        html += f'<td>{nameHtml}</td>'
        html += f'<td>{imgFullHtml}</td>'
        html += f'<td>{imgGlycanHtml}</td>'
        html += f'<td>{imgAglyconHtml}</td>'
        html += f'<td class="seq">{seqHtml}</td>'
        html += f'<td>{bondHtml}</td>'
        html += f'<td>{modHtml}</td>'
        html += f'<td>{orgHtml}</td>'
        html += f'<td><small>{classHtml}</small></td>'
        html += f'<td style="text-align:center">{totalSugar}</td>'
        html += f'<td style="text-align:center">{maxChain}</td>'
        html += f'<td>{doiHtml}</td>'
        html += f'<td>{bioHtml}</td>'
        html += '</tr>\n'

    html += '</tbody>\n</table>\n</body>\n</html>'
    print(f"\nHTML rendered in {time.time()-t0:.1f}s")
    return html


if __name__ == "__main__":
    print("Loading Full Unified dataset...")
    dfPruned = pd.read_csv(os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched_Final.csv"), low_memory=False)

    # Filter Saponins across all three classification columns
    mask = (
        dfPruned['np_classifier_class'].str.contains('saponin', case=False, na=False) |
        dfPruned['Detailed_NP_Class'].str.contains('saponin', case=False, na=False) |
        dfPruned['Super_Scaffold_Class'].str.contains('saponin', case=False, na=False)
    )
    saponinDf = dfPruned[mask].copy()

    print(f"Identified {len(saponinDf):,} Saponin molecules.")

    # Save dedicated sub-database CSV
    saponinDf.to_csv(SAPONIN_CSV, index=False)
    print(f"Saved: {SAPONIN_CSV}")

    # Generate HTML report (top 1000 to keep file size manageable)
    print("Generating HTML Visual Report (top 1000 rows)...")
    htmlContent = buildSaponinHtmlReport(saponinDf, max_rows=1000)

    outPath = os.path.join(REPORT_DIR, "GlycoNP_Saponin_Visual_Report.html")
    with open(outPath, "w", encoding="utf-8") as f:
        f.write(htmlContent)
    print(f"Saved HTML Report: {outPath}")
