#!/usr/bin/env python3
"""
==========================================================================
  [EN] GlycoNP Saponin Database — Universal Statistical Charts Generator
       Aggregates parsing logic to generate 40+ Plotly interactive HTMLs 
       and corresponding high-resolution static PNGs.
  [CN] GlycoNP 皂苷子库 — 全局可视化学科研级图谱生成器 
       集成了 40余张 高阶互动HTML 图谱生成的全域脚本。
       
  Input:  reports/GlycoNP_Saponin_DB.csv (Version-agnostic)
  Output: reports/saponin_figures/
==========================================================================
"""
import os
import re
import sys
import json
import warnings
import traceback
from pathlib import Path
from collections import Counter
from typing import List, Dict, Tuple, Optional

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

warnings.filterwarnings("ignore")

# -- Paths ----------------------------------------------------------------
BASE_DIR = Path(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
CSV_PATH = BASE_DIR / "reports" / "GlycoNP_Saponin_DB.csv"
OUT_DIR = BASE_DIR / "reports" / "saponin_figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# -- Global Plotly layout (white background, professional) ----------------
LAYOUT_DEFAULTS = dict(
    template="plotly_white",
    font=dict(family="Arial", size=13, color="#222"),
    paper_bgcolor="white",
    plot_bgcolor="white",
    margin=dict(l=80, r=40, t=80, b=80),
    title_font=dict(size=17, family="Arial", color="#222"),
)

# Viridis-style sequential colorscale for heatmaps (matches reference image)
HEATMAP_COLORSCALE = "Viridis"
# Qualitative palette for categorical data
QUAL_COLORS = px.colors.qualitative.Set2 + px.colors.qualitative.Pastel1

# -- Modification abbreviations → full names (for chart labels) -----------
# These abbreviations come from the Glycan_Modifications column in the DB.
# Adding full names makes charts accessible to non-specialists.
MOD_FULL_NAMES = {
    "O-Ac":     "O-Acetyl",
    "O-Me":     "O-Methyl",
    "NAc":      "N-Acetyl",
    "O-SO3":    "O-Sulfate",
    "Tig":      "Tigloyl (2-MeBut-2-enoyl)",
    "pCou":     "p-Coumaroyl",
    "Fer":      "Feruloyl",
    "Cin":      "Cinnamoyl",
    "Ben":      "Benzoyl",
    "Acyl-10C": "Acyl (C10, decanoyl)",
    "Acyl-30C": "Acyl (C30, long-chain)",
    "Acyl-16C": "Acyl (C16, palmitoyl)",
    "Acyl-18C": "Acyl (C18, stearoyl)",
    "Ang":      "Angeloyl",
    "MeBut":    "2-Methylbutanoyl",
    "O-Pyr":    "O-Pyruvyl",
    "N-Me":     "N-Methyl",
    "dOx":      "Deoxy",
}


def expandModName(tag: str) -> str:
    """Expand a modification abbreviation to its full chemist-readable name."""
    return MOD_FULL_NAMES.get(tag, tag)


# =========================================================================
# ETL LAYER — Load, clean, extract features
# =========================================================================

def extractMonosaccharideList(seq: str) -> List[str]:
    """Extract unique monosaccharide names from a Sugar_Sequence string.
    E.g. 'L-Rha-(b1-2)-D-GlcA' -> ['L-Rha', 'D-GlcA']
    """
    if pd.isna(seq) or not seq:
        return []
    s = str(seq)
    # Remove linkage annotations like (a1-2), (b1-4)
    s = re.sub(r"\([ab]\d-\d\)", " ", s)
    s = s.replace(";", " ").replace("[", " ").replace("]", " ")
    pattern = r"(?<![A-Za-z])([DL]-(?:d|me|6d)?[A-Z][a-z]+[A-Za-z]*|Hex|Pen|dHex|Oct|Hept|HexA)(?![a-z])"
    return re.findall(pattern, s)


def extractModTags(mods: str) -> List[str]:
    """Extract modification tags from Glycan_Modifications.
    E.g. 'L-Qui_1(*O-Ac,*O-Me)' -> ['O-Ac', 'O-Me']
    """
    if pd.isna(mods) or not mods:
        return []
    return re.findall(r"\*([A-Za-z\d\-]+)", str(mods))


def parseSequenceToEdges(seq: str) -> List[Dict]:
    """Parse Sugar_Sequence into structured edges for linkage analysis.
    E.g. 'D-Glc-(b1-4)-L-Rha' -> [{donor:'D-Glc', acceptor:'L-Rha', anomer:'beta', pos:'1->4'}]
    """
    if pd.isna(seq) or not seq:
        return []
    edges = []
    pattern = r"([A-Za-z\d\-]+)-\(([ab])(\d)-(\d)\)-([A-Za-z\d\-]+)"
    for match in re.finditer(pattern, str(seq)):
        donor, anomerChar, posFrom, posTo, acceptor = match.groups()
        anomer = "α" if anomerChar == "a" else "β"
        edges.append({"donor": donor, "acceptor": acceptor,
                      "anomer": anomer, "pos": f"{posFrom}→{posTo}"})
    return edges


def extractMotifs(seq: str, windowSize: int = 2) -> List[str]:
    """Sliding-window motif extraction from sugar sequences."""
    if pd.isna(seq) or not seq:
        return []
    chains = [c.strip() for c in str(seq).split(";")]
    motifs = []
    for chain in chains:
        chain = re.sub(r"^\[.*?\]-", "", chain)
        sugarNames = re.split(r"-\([ab]\d-\d\)-", chain)
        linkages = re.findall(r"\([ab]\d-\d\)", chain)
        if len(sugarNames) < windowSize:
            continue
        for i in range(len(sugarNames) - windowSize + 1):
            fragment = sugarNames[i]
            for j in range(windowSize - 1):
                if i + j < len(linkages):
                    fragment += f"-{linkages[i + j]}-{sugarNames[i + j + 1]}"
            motifs.append(fragment)
    return motifs


def loadAndTransform() -> pd.DataFrame:
    """Load Saponin CSV, filter zero-sugar rows, and compute derived columns."""
    print(f"Loading {CSV_PATH}...")
    df = pd.read_csv(CSV_PATH, low_memory=False)
    totalRaw = len(df)

    # Cast numeric columns
    for col in ["Total_Sugar_Count", "Max_Chain_Length", "exact_molecular_weight",
                "alogp", "topological_polar_surface_area"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Filter: exclude Total_Sugar_Count == 0
    df = df[df["Total_Sugar_Count"] > 0].copy()
    print(f"  Raw: {totalRaw:,} -> Filtered (sugar>0): {len(df):,}")

    # Alias: use Consensus_Sugar_Sequence as Sugar_Sequence for compatibility
    if "Consensus_Sugar_Sequence" in df.columns and "Sugar_Sequence" not in df.columns:
        df["Sugar_Sequence"] = df["Consensus_Sugar_Sequence"]

    # Derived columns
    df["Mono_List"] = df["Sugar_Sequence"].apply(extractMonosaccharideList)
    df["Mod_Tags"] = df["Glycan_Modifications"].apply(extractModTags)
    df["NP_Pathway"] = df["np_classifier_pathway"].fillna("Unknown")
    df["NP_Class"] = df["np_classifier_class"].fillna("Unknown")
    df["NP_Superclass"] = df["np_classifier_superclass"].fillna("Unknown")
    df["Organism_Type"] = df["Organism_Type"].fillna("Unknown")

    # Saponin-specific: classify into Steroidal vs Triterpenoid
    def classifySaponinType(row):
        sc = str(row.get("Super_Scaffold_Class", ""))
        nc = str(row.get("np_classifier_superclass", ""))
        combined = (sc + " " + nc).lower()
        if "steroid" in combined:
            return "Steroidal"
        elif "triterpen" in combined:
            return "Triterpenoid"
        return "Other"
    df["Saponin_Type"] = df.apply(classifySaponinType, axis=1)

    print(f"  ETL complete. {len(df):,} molecules ready.\n")
    return df


def saveHtml(fig, filename: str, height: int = 650, width: int = 1100):
    """Save a Plotly figure as both interactive HTML and high-res PNG.

    PNG export uses kaleido; if not installed, only HTML is saved.
    scale=3 yields ~3300×1950 px images — crisp enough for reports.
    Extra bottom/left margin added in PNG to prevent label clipping.
    """
    curMargin = fig.layout.margin
    customL = getattr(curMargin, "l", None)
    customR = getattr(curMargin, "r", None)
    customT = getattr(curMargin, "t", None)
    customB = getattr(curMargin, "b", None)

    fig.update_layout(**LAYOUT_DEFAULTS)
    
    # Restore custom margins if they were set, otherwise use layout defaults (or fallback values)
    fig.update_layout(margin=dict(
        l=customL if customL is not None else 80,
        r=customR if customR is not None else 40,
        t=customT if customT is not None else 80,
        b=customB if customB is not None else 80,
    ))

    if height:
        fig.update_layout(height=height)
    if width:
        fig.update_layout(width=width)

    # -- HTML --
    outPath = OUT_DIR / filename
    fig.write_html(str(outPath), include_plotlyjs="cdn")
    print(f"  -> Saved: {outPath.name}")

    # -- PNG (high-resolution static export) --
    pngName = filename.replace(".html", ".png")
    pngPath = OUT_DIR / pngName
    try:
        # 保留各图表自定义的 margin，只在默认值更大时才扩展
        # Preserve per-chart custom margins; only expand if defaults are larger
        curMargin = fig.layout.margin
        fig.update_layout(margin=dict(
            l=max(getattr(curMargin, "l", 80) or 80, 100),
            r=max(getattr(curMargin, "r", 40) or 40, 60),
            t=max(getattr(curMargin, "t", 80) or 80, 90),
            b=max(getattr(curMargin, "b", 80) or 80, 100),
        ))
        fig.write_image(str(pngPath), format="png", scale=3,
                        width=width, height=height, engine="kaleido")
        print(f"  -> Saved: {pngName}")
    except Exception as e:
        print(f"  [WARN] PNG export failed ({e}). Install kaleido: pip install -U kaleido")


# =========================================================================
# A1: Monosaccharide Frequency Top-20 Bar Chart
# =========================================================================

def chartA1MonoFrequency(df: pd.DataFrame):
    """Horizontal bar chart of the 20 most frequent monosaccharide types."""
    print("[A1] Monosaccharide Frequency Top-20...")
    allSugars = []
    for lst in df["Mono_List"]:
        if isinstance(lst, list):
            allSugars.extend(lst)
    sc = Counter(allSugars).most_common(20)
    names, counts = zip(*sc) if sc else ([], [])

    fig = go.Figure(go.Bar(
        y=list(reversed(names)), x=list(reversed(counts)),
        orientation="h",
        marker_color=QUAL_COLORS[:len(names)],
        text=list(reversed(counts)), textposition="outside",
    ))
    fig.update_layout(
        title=f"Top 20 Monosaccharide Frequency (N={sum(counts):,} occurrences)",
        xaxis_title="Occurrence Count",
        yaxis_title="Monosaccharide",
        height=600,
    )
    saveHtml(fig, "A1_mono_frequency.html", height=600)


# =========================================================================
# A2: Macro Statistics — Sugar Count + Chain Length + Branching
# =========================================================================

def chartA2MacroStats(df: pd.DataFrame):
    """Three-panel macro statistics: sugar count, chain length, branching."""
    print("[A2] Macro Statistics...")
    fig = make_subplots(rows=1, cols=3,
        subplot_titles=["Sugar Count Distribution", "Max Chain Length", "Linear vs Branched"],
        specs=[[{"type": "histogram"}, {"type": "histogram"}, {"type": "pie"}]])

    # Sugar count histogram
    sugarCounts = df["Total_Sugar_Count"].dropna().astype(int)
    fig.add_trace(go.Histogram(x=sugarCounts, marker_color="#4C78A8",
                               name="Sugar Count", nbinsx=15), row=1, col=1)

    # Chain length histogram
    chainLens = df["Max_Chain_Length"].dropna().astype(int)
    fig.add_trace(go.Histogram(x=chainLens, marker_color="#E45756",
                               name="Chain Length", nbinsx=15), row=1, col=2)

    # Branching pie
    hasBranch = df["Sugar_Sequence"].str.contains(r"\[", na=False)
    branchCounts = hasBranch.value_counts()
    labels = ["Linear", "Branched"]
    values = [branchCounts.get(False, 0), branchCounts.get(True, 0)]
    fig.add_trace(go.Pie(labels=labels, values=values,
                         marker_colors=["#72B7B2", "#B279A2"]), row=1, col=3)

    fig.update_layout(
        title=f"Macro Statistics Overview (N={len(df):,}, 0-sugar excluded)",
        showlegend=False, height=450,
    )
    saveHtml(fig, "A2_macro_statistics.html", height=450)


# =========================================================================
# A3: Monosaccharide × Organism Type Heatmap
# =========================================================================

def chartA3MonoVsOrganism(df: pd.DataFrame):
    """Proportion heatmap: monosaccharide vs organism type (column-normalized)."""
    print("[A3] Monosaccharide vs Organism Type...")
    rows = []
    for _, r in df.iterrows():
        orgType = r.get("Organism_Type", "Unknown")
        if orgType == "Unknown":
            continue
        for sugar in set(r.get("Mono_List", [])):
            rows.append({"Sugar": sugar, "Organism": orgType})
    if not rows:
        print("  [SKIP] No data"); return
    mdf = pd.DataFrame(rows)
    topSugars = mdf["Sugar"].value_counts().head(15).index.tolist()
    topOrgs = [o for o in ["Plant", "Animal", "Bacteria", "Fungi", "Marine"]
               if o in mdf["Organism"].unique()]
    if not topOrgs:
        topOrgs = mdf["Organism"].value_counts().head(5).index.tolist()
    mdf = mdf[mdf["Sugar"].isin(topSugars) & mdf["Organism"].isin(topOrgs)]
    pivot = mdf.groupby(["Sugar", "Organism"]).size().unstack(fill_value=0)

    # Column normalize: % of molecules in each organism containing this sugar
    orgTotals = df[df["Organism_Type"] != "Unknown"].groupby("Organism_Type").size()
    # Add N= sample sizes to column labels
    xLabels = {}
    for col in pivot.columns:
        n = int(orgTotals.get(col, 0))
        if n > 0:
            pivot[col] = pivot[col] / n * 100
        xLabels[col] = f"{col} (N={n:,})"
    pivot = pivot.reindex(index=topSugars).fillna(0)
    pivot.columns = [xLabels.get(c, c) for c in pivot.columns]

    fig = px.imshow(pivot.values, x=pivot.columns.tolist(), y=pivot.index.tolist(),
                    color_continuous_scale=HEATMAP_COLORSCALE, aspect="auto",
                    labels=dict(color="Col %"),
                    text_auto=".1f")
    fig.update_layout(
        title="Monosaccharide × Organism Type<br>"
              "<sub>Column-normalized: within each organism, % of molecules containing this sugar</sub>",
        xaxis_title="Organism Type", yaxis_title="Monosaccharide",
    )
    saveHtml(fig, "A3_mono_vs_organism.html")


# =========================================================================
# A4: Monosaccharide × NP Superclass Heatmap
# =========================================================================

def chartA4MonoVsSuperclass(df: pd.DataFrame):
    """Proportion heatmap: monosaccharide vs NP superclass."""
    print("[A4] Monosaccharide vs NP Superclass...")
    rows = []
    for _, r in df.iterrows():
        sc = r.get("NP_Superclass", "Unknown")
        if sc == "Unknown":
            continue
        for sugar in set(r.get("Mono_List", [])):
            rows.append({"Sugar": sugar, "Superclass": sc})
    if not rows:
        print("  [SKIP] No data"); return
    mdf = pd.DataFrame(rows)
    topSugars = mdf["Sugar"].value_counts().head(12).index.tolist()
    # Only keep saponin-relevant superclasses
    validSc = ["Steroids", "Triterpenoids"]
    mdf = mdf[mdf["Superclass"].isin(validSc)]
    mdf = mdf[mdf["Sugar"].isin(topSugars)]
    pivot = mdf.groupby(["Sugar", "Superclass"]).size().unstack(fill_value=0)
    colTotals = pivot.sum(axis=0)
    pivotNorm = pivot.div(colTotals, axis=1) * 100
    # Add N= to column labels
    pivotNorm.columns = [f"{c} (N={int(colTotals[c]):,})" for c in pivotNorm.columns]

    fig = px.imshow(pivotNorm.values, x=pivotNorm.columns.tolist(),
                    y=pivotNorm.index.tolist(),
                    color_continuous_scale=HEATMAP_COLORSCALE, aspect="auto",
                    labels=dict(color="Col %"), text_auto=".1f")
    fig.update_layout(
        title="Monosaccharide × NP Superclass (Steroids vs Triterpenoids)<br>"
              "<sub>Column-normalized: within each superclass, % of molecules with this sugar</sub>",
        xaxis_title="NP Superclass", yaxis_title="Monosaccharide",
    )
    saveHtml(fig, "A4_mono_vs_superclass.html")


# =========================================================================
# A5: Sugar Sequence × Scaffold Heatmap
# =========================================================================

def chartA5SeqVsScaffold(df: pd.DataFrame):
    """Sugar Sequence × Murcko Scaffold co-occurrence heatmap.

    X 轴按实际计算的 Murcko_Scaffold (纯碳骨架 SMILES) 分组，
    确保分组基于化学结构相似性而非文本标签。
    图片直接从 Murcko_Scaffold SMILES 渲染，保证数据一致性。

    X-axis groups by the actual Murcko_Scaffold column (pure-carbon
    framework SMILES), ensuring classification is based on structural
    similarity rather than text labels. Images are rendered directly
    from Murcko_Scaffold SMILES for data consistency.
    """
    print("[A5] Sugar Sequence vs Murcko Scaffold (pure-carbon)...")
    subDf = df[df["Sugar_Sequence"].notna()].copy()
    subDf = subDf[subDf["Murcko_Scaffold"].notna() & (subDf["Murcko_Scaffold"] != "")]
    if subDf.empty:
        print("  [SKIP] No data"); return

    def truncLabel(s, maxLen=40):
        return s if len(str(s)) <= maxLen else str(s)[:maxLen-3] + "..."

    # 用 Murcko_Scaffold 列直接分组 (Group by Murcko_Scaffold column directly)
    topScaffolds = subDf["Murcko_Scaffold"].value_counts().head(8).index.tolist()
    topSeqs = subDf["Sugar_Sequence"].value_counts().head(15).index.tolist()
    subDf = subDf[subDf["Murcko_Scaffold"].isin(topScaffolds) &
                  subDf["Sugar_Sequence"].isin(topSeqs)]
    subDf["Seq_Label"] = subDf["Sugar_Sequence"].apply(truncLabel)

    scaffoldImages = {}
    scaffoldRingInfo = {}
    # 预计算全量样本数 (Pre-compute full-dataset scaffold counts)
    totalWithScaffold = len(df[df["Murcko_Scaffold"].notna() & (df["Murcko_Scaffold"] != "")])
    fullScaffoldN = df[df["Murcko_Scaffold"].notna()]["Murcko_Scaffold"].value_counts()

    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
        from rdkit.Chem import rdDepictor
        from PIL import ImageDraw, ImageFont
        import io, base64

        try:
            rdDepictor.SetPreferCoordGen(True)
        except Exception:
            pass

        # 尝试加载 Arial 字体 (Try to load Arial font)
        try:
            arialFont = ImageFont.truetype("arial.ttf", 16)
        except Exception:
            arialFont = ImageFont.load_default()

        for smi in topScaffolds:
            mol_scaffold = Chem.MolFromSmiles(smi)
            if mol_scaffold:
                nRings = mol_scaffold.GetRingInfo().NumRings()
                nAtoms = mol_scaffold.GetNumHeavyAtoms()
                scaffoldRingInfo[smi] = (nRings, nAtoms)
            else:
                scaffoldRingInfo[smi] = (0, 0)
                continue

            # 递归剪除 degree-1 原子 (Prune degree-1 dangling atoms)
            rw = Chem.RWMol(mol_scaffold)
            while True:
                to_remove = [a.GetIdx() for a in rw.GetAtoms() if a.GetDegree() == 1]
                if not to_remove:
                    break
                for idx in sorted(to_remove, reverse=True):
                    rw.RemoveAtom(idx)
            try:
                Chem.SanitizeMol(rw)
                drawMol = rw
            except Exception:
                drawMol = mol_scaffold

            try:
                rdDepictor.Compute2DCoords(drawMol)
            except Exception:
                pass

            # 生成骨架图片 (Generate scaffold image)
            img = Draw.MolToImage(drawMol, size=(250, 160))

            # 把样本数 + 占比直接写到图片底部 (Burn N + % into image)
            from PIL import Image
            n = fullScaffoldN.get(smi, 0)
            pct = n / totalWithScaffold * 100 if totalWithScaffold > 0 else 0
            # 拼接一个带文字的新图片 (Create composite image with text footer)
            textLine = f"N={n:,} ({pct:.1f}%)"
            footerH = 28
            composite = Image.new("RGB", (250, 160 + footerH), "white")
            composite.paste(img, (0, 0))
            draw = ImageDraw.Draw(composite)
            # 居中写文字 (Center-align text)
            bbox = draw.textbbox((0, 0), textLine, font=arialFont)
            textW = bbox[2] - bbox[0]
            draw.text(((250 - textW) // 2, 162), textLine, fill="black", font=arialFont)

            buf = io.BytesIO()
            composite.save(buf, format="PNG")
            scaffoldImages[smi] = base64.b64encode(buf.getvalue()).decode()

        print(f"    Rendered {len(scaffoldImages)} purified Murcko scaffold images (with labels).")
    except ImportError:
        print("    [WARN] RDKit not available.")

    # 构建 X 轴标签 (Build simple numeric labels for pivot)
    scaffoldLabels = {}
    for i, smi in enumerate(topScaffolds, 1):
        scaffoldLabels[smi] = f"#{i}"
    subDf["Scaffold_Label"] = subDf["Murcko_Scaffold"].map(scaffoldLabels)

    pivot = subDf.groupby(["Seq_Label", "Scaffold_Label"]).size().unstack(fill_value=0)
    orderedCols = [scaffoldLabels[s] for s in topScaffolds if scaffoldLabels[s] in pivot.columns]
    pivot = pivot.reindex(columns=orderedCols, fill_value=0)

    fig = go.Figure(data=go.Heatmap(
        z=pivot.values,
        x=pivot.columns.tolist(), y=pivot.index.tolist(),
        colorscale=HEATMAP_COLORSCALE,
        text=[[str(int(v)) for v in row] for row in pivot.values],
        texttemplate="%{text}", textfont=dict(size=11, family="Arial"),
        colorbar=dict(title="Count"),
    ))

    # 在 X 轴下方嵌入骨架图片（样本数已烧入图片中）
    if scaffoldImages:
        nCols = len(orderedCols)
        for smi in topScaffolds:
            if smi not in scaffoldImages:
                continue
            label = scaffoldLabels.get(smi, "")
            if label not in orderedCols:
                continue
            colIdx = orderedCols.index(label)
            xFrac = (colIdx + 0.5) / nCols
            fig.add_layout_image(
                source=f"data:image/png;base64,{scaffoldImages[smi]}",
                x=xFrac, y=-0.04,
                xref="paper", yref="paper",
                xanchor="center", yanchor="top",
                sizex=0.9 / nCols, sizey=0.28,
            )

    fig.update_layout(
        title=f"Sugar Sequences × Murcko Scaffold (Top {len(orderedCols)})<br>"
              "<sub>X-axis: Pure-carbon Murcko scaffolds from Aglycon SMILES</sub>",
        xaxis_title="", yaxis_title="Sugar Sequence",
        font=dict(family="Arial"),
        height=1200, width=1400,
        margin=dict(l=120, r=50, t=90, b=420),
    )
    fig.update_xaxes(showticklabels=False)
    saveHtml(fig, "A5_seq_vs_scaffold.html", height=1200, width=1400)




# =========================================================================
# A6: Modification × Organism Heatmap
# =========================================================================

def chartA6ModVsOrganism(df: pd.DataFrame):
    """Modification type vs organism type heatmap."""
    print("[A6] Modification vs Organism...")
    rows = []
    for _, r in df.iterrows():
        orgType = r.get("Organism_Type", "Unknown")
        if orgType == "Unknown":
            continue
        for mod in r.get("Mod_Tags", []):
            rows.append({"Modification": mod, "Organism": orgType})
    if not rows:
        print("  [SKIP] No data"); return
    mdf = pd.DataFrame(rows)
    topMods = mdf["Modification"].value_counts().head(10).index.tolist()
    topOrgs = [o for o in ["Plant", "Animal", "Bacteria", "Fungi", "Marine"]
               if o in mdf["Organism"].unique()]
    mdf = mdf[mdf["Modification"].isin(topMods) & mdf["Organism"].isin(topOrgs)]
    pivot = mdf.groupby(["Modification", "Organism"]).size().unstack(fill_value=0)
    colTotals = pivot.sum(axis=0)
    pivotNorm = pivot.div(colTotals, axis=1) * 100
    # Add N= to column labels
    pivotNorm.columns = [f"{c} (N={int(colTotals[c]):,})" for c in pivotNorm.columns]
    # Expand abbreviations to full names for readability
    pivotNorm.index = [expandModName(m) for m in pivotNorm.index]

    fig = px.imshow(pivotNorm.values, x=pivotNorm.columns.tolist(),
                    y=pivotNorm.index.tolist(),
                    color_continuous_scale=HEATMAP_COLORSCALE, aspect="auto",
                    labels=dict(color="Col %"), text_auto=".1f")
    fig.update_layout(
        title="Modification × Organism Type<br>"
              "<sub>Column-normalized: within each organism, % of modified bonds with this tag</sub>",
        xaxis_title="Organism Type", yaxis_title="Modification",
    )
    saveHtml(fig, "A6_mod_vs_organism.html")


# =========================================================================
# A7: Modification × NP Superclass Heatmap
# =========================================================================

def chartA7ModVsSuperclass(df: pd.DataFrame):
    """Modification type vs NP superclass heatmap."""
    print("[A7] Modification vs NP Superclass...")
    rows = []
    for _, r in df.iterrows():
        sc = r.get("NP_Superclass", "Unknown")
        if sc == "Unknown":
            continue
        for mod in r.get("Mod_Tags", []):
            rows.append({"Modification": mod, "Superclass": sc})
    if not rows:
        print("  [SKIP] No data"); return
    mdf = pd.DataFrame(rows)
    topMods = mdf["Modification"].value_counts().head(10).index.tolist()
    # Only saponin-relevant superclasses
    validSc = ["Steroids", "Triterpenoids"]
    mdf = mdf[mdf["Superclass"].isin(validSc) & mdf["Modification"].isin(topMods)]
    pivot = mdf.groupby(["Modification", "Superclass"]).size().unstack(fill_value=0)
    colTotals = pivot.sum(axis=0)
    pivotNorm = pivot.div(colTotals, axis=1) * 100
    # Add N= to column labels
    pivotNorm.columns = [f"{c} (N={int(colTotals[c]):,})" for c in pivotNorm.columns]
    # Expand abbreviations to full names for readability
    pivotNorm.index = [expandModName(m) for m in pivotNorm.index]

    fig = px.imshow(pivotNorm.values, x=pivotNorm.columns.tolist(),
                    y=pivotNorm.index.tolist(),
                    color_continuous_scale=HEATMAP_COLORSCALE, aspect="auto",
                    labels=dict(color="Col %"), text_auto=".1f")
    fig.update_layout(
        title="Modification × NP Superclass (Steroids vs Triterpenoids)<br>"
              "<sub>Column-normalized: within each superclass, % of modifications of this type</sub>",
        xaxis_title="NP Superclass", yaxis_title="Modification",
    )
    saveHtml(fig, "A7_mod_vs_superclass.html")


# =========================================================================
# A8: Linkage Type Multi-Panel (Organism / Superclass)
# Uses COLUMN-normalization: "within each category, what % of bonds
# are each type?" — more informative than row-normalization for
# datasets dominated by one category (e.g. Plant).
# Also shows raw counts in hover and N= sample sizes on axis labels.
# Pathway chart is dropped because saponins are ~100% Terpenoids.
# =========================================================================

def chartA8LinkageMulti(df: pd.DataFrame):
    """Column-normalized heatmaps of glycosidic bond distribution."""
    print("[A8] Linkage Multi-Panel (column-normalized)...")

    # Clean up stale Pathway chart from previous versions
    staleFile = OUT_DIR / "A8a_linkage_vs_pathway.html"
    stalePng = OUT_DIR / "A8a_linkage_vs_pathway.png"
    for f in [staleFile, stalePng]:
        if f.exists():
            f.unlink()
            print(f"  [CLEANUP] Deleted stale: {f.name}")

    # Parse all bond details to get linkage types
    bondTypes = []
    for _, r in df.iterrows():
        bondRaw = r.get("Glycan-Aglycone_Bond_Detail", "[]")
        try:
            bonds = json.loads(str(bondRaw)) if pd.notna(bondRaw) and str(bondRaw) not in ("nan","","[]") else []
        except Exception:
            bonds = []
        for bd in bonds:
            bond = bd.get("bond", "")
            if bond:
                bondTypes.append({
                    "Bond": bond,
                    "Organism": r.get("Organism_Type", "Unknown"),
                    # For saponins, only Steroids/Triterpenoids are meaningful
                    "Superclass": r.get("NP_Superclass", "Unknown"),
                })
    if not bondTypes:
        print("  [SKIP] No bond data"); return

    bdf = pd.DataFrame(bondTypes)

    # Show ALL bond types found — user asked about N-bond, S-bond
    bondDist = bdf["Bond"].value_counts()
    print(f"    All bond types found (N={len(bdf):,}):")
    for bondType, cnt in bondDist.items():
        print(f"      {bondType:20s}  {cnt:,}")
    hasNbond = any("N" in b for b in bondDist.index)
    hasSbond = any("S" in b for b in bondDist.index)
    print(f"    N-linked bonds present: {hasNbond}")
    print(f"    S-linked bonds present: {hasSbond}")

    # Keep ALL bond types (don't filter) so nothing is hidden
    # Only filter if there are too many (>12) for readability
    if len(bondDist) > 12:
        topBonds = bondDist.head(12).index.tolist()
        bdf = bdf[bdf["Bond"].isin(topBonds)]

    # --- A8b: Organism Type ---
    for dim, dimName, fname in [
        ("Organism", "Organism Type", "A8b_linkage_vs_organism.html"),
        ("Superclass", "NP Superclass", "A8c_linkage_vs_superclass.html"),
    ]:
        sub = bdf[bdf[dim] != "Unknown"]
        # Only keep saponin-relevant superclasses (skip Alkaloids etc.)
        if dim == "Superclass":
            validSc = ["Steroids", "Triterpenoids"]
            sub = sub[sub[dim].isin(validSc)]
        topDim = sub[dim].value_counts().head(6).index.tolist()
        sub = sub[sub[dim].isin(topDim)]
        if sub.empty:
            continue

        # Raw counts pivot
        pivRaw = sub.groupby(["Bond", dim]).size().unstack(fill_value=0)
        # Column totals for N= labels and normalization
        colTotals = pivRaw.sum(axis=0)
        # Column-normalize: "within each column, what % of bonds are this type?"
        pivPct = pivRaw.div(colTotals, axis=1) * 100

        # Build axis labels with N= sample size
        xLabels = [f"{col}\n(N={colTotals[col]:,})" for col in pivPct.columns]
        rowTotals = pivRaw.sum(axis=1)
        yLabels = [f"{idx} (N={rowTotals[idx]:,})" for idx in pivPct.index]

        # Build custom hover text showing both % and raw count
        hoverText = []
        for i, bondType in enumerate(pivPct.index):
            row = []
            for j, category in enumerate(pivPct.columns):
                rawVal = int(pivRaw.iloc[i, j])
                pctVal = pivPct.iloc[i, j]
                row.append(f"{bondType}<br>{category}: {pctVal:.1f}% ({rawVal:,} bonds)")
            hoverText.append(row)

        fig = go.Figure(data=go.Heatmap(
            z=pivPct.values,
            x=xLabels, y=yLabels,
            colorscale=HEATMAP_COLORSCALE,
            text=[[f"{v:.0f}%" for v in row] for row in pivPct.values],
            texttemplate="%{text}", textfont=dict(size=12),
            hovertext=hoverText, hoverinfo="text",
            colorbar=dict(title="Col %"),
        ))
        fig.update_layout(
            title=f"Glycosidic Bond × {dimName} — Column-Normalized<br>"
                  f"<sub>Showing: within each {dimName}, what % of bonds are each type?</sub>",
            xaxis_title=dimName, yaxis_title="Bond Type",
            height=550,
        )
        saveHtml(fig, fname, height=550)


# =========================================================================
# D1: Steroidal vs Triterpenoid Saponin — Sunburst
# =========================================================================

def chartD1SaponinTypeSunburst(df: pd.DataFrame):
    """Sunburst chart: Saponin_Type -> Detailed_NP_Class hierarchy."""
    print("[D1] Saponin Type Sunburst...")
    subDf = df[df["Saponin_Type"].isin(["Steroidal", "Triterpenoid"])].copy()
    subDf["Detail"] = subDf["Detailed_NP_Class"].fillna("Unclassified").apply(
        lambda s: s if len(str(s)) <= 30 else str(s)[:27] + "...")

    fig = px.sunburst(subDf, path=["Saponin_Type", "Detail"],
                      color="Saponin_Type",
                      color_discrete_map={"Steroidal": "#4C78A8", "Triterpenoid": "#E45756"})
    fig.update_layout(title=f"Steroidal vs Triterpenoid Saponin Classification (N={len(subDf):,})")
    saveHtml(fig, "D1_saponin_type_sunburst.html", height=700)


# =========================================================================
# D2: Detailed NP Class Top-15 Bar
# =========================================================================

def chartD2DetailedClassBar(df: pd.DataFrame):
    """Top 15 detailed NP class frequency bar chart."""
    print("[D2] Detailed NP Class Top-15...")
    classCounts = df["Detailed_NP_Class"].dropna().value_counts().head(15)

    fig = go.Figure(go.Bar(
        y=classCounts.index[::-1], x=classCounts.values[::-1],
        orientation="h", marker_color="#4C78A8",
        text=classCounts.values[::-1], textposition="outside",
    ))
    fig.update_layout(
        title=f"Top 15 Detailed NP Class Distribution (N={classCounts.sum():,})",
        xaxis_title="Compound Count", yaxis_title="Detailed NP Class",
        height=550,
    )
    saveHtml(fig, "D2_detailed_class_bar.html", height=550)


# =========================================================================
# D3: Source Family (LOTUS_family) Top-15 Bar
# =========================================================================

def chartD3FamilyBar(df: pd.DataFrame):
    """Top 15 plant/organism families from LOTUS_family."""
    print("[D3] Source Family Top-15...")
    famCounts = df["LOTUS_family"].dropna()
    famCounts = famCounts[famCounts != ""].value_counts().head(15)
    if famCounts.empty:
        print("  [SKIP] No family data"); return

    fig = go.Figure(go.Bar(
        y=famCounts.index[::-1], x=famCounts.values[::-1],
        orientation="h", marker_color="#E45756",
        text=famCounts.values[::-1], textposition="outside",
    ))
    fig.update_layout(
        title=f"Top 15 Source Families (N={famCounts.sum():,})",
        xaxis_title="Compound Count", yaxis_title="Family",
        height=550,
    )
    saveHtml(fig, "D3_source_family_bar.html", height=550)


# =========================================================================
# D4: Molecular Weight vs aLogP Scatter (colored by Saponin Type)
# =========================================================================

def chartD4MwVsAlogp(df: pd.DataFrame):
    """Scatter plot of MW vs aLogP, colored by steroidal/triterpenoid type."""
    print("[D4] MW vs aLogP Scatter...")
    subDf = df[df["exact_molecular_weight"].notna() & df["alogp"].notna()].copy()
    subDf = subDf[subDf["Saponin_Type"].isin(["Steroidal", "Triterpenoid"])]

    fig = px.scatter(subDf, x="exact_molecular_weight", y="alogp",
                     color="Saponin_Type", opacity=0.4,
                     color_discrete_map={"Steroidal": "#4C78A8", "Triterpenoid": "#E45756"},
                     labels={"exact_molecular_weight": "Molecular Weight (Da)",
                             "alogp": "aLogP", "Saponin_Type": "Type"},
                     hover_data=["name", "Detailed_NP_Class"])
    fig.update_layout(
        title=f"Physicochemical Space: MW vs aLogP (N={len(subDf):,})",
        height=600,
    )
    saveHtml(fig, "D4_mw_vs_alogp.html", height=600)


# =========================================================================
# Import Part 2 (A9-A17, B1-B3, C1, D5-D7) and Part 3 (S1-S14) functions
# =========================================================================
import importlib.util as _ilu

def _loadModule(name, path):
    """动态加载模块 (Dynamically load module from file path)。"""
    spec = _ilu.spec_from_file_location(name, path)
    mod = _ilu.module_from_spec(spec)
    # 将当前模块的全局变量注入子模块 (Inject globals into sub-module)
    for key in list(globals().keys()):
        if not key.startswith("_"):
            setattr(mod, key, globals()[key])
    spec.loader.exec_module(mod)
    return mod

_p2path = Path(__file__).parent / "generate_saponin_charts_part2.py"
_p3path = Path(__file__).parent / "generate_saponin_charts_part3.py"

_part2 = _loadModule("saponin_part2", str(_p2path))
_part3 = _loadModule("saponin_part3", str(_p3path))


# =========================================================================
# MAIN: Run ALL chart generators (A1-A17, B1-B3, C1, D1-D7, S1-S14)
# =========================================================================

def main():
    print("=" * 70)
    print("  GlycoNP Saponin Database — Chart Generation (40 charts)")
    print("=" * 70)
    df = loadAndTransform()

    chartFunctions = [
        # ── A 系列: 核心糖链分析 (Core Glycan Analysis) ──
        chartA1MonoFrequency,
        chartA2MacroStats,
        chartA3MonoVsOrganism,
        chartA4MonoVsSuperclass,
        chartA5SeqVsScaffold,
        chartA6ModVsOrganism,
        chartA7ModVsSuperclass,
        chartA8LinkageMulti,
        # ── Part 2: A9-A17 ──
        _part2.chartA9GlycanOrgNetwork,
        _part2.chartA10MotifNetwork,
        _part2.chartA11LinkageCompass,
        _part2.chartA14CoOccurrence,
        _part2.chartA15ScaffoldPCA,
        _part2.chartA16ClassVsSugarCount,
        _part2.chartA17AlphaBetaRatio,
        # ── B 系列: 交叉视图 (Cross-views) ──
        _part2.chartB1Sankey,
        _part2.chartB2SequenceHeatmap,
        _part2.chartB3KingdomMods,
        # ── C 系列: 报告 (Reports) ──
        _part2.reportC1Summary,
        # ── D 系列: 皂苷专题 (Saponin Specifics) ──
        _part2.chartD1SaponinTypeSunburst,
        _part2.chartD2DetailedClassBar,
        _part2.chartD3FamilyBar,
        chartD4MwVsAlogp,
        _part2.chartD5LiteratureCoverage,
        _part2.chartD6Bioactivity,
        _part2.chartD7TpsaVsMw,
        # ── S 系列: 糖链深度分析 Part 3 (Sugar Chain Deep Analysis) ──
        _part3.chartS1DisaccharideFreq,
        _part3.chartS2TrisaccharideFreq,
        _part3.chartS3TetrasaccharideFreq,
        _part3.chartS4TerminalSugar,
        _part3.chartS5ReducingEndSugar,
        _part3.chartS6BranchPointSugar,
        _part3.chartS7MonoBiDesmosidic,
        _part3.chartS8LinkagePositionRadar,
        _part3.chartS9RareSugarDist,
        _part3.chartS10SugarCooccurrence,
        _part3.chartS11PolysaccharideMotifMatch,
        _part3.chartS12SugarDiversityIndex,
        _part3.chartS13ChainLenVsDetailedClass,
        _part3.chartS14SugarFingerprintRadar,
    ]

    successCount = 0
    for func in chartFunctions:
        try:
            func(df)
            successCount += 1
        except Exception as e:
            print(f"  [ERROR] {func.__name__}: {e}")
            traceback.print_exc()

    print(f"\n{'=' * 70}")
    print(f"  Done! {successCount}/{len(chartFunctions)} charts generated.")
    print(f"  Output: {OUT_DIR}")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()


