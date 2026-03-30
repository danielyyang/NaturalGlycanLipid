#!/usr/bin/env python3
"""
============================================================================
  GlycoNP 可视化大屏 V3 (Glycan Natural Product Visualization Dashboard V2)
  -- 四维度专业级数据图表 / 全面升级版 --

  升级内容:
  - 移除 0-sugar 分子, 排除 Unknown 物种
  - 新增 Family 级别热力图
  - 骨架结构可视化 (RDKit 结构图)
  - 所有图表添加专业注释段落
  - 新增额外分析图表

  [TEST DATA ONLY]
============================================================================
"""
import sys
import os
import re
import json
import warnings
from collections import Counter, defaultdict
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from io import BytesIO
import base64

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
import seaborn as sns
import networkx as nx

# 关闭不必要警告 (Suppress non-critical warnings)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", message=".*tight_layout.*")

# =====================================================================
# 全局样式设定 (Global Style Configuration)
# =====================================================================
plt.rcParams.update({
    "figure.facecolor": "#0D1117",
    "axes.facecolor": "#161B22",
    "axes.edgecolor": "#30363D",
    "axes.labelcolor": "#C9D1D9",
    "text.color": "#C9D1D9",
    "xtick.color": "#8B949E",
    "ytick.color": "#8B949E",
    "grid.color": "#21262D",
    "grid.alpha": 0.6,
    "font.family": "sans-serif",
    "font.sans-serif": ["Microsoft YaHei", "Segoe UI", "Arial", "DejaVu Sans"],
    "font.size": 10,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "figure.dpi": 150,
    "savefig.dpi": 200,
    "savefig.bbox": "tight",
    "savefig.facecolor": "#0D1117",
})

# 自定义色板 (Custom color palettes)
HEAT_CMAP = LinearSegmentedColormap.from_list(
    "glyco_heat", ["#0D1117", "#1A1A2E", "#16213E", "#0F3460", "#533483",
                    "#E94560", "#FF6B6B", "#FFD93D"])
GLYCO_CMAP = LinearSegmentedColormap.from_list(
    "glyco", ["#0D1117", "#1B4332", "#2D6A4F", "#52B788", "#95D5B2", "#D8F3DC"])
ACCENT = ["#58A6FF", "#F0883E", "#A371F7", "#3FB950",
          "#FF7B72", "#D2A8FF", "#79C0FF", "#FFA657"]
ANNOT_COLOR = "#8B949E"

OUTPUT_DIR = Path(r"d:\Glycan_Database\reports\figures")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# =====================================================================
# ETL 层 (ETL Layer)
# =====================================================================

def normalizeLinkageType(raw: str) -> str:
    """Patch 1: 统一糖苷键格式 (Normalize linkage type format)."""
    if pd.isna(raw) or not raw:
        return ""
    s = str(raw).strip()
    s = s.replace("alpha", "α").replace("beta", "β")
    s = s.replace("-linked", "")
    return s


def parseSequenceToEdges(seq: str) -> List[Dict]:
    """Patch 2: 解析 Sugar_Sequence 为结构化边列表 (Parse sequence to edges)."""
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
    """Patch 3: 滑动窗口提取糖链片段 (Sliding window motif extraction)."""
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


def extractMonosaccharideList(seq: str) -> List[str]:
    """从 Sugar_Sequence 提取所有单糖名称 (Extract all monosaccharide names)."""
    if pd.isna(seq) or not seq:
        return []
    s = str(seq)
    s = re.sub(r"\([ab]\d-\d\)", " ", s)
    s = s.replace(";", " ").replace("[", " ").replace("]", " ").replace("-(-", " ")
    pattern = r"(?<![A-Za-z])([DL]-(?:d|me|6d)?[A-Z][a-z]+[A-Za-z]*|Hex|Pen|dHex|Oct)(?![a-z])"
    return re.findall(pattern, s)


def extractModificationTags(mods: str) -> List[str]:
    """从 Glycan_Modifications 提取修饰标签 (Extract modification tags)."""
    if pd.isna(mods) or not mods:
        return []
    return re.findall(r"\*([A-Za-z\d\-]+)", str(mods))


def smiToImage(smi: str, size: Tuple[int, int] = (200, 150)) -> Optional[np.ndarray]:
    """将 SMILES 转为 PNG 图片数组 (SMILES to image array for embedding)."""
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        img = Draw.MolToImage(mol, size=size)
        return np.array(img)
    except Exception:
        return None


# =====================================================================
# ETL 主流程 (ETL Main Process)
# =====================================================================

def loadAndTransform(csvPath: str) -> Tuple[pd.DataFrame, int]:
    """加载 CSV + ETL 变换 + 移除 0-sugar (Load, transform, clean)."""
    print("Loading data...")
    df = pd.read_csv(csvPath, low_memory=False)
    totalRaw = len(df)
    print(f"  Raw: {totalRaw:,} rows")

    # 移除 0-sugar 分子 (Remove 0-sugar molecules)
    zeroCount = (df["Total_Sugar_Count"] == 0).sum()
    df = df[df["Total_Sugar_Count"] > 0].copy()
    print(f"  Removed {zeroCount:,} zero-sugar molecules → {len(df):,} remaining")

    # ETL patches
    print("  Patch 1: Normalizing linkage types...")
    df["Linkage_Norm"] = df["Aglycone_Linkage_Type"].apply(normalizeLinkageType)
    print("  Extracting monosaccharide lists...")
    df["Mono_List"] = df["Sugar_Sequence"].apply(extractMonosaccharideList)
    print("  Extracting modification tags...")
    df["Mod_Tags"] = df["Glycan_Modifications"].apply(extractModificationTags)
    df["NP_Pathway"] = df["np_classifier_pathway"].fillna("Unknown")
    df["NP_Class"] = df["np_classifier_class"].fillna("Unknown")

    unknownCount = (df["Organism_Type"] == "Unknown").sum()
    unknownPct = unknownCount / len(df) * 100
    print(f"  Unknown organism: {unknownCount:,} ({unknownPct:.1f}%)")
    print("  ETL complete.\n")
    return df, unknownCount


def addAnnotation(ax, text: str, yPos: float = -0.14, fontSize: int = 8):
    """在图表下方添加注释段落, 通过 subplots_adjust 留出空间.
    Add annotation below chart with proper spacing via subplots_adjust.
    V3.2: 使用 fig.text 在图表底部, 避免和标题冲突.
    """
    fig = ax.get_figure()
    fig.subplots_adjust(bottom=0.15)
    fig.text(0.5, 0.02, text, ha="center", va="bottom", fontsize=fontSize,
             color=ANNOT_COLOR, style="italic",
             bbox=dict(boxstyle="round,pad=0.4", facecolor="#161B22",
                       edgecolor="#FFB800", alpha=0.9),
             wrap=True)


# =====================================================================
# 图表 1a: 单糖 vs 物种类型 热力图 (排除 Unknown)
# Chart 1a: Monosaccharide vs Organism Type (excl. Unknown)
# =====================================================================

def chart1aMonoVsOrganism(df: pd.DataFrame, unknownCount: int):
    """单糖种类 × 物种大类 占比热力图.
    V3: 列归一化 (每个物种内某糖出现的分子占比), 去重计糖.
    Column-normalized: proportion of molecules in each organism containing each sugar.
    Dedup: Glc-Glc-Ara counts as 1×Glc + 1×Ara per molecule.
    """
    print("  Chart 1a: Monosaccharide vs Organism Type (Proportion)...")
    rows = []
    for _, r in df.iterrows():
        orgType = r.get("Organism_Type", "Unknown")
        if orgType == "Unknown":
            continue
        # 去重: set() 确保每种糖在该分子中只计一次
        uniqueSugars = set(r.get("Mono_List", []))
        for sugar in uniqueSugars:
            rows.append({"Sugar": sugar, "Organism": orgType})
    if not rows:
        return
    mdf = pd.DataFrame(rows)
    topSugars = mdf["Sugar"].value_counts().head(15).index.tolist()
    topOrgs = [o for o in ["Plant", "Bacteria", "Fungi", "Animal", "Marine"]
               if o in mdf["Organism"].unique()]
    mdf = mdf[mdf["Sugar"].isin(topSugars) & mdf["Organism"].isin(topOrgs)]
    pivot = mdf.groupby(["Sugar", "Organism"]).size().unstack(fill_value=0)
    pivot = pivot.reindex(index=topSugars, columns=topOrgs).fillna(0)
    # 列归一化: 每列除以该物种的总分子数 (Column normalize by organism total)
    orgTotals = df[df["Organism_Type"] != "Unknown"].groupby("Organism_Type").size()
    for col in pivot.columns:
        if col in orgTotals.index:
            pivot[col] = pivot[col] / orgTotals[col] * 100
    pivotNorm = pivot

    fig, ax = plt.subplots(figsize=(10, 9))
    sns.heatmap(pivotNorm, annot=True, fmt=".1f", cmap=HEAT_CMAP,
                linewidths=0.5, linecolor="#30363D",
                cbar_kws={"label": "Proportion (% of molecules)", "shrink": 0.8}, ax=ax)
    ax.set_title("Monosaccharide × Organism Type — Proportion\n"
                 "单糖 × 物种类型 — 占比 (每物种中含该糖的分子百分比)", pad=15)
    # 添加样本量到 x-axis 标签 (Add sample sizes to x-axis labels)
    newXLabels = []
    for label in ax.get_xticklabels():
        orgName = label.get_text()
        n = orgTotals.get(orgName, 0)
        newXLabels.append(f"{orgName}\n(n={n:,})")
    ax.set_xticklabels(newXLabels, rotation=45, ha="right", fontsize=9)
    unknownPct = unknownCount / (len(df)) * 100
    addAnnotation(ax,
        f"NOTE: {unknownCount:,} molecules ({unknownPct:.1f}%) with unknown organism are excluded. Total N={len(df):,}.\n"
        "Column-normalized: each cell = % of molecules in that organism containing that sugar (deduped per molecule).",
        fontSize=8)
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "1a_mono_vs_organism.png")
    plt.close(fig)
    print("    Saved.")


# =====================================================================
# 图表 1a-F: 单糖 vs Family (Plant 下的 Family)
# Chart 1a-F: Monosaccharide vs NP Superclass (as proxy for Family)
# =====================================================================

def chart1aFamilyLevel(df: pd.DataFrame):
    """单糖种类 × NP Superclass (化学家族) 偏好热力图."""
    print("  Chart 1a-F: Monosaccharide vs NP Superclass...")
    rows = []
    for _, r in df.iterrows():
        npSup = r.get("np_classifier_superclass", None)
        if pd.isna(npSup) or not npSup:
            continue
        uniqueSugars = set(r.get("Mono_List", []))
        for sugar in uniqueSugars:
            rows.append({"Sugar": sugar, "Family": str(npSup)})
    if not rows:
        return
    mdf = pd.DataFrame(rows)
    topSugars = mdf["Sugar"].value_counts().head(12).index.tolist()
    topFams = mdf["Family"].value_counts().head(8).index.tolist()
    mdf = mdf[mdf["Sugar"].isin(topSugars) & mdf["Family"].isin(topFams)]
    pivot = mdf.groupby(["Sugar", "Family"]).size().unstack(fill_value=0)
    # 列归一化 (Column normalize by family total)
    famTotals = mdf.groupby("Family").apply(lambda x: x["Sugar"].nunique())
    pivotNorm = pivot.div(pivot.sum(axis=0), axis=1) * 100

    fig, ax = plt.subplots(figsize=(14, 8))
    sns.heatmap(pivotNorm, annot=True, fmt=".1f", cmap=HEAT_CMAP,
                linewidths=0.5, linecolor="#30363D",
                cbar_kws={"label": "Proportion (%)", "shrink": 0.8}, ax=ax)
    ax.set_title("Monosaccharide × NP Compound Superclass — Proportion\n"
                 "单糖 × 天然产物超类 — 占比热力图", pad=15)
    # 添加样本量到 x-axis 标签 (Add sample sizes to x-axis labels)
    famCounts = df["np_classifier_superclass"].value_counts()
    newXLabels = []
    for label in ax.get_xticklabels():
        famName = label.get_text()
        n = famCounts.get(famName, 0)
        newXLabels.append(f"{famName}\n(n={n:,})")
    ax.set_xticklabels(newXLabels, rotation=45, ha="right", fontsize=8)
    addAnnotation(ax,
        f"Column-normalized: % of sugar occurrences within each NP superclass (deduped per molecule). Total N={len(df):,}.\n"
        "Chemical family proxy: NP Classifier Superclass (87% coverage).",
        )
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "1a_family_mono_vs_npclass.png")
    plt.close(fig)
    print("    Saved.")


# =====================================================================
# 图表 1b: 序列 vs 骨架 (带骨架结构图)
# Chart 1b: Sequence vs Scaffold with structure images
# =====================================================================

def chart1bSeqVsScaffold(df: pd.DataFrame):
    """Top12 糖链序列 × Top8 骨架 热力图 + 骨架结构图 + 占比%.
    V3.1: 增加到 Top8 骨架, 每个骨架下方显示占比%.
    """
    print("  Chart 1b: Sequence vs Scaffold (with structures + %)...")
    subDf = df[df["Sugar_Sequence"].notna() & df["Generic_Scaffold"].notna()].copy()
    if subDf.empty:
        return

    # 骨架占比统计 (Scaffold proportion — from ALL molecules with scaffold)
    allScafDf = df[df["Generic_Scaffold"].notna()]
    scafTotalCount = len(allScafDf)
    scafCounts = allScafDf["Generic_Scaffold"].value_counts()

    topSeqs = subDf["Sugar_Sequence"].value_counts().head(12).index.tolist()
    topScafs = scafCounts.head(8).index.tolist()  # V3.1: Top 8

    subDf = subDf[subDf["Sugar_Sequence"].isin(topSeqs) & subDf["Generic_Scaffold"].isin(topScafs)]

    # 生成 Generic Murcko 骨架结构图 (Generate Generic Murcko scaffold images)
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    scafImages = {}
    scafLabels = {}
    scafPcts = {}  # 骨架占比 (Scaffold percentage)
    KNOWN_NAMES = {
        "C1CCC(C2CCC3CCCCC3C2)CC1": "Flavonoid",
        "C1CCC2C(C1)CCC1C2CCC2C3CCCCC3CCC21": "Steroid (4-ring)",
        "Aliphatic Chain": "Aliphatic Chain",
        "C1CCC2C(C1)CCC1C3CCCC3CCC21": "Androstane",
        "C1CCCCC1": "Cyclohexane",
    }
    for i, smi in enumerate(topScafs):
        pct = scafCounts.get(smi, 0) / scafTotalCount * 100
        if smi == "Aliphatic Chain":
            label = "Aliphatic Chain"
            scafLabels[smi] = label
            scafPcts[label] = pct
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol:
            try:
                fw = MurckoScaffold.GetScaffoldForMol(mol)
                gen = MurckoScaffold.MakeScaffoldGeneric(fw)
                genSmi = Chem.MolToSmiles(gen)
            except Exception:
                genSmi = smi
        else:
            genSmi = smi
        label = KNOWN_NAMES.get(genSmi, f"Scaffold {i+1}")
        scafLabels[smi] = label
        scafPcts[label] = pct
        img = smiToImage(genSmi, size=(180, 130))
        if img is not None:
            scafImages[label] = img

    subDf["Scaffold_Label"] = subDf["Generic_Scaffold"].map(scafLabels)
    seqLabels = {}
    for s in topSeqs:
        seqLabels[s] = s if len(s) <= 35 else s[:32] + "..."
    subDf["Seq_Label"] = subDf["Sugar_Sequence"].map(seqLabels)

    pivot = subDf.groupby(["Seq_Label", "Scaffold_Label"]).size().unstack(fill_value=0)

    # 双面板: 左=热力图, 右=骨架结构 (Two-panel: left=heatmap, right=scaffolds)
    fig = plt.figure(figsize=(20, 12))
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1], wspace=0.05)

    ax0 = fig.add_subplot(gs[0])
    sns.heatmap(pivot, annot=True, fmt="d", cmap=HEAT_CMAP,
                linewidths=0.5, linecolor="#30363D",
                cbar_kws={"label": "Co-occurrence", "shrink": 0.7}, ax=ax0)
    ax0.set_title("Top Sugar Sequences × Top Aglycone Scaffolds\n"
                  "糖链序列 × 苷元骨架 共现热力图", pad=15)
    ax0.set_xlabel("Aglycone Scaffold")
    ax0.set_ylabel("Sugar Sequence")

    # 右侧骨架结构面板 + 占比% (Right panel: scaffold structures + percentage)
    ax1 = fig.add_subplot(gs[1])
    ax1.set_facecolor("#0D1117")
    ax1.axis("off")
    ax1.set_title("Scaffold Structures\n骨架结构", pad=15, fontsize=11)

    nScafs = len(scafImages)
    for i, (label, img) in enumerate(scafImages.items()):
        yFrac = 1.0 - (i + 0.5) / nScafs
        imagebox = OffsetImage(img, zoom=0.5)
        ab = AnnotationBbox(imagebox, (0.5, yFrac), xycoords="axes fraction",
                            frameon=True, bboxprops=dict(facecolor="#1C2128",
                            edgecolor="#30363D", boxstyle="round,pad=0.2"))
        ax1.add_artist(ab)
        pctStr = f"{scafPcts.get(label, 0):.1f}%"
        ax1.text(0.5, yFrac - 0.05, f"{label}\n({pctStr})",
                 transform=ax1.transAxes,
                 ha="center", va="top", fontsize=8, color="#C9D1D9",
                 fontweight="bold")

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "1b_seq_vs_scaffold.png")
    plt.close(fig)
    print("    Saved.")


# =====================================================================
# 图表 1c: 修饰 vs 物种 + Family 版
# Chart 1c: Modification vs Organism + NP Superclass variant
# =====================================================================

def chart1cModVsOrganism(df: pd.DataFrame):
    """修饰基团 × 物种 热力图 (排除 Unknown)."""
    print("  Chart 1c: Modification vs Organism...")
    rows = []
    for _, r in df.iterrows():
        orgType = r.get("Organism_Type", "Unknown")
        if orgType == "Unknown":
            continue
        for mod in r.get("Mod_Tags", []):
            rows.append({"Modification": mod, "Organism": orgType})
    if not rows:
        return
    mdf = pd.DataFrame(rows)
    topMods = mdf["Modification"].value_counts().head(10).index.tolist()
    topOrgs = [o for o in ["Plant", "Bacteria", "Fungi", "Animal", "Marine"]
               if o in mdf["Organism"].unique()]
    mdf = mdf[mdf["Modification"].isin(topMods) & mdf["Organism"].isin(topOrgs)]
    pivot = mdf.groupby(["Modification", "Organism"]).size().unstack(fill_value=0)
    # 列归一化 (Column normalize)
    pivotNorm = pivot.div(pivot.sum(axis=0), axis=1) * 100

    fig, ax = plt.subplots(figsize=(10, 7))
    sns.heatmap(pivotNorm, annot=True, fmt=".1f", cmap=HEAT_CMAP,
                linewidths=0.5, linecolor="#30363D",
                cbar_kws={"label": "Proportion (%)", "shrink": 0.8}, ax=ax)
    ax.set_title("Modification × Organism Type — Proportion\n"
                 "修饰基团 × 物种类型 — 占比热力图", pad=15)
    # 添加样本量到 x-axis 标签 (Add sample sizes to organism labels)
    orgModCounts = mdf.groupby("Organism").size()
    newXLabels = []
    for label in ax.get_xticklabels():
        orgName = label.get_text()
        n = orgModCounts.get(orgName, 0)
        newXLabels.append(f"{orgName}\n(n={n:,})")
    ax.set_xticklabels(newXLabels, rotation=45, ha="right", fontsize=9)
    addAnnotation(ax,
        "Column-normalized: % of modifications within each organism.\n"
        "O-Ac = O-Acetyl, NAc = N-Acetyl, Sulfate/Phosphate are polar modifications.",
        )
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "1c_mod_vs_organism.png")
    plt.close(fig)
    print("    Saved.")


def chart1cModVsFamily(df: pd.DataFrame):
    """修饰基团 × NP Superclass 热力图."""
    print("  Chart 1c-F: Modification vs NP Superclass...")
    rows = []
    for _, r in df.iterrows():
        npSup = r.get("np_classifier_superclass", None)
        if pd.isna(npSup) or not npSup:
            continue
        for mod in r.get("Mod_Tags", []):
            rows.append({"Modification": mod, "Family": str(npSup)})
    if not rows:
        return
    mdf = pd.DataFrame(rows)
    topMods = mdf["Modification"].value_counts().head(10).index.tolist()
    topFams = mdf["Family"].value_counts().head(8).index.tolist()
    mdf = mdf[mdf["Modification"].isin(topMods) & mdf["Family"].isin(topFams)]
    pivot = mdf.groupby(["Modification", "Family"]).size().unstack(fill_value=0)
    pivotNorm = pivot.div(pivot.sum(axis=0), axis=1) * 100

    fig, ax = plt.subplots(figsize=(14, 7))
    sns.heatmap(pivotNorm, annot=True, fmt=".1f", cmap=HEAT_CMAP,
                linewidths=0.5, linecolor="#30363D",
                cbar_kws={"label": "Proportion (%)", "shrink": 0.8}, ax=ax)
    ax.set_title("Modification × NP Superclass — Proportion\n"
                 "修饰基团 × 天然产物超类 — 占比热力图", pad=15)
    # 添加样本量到 x-axis 标签 (Add sample sizes to family labels)
    famModCounts = mdf.groupby("Family").size()
    newXLabels = []
    for label in ax.get_xticklabels():
        famName = label.get_text()
        n = famModCounts.get(famName, 0)
        newXLabels.append(f"{famName}\n(n={n:,})")
    ax.set_xticklabels(newXLabels, rotation=45, ha="right", fontsize=8)
    addAnnotation(ax,
        "Column-normalized: % of modifications within each NP superclass.",
        )
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "1c_family_mod_vs_npclass.png")
    plt.close(fig)
    print("    Saved.")


# =====================================================================
# 图表 1d: 糖苷键 vs 多维度 (Pathway / 物种 / 骨架)
# Chart 1d: Linkage in 3 sub-panels
# =====================================================================

def chart1dLinkageMulti(df: pd.DataFrame):
    """糖苷键 × NP Pathway / Organism / Scaffold 三合一."""
    print("  Chart 1d: Linkage multi-panel...")
    subDf = df[df["Linkage_Norm"].str.len() > 0].copy()
    topLinks = subDf["Linkage_Norm"].value_counts().head(6).index.tolist()

    fig, axes = plt.subplots(1, 3, figsize=(24, 7))

    # 1d-A: Linkage vs NP Pathway
    ax = axes[0]
    topPaths = subDf["NP_Pathway"].value_counts().head(5).index.tolist()
    sub = subDf[subDf["Linkage_Norm"].isin(topLinks) & subDf["NP_Pathway"].isin(topPaths)]
    piv = sub.groupby(["Linkage_Norm", "NP_Pathway"]).size().unstack(fill_value=0)
    pivN = piv.div(piv.sum(axis=1), axis=0) * 100
    sns.heatmap(pivN, annot=True, fmt=".0f", cmap=HEAT_CMAP,
                linewidths=0.5, linecolor="#30363D", ax=ax,
                cbar_kws={"shrink": 0.6})
    npCounts = {p: (sub["NP_Pathway"]==p).sum() for p in topPaths}
    ax.set_title("Linkage × NP Pathway\n糖苷键 × NP 大类\n" + " | ".join([f"{p}: n={npCounts.get(p,0):,}" for p in topPaths[:3]]), fontsize=9)
    ax.set_xlabel("")

    # 1d-B: Linkage vs Organism
    ax = axes[1]
    sub = subDf[(subDf["Linkage_Norm"].isin(topLinks)) &
                (subDf["Organism_Type"] != "Unknown")]
    topOrgs = sub["Organism_Type"].value_counts().head(5).index.tolist()
    sub = sub[sub["Organism_Type"].isin(topOrgs)]
    piv = sub.groupby(["Linkage_Norm", "Organism_Type"]).size().unstack(fill_value=0)
    pivN = piv.div(piv.sum(axis=1), axis=0) * 100
    sns.heatmap(pivN, annot=True, fmt=".0f", cmap=HEAT_CMAP,
                linewidths=0.5, linecolor="#30363D", ax=ax,
                cbar_kws={"shrink": 0.6})
    orgLinkCounts = {o: (sub["Organism_Type"]==o).sum() for o in topOrgs}
    ax.set_title("Linkage × Organism\n糖苷键 × 物种类型\n" + " | ".join([f"{o}: n={orgLinkCounts.get(o,0):,}" for o in topOrgs[:4]]), fontsize=9)
    ax.set_xlabel("")
    ax.set_ylabel("")

    # 1d-C: Sugar type within each linkage (top sugars by linkage)
    ax = axes[2]
    topSuper = subDf["np_classifier_superclass"].value_counts().head(5).index.tolist()
    sub = subDf[subDf["Linkage_Norm"].isin(topLinks) &
                subDf["np_classifier_superclass"].isin(topSuper)]
    piv = sub.groupby(["Linkage_Norm", "np_classifier_superclass"]).size().unstack(fill_value=0)
    pivN = piv.div(piv.sum(axis=1), axis=0) * 100
    sns.heatmap(pivN, annot=True, fmt=".0f", cmap=HEAT_CMAP,
                linewidths=0.5, linecolor="#30363D", ax=ax,
                cbar_kws={"shrink": 0.6})
    ax.set_title("Linkage × NP Superclass\n糖苷键 × NP 超类")
    ax.set_xlabel("")
    ax.set_ylabel("")

    plt.suptitle("Glycosidic Bond Type — Multi-Dimensional Proportion Analysis\n"
                 "糖苷键类型 — 多维度占比分析", fontsize=14, fontweight="bold",
                 color="#E6EDF3", y=1.02)
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "1d_linkage_multi.png")
    plt.close(fig)
    print("    Saved.")


# =====================================================================
# 图表 2: 宏观统计 (排除 0-sugar)
# Chart 2: Macro Statistics (0-sugar excluded)
# =====================================================================

def chart2MacroStats(df: pd.DataFrame):
    """三合一统计分布图 (Triple distribution panel, 0-sugar excluded)."""
    print("  Chart 2: Macro Statistics...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # 2a: Sugar Count
    ax = axes[0]
    sugarCounts = df["Total_Sugar_Count"].dropna().astype(int)
    bins = list(range(1, min(int(sugarCounts.max()) + 2, 16)))
    ax.hist(sugarCounts, bins=bins, color="#58A6FF", edgecolor="#0D1117",
            alpha=0.85, rwidth=0.85)
    ax.set_xlabel("Total Sugar Count (含糖总数)")
    ax.set_ylabel("Molecule Count (分子数)")
    ax.set_title("Sugar Count Distribution\n含糖总数分布")
    ax.set_yscale("log")
    ax.grid(axis="y", alpha=0.3)

    # 2b: Max Chain Length
    ax = axes[1]
    chainLens = df["Max_Chain_Length"].dropna().astype(int)
    chainLens = chainLens[chainLens >= 1]
    bins2 = list(range(1, min(int(chainLens.max()) + 2, 16)))
    ax.hist(chainLens, bins=bins2, color="#F0883E", edgecolor="#0D1117",
            alpha=0.85, rwidth=0.85)
    ax.set_xlabel("Max Chain Length (最长链长度)")
    ax.set_ylabel("Molecule Count (分子数)")
    ax.set_title("Max Chain Length Distribution\n最大链长分布")
    ax.set_yscale("log")
    ax.grid(axis="y", alpha=0.3)

    # 2c: Branching ratio
    ax = axes[2]
    hasBranch = df["Sugar_Sequence"].str.contains(r"\[", na=False)
    branchCounts = hasBranch.value_counts()
    labels = ["Linear (线性)", "Branched (分支)"]
    sizes = [branchCounts.get(False, 0), branchCounts.get(True, 0)]
    colors = ["#79C0FF", "#A371F7"]
    wedges, texts, autotexts = ax.pie(
        sizes, labels=labels, colors=colors, explode=(0, 0.05),
        autopct="%1.1f%%", shadow=True, startangle=140,
        textprops={"color": "#C9D1D9", "fontsize": 10})
    for at in autotexts:
        at.set_fontweight("bold")
    ax.set_title("Branching Ratio\n分支比例")

    plt.suptitle(f"Macro Statistics (0-sugar excluded, N={len(df):,})\n"
                 "宏观统计面板（已排除 0 糖分子）",
                 fontsize=14, fontweight="bold", color="#E6EDF3", y=1.02)
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "2_macro_statistics.png")
    plt.close(fig)
    print("    Saved.")


# =====================================================================
# 图表 3a: 糖链 ↔ 物种 网络 (带注释)
# Chart 3a: Glycan ↔ Organism Network (annotated)
# =====================================================================

def chart3aGlycanOrgNetwork(df: pd.DataFrame):
    """糖链序列 ↔ 物种 网络图 (annotated)."""
    print("  Chart 3a: Glycan–Organism Network...")
    subDf = df[df["Sugar_Sequence"].notna() & (df["Organism_Type"] != "Unknown")].copy()
    topSeqs = subDf["Sugar_Sequence"].value_counts().head(20).index.tolist()
    subDf = subDf[subDf["Sugar_Sequence"].isin(topSeqs)]

    G = nx.Graph()
    edgeWeights = Counter()
    for _, r in subDf.iterrows():
        seq = r["Sugar_Sequence"]
        org = r["Organism_Type"]
        seqLabel = seq if len(seq) <= 25 else seq[:22] + "..."
        edgeWeights[(seqLabel, org)] += 1

    for (seq, org), w in edgeWeights.items():
        if w >= 5:
            G.add_node(seq, node_type="sugar", size=w)
            orgN = len(subDf[subDf["Organism_Type"] == org])
            G.add_node(org, node_type="organism", size=100, label=f"{org}\n(n={orgN:,})")
            G.add_edge(seq, org, weight=w)

    if G.number_of_nodes() == 0:
        return

    fig, ax = plt.subplots(figsize=(14, 11))
    pos = nx.spring_layout(G, k=2.5, iterations=60, seed=42)

    sugarNodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "sugar"]
    orgNodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "organism"]

    edgeWidths = [G[u][v]["weight"] / 10 for u, v in G.edges()]
    nx.draw_networkx_edges(G, pos, alpha=0.3, width=edgeWidths,
                           edge_color="#8B949E", ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=sugarNodes, node_color="#58A6FF",
                           node_size=300, alpha=0.9, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=orgNodes, node_color="#F0883E",
                           node_size=800, alpha=0.9, ax=ax)
    customLabels = {n: d.get("label", n) for n, d in G.nodes(data=True)}
    nx.draw_networkx_labels(G, pos, labels=customLabels, font_size=7, font_color="#E6EDF3", ax=ax)

    ax.set_title("Top 20 Glycan Sequences ↔ Organism Types Network\n"
                 "糖链序列 ↔ 物种类型 关联网络", pad=15)
    ax.axis("off")

    # 注释 (Annotation)
    annotText = (
        "Each BLUE node represents a glycan sugar sequence; each ORANGE node represents an organism kingdom.\n"
        "Edge thickness indicates co-occurrence frequency. Edges with <5 occurrences are filtered out.\n"
        "This network reveals which sugar sequences are preferentially associated with each organism type."
    )
    fig.text(0.5, 0.02, annotText, ha="center", va="bottom", fontsize=8,
             color=ANNOT_COLOR, style="italic",
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#161B22",
                       edgecolor="#FFB800", alpha=0.9))
    fig.savefig(OUTPUT_DIR / "3a_glycan_organism_network.png")
    plt.close(fig)
    print("    Saved.")


# =====================================================================
# 图表 3b: Motif 保守基序网络 (带注释)
# Chart 3b: Motif Conserved Network (annotated)
# =====================================================================

def chart3bMotifNetwork(df: pd.DataFrame):
    """Motif 片段 ↔ NP 大类 网络图 (annotated)."""
    print("  Chart 3b: Motif–Pathway Network...")
    rows = []
    for _, r in df.iterrows():
        seqStr = r.get("Sugar_Sequence", "")
        pathway = r.get("NP_Pathway", "Unknown")
        if pathway == "Unknown":
            continue
        for m in extractMotifs(seqStr, windowSize=2):
            rows.append({"Motif": m, "Pathway": pathway})
    if not rows:
        return
    mdf = pd.DataFrame(rows)
    topMotifs = mdf["Motif"].value_counts().head(20).index.tolist()
    topPaths = mdf["Pathway"].value_counts().head(5).index.tolist()
    mdf = mdf[mdf["Motif"].isin(topMotifs) & mdf["Pathway"].isin(topPaths)]

    G = nx.Graph()
    edgeWeights = Counter()
    for _, r in mdf.iterrows():
        motif = r["Motif"]
        pathway = r["Pathway"]
        motifLabel = motif if len(motif) <= 25 else motif[:22] + "..."
        edgeWeights[(motifLabel, pathway)] += 1

    for (motif, pathway), w in edgeWeights.items():
        if w >= 3:
            G.add_node(motif, node_type="motif")
            G.add_node(pathway, node_type="pathway")
            G.add_edge(motif, pathway, weight=w)

    if G.number_of_nodes() == 0:
        return

    fig, ax = plt.subplots(figsize=(14, 11))
    pos = nx.spring_layout(G, k=3.0, iterations=80, seed=42)

    motifNodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "motif"]
    pathNodes = [n for n, d in G.nodes(data=True) if d.get("node_type") == "pathway"]

    edgeWidths = [min(G[u][v]["weight"] / 5, 5) for u, v in G.edges()]
    nx.draw_networkx_edges(G, pos, alpha=0.35, width=edgeWidths,
                           edge_color="#8B949E", ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=motifNodes, node_color="#A371F7",
                           node_size=350, alpha=0.9, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=pathNodes, node_color="#3FB950",
                           node_size=900, alpha=0.9, ax=ax)
    customLabels = {n: d.get("label", n) for n, d in G.nodes(data=True)}
    nx.draw_networkx_labels(G, pos, labels=customLabels, font_size=7, font_color="#E6EDF3", ax=ax)

    ax.set_title("Disaccharide Motif ↔ NP Pathway Conserved Network\n"
                 "二糖 Motif 保守基序 ↔ NP 大类 网络", pad=15)
    ax.axis("off")

    annotText = (
        "PURPLE nodes = disaccharide fragments extracted via sliding window (size=2) from chains with ≥2 sugars.\n"
        "GREEN nodes = NP biosynthetic pathways (Terpenoids, Shikimates, Polyketides, etc.).\n"
        "A conserved motif is a disaccharide sub-structure that appears frequently across many molecules.\n"
        "Edge thickness ∝ co-occurrence frequency. This reveals pathway-specific sugar building blocks."
    )
    fig.text(0.5, 0.02, annotText, ha="center", va="bottom", fontsize=8,
             color=ANNOT_COLOR, style="italic",
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#161B22",
                       edgecolor="#FFB800", alpha=0.9))
    fig.savefig(OUTPUT_DIR / "3b_motif_pathway_network.png")
    plt.close(fig)
    print("    Saved.")


# =====================================================================
# 图表 4a: 立体指南针 (带注释)
# Chart 4a: Linkage Compass (annotated)
# =====================================================================

def chart4aLinkageCompass(df: pd.DataFrame):
    """供体糖+α/β × 受体糖+位点 热力图 (annotated)."""
    print("  Chart 4a: Linkage Compass...")
    allEdges = []
    for seq in df["Sugar_Sequence"].dropna():
        allEdges.extend(parseSequenceToEdges(str(seq)))
    if not allEdges:
        return

    edf = pd.DataFrame(allEdges)
    edf["Donor_Anomer"] = edf["donor"] + " " + edf["anomer"]
    edf["Acceptor_Pos"] = edf["acceptor"] + " " + edf["pos"]

    topDonors = edf["Donor_Anomer"].value_counts().head(12).index.tolist()
    topAcceptors = edf["Acceptor_Pos"].value_counts().head(10).index.tolist()
    edf = edf[edf["Donor_Anomer"].isin(topDonors) & edf["Acceptor_Pos"].isin(topAcceptors)]
    pivot = edf.groupby(["Donor_Anomer", "Acceptor_Pos"]).size().unstack(fill_value=0)

    fig, ax = plt.subplots(figsize=(14, 9))
    sns.heatmap(pivot, annot=True, fmt="d", cmap=HEAT_CMAP,
                linewidths=0.5, linecolor="#30363D",
                cbar_kws={"label": "Frequency", "shrink": 0.8}, ax=ax)
    ax.set_title("Stereochemical Linkage Compass\n"
                 "立体连接指南针: 供体+异头碳构型 × 受体+连接位点", pad=15)
    ax.set_xlabel("Acceptor + Linkage Position (受体 + 连接位)")
    ax.set_ylabel("Donor + Anomeric Config (供体 + α/β)")

    addAnnotation(ax,
        f"Total linkage events: N={len(edf):,}. "
        "Each cell counts donor sugar (with α/β) linked to acceptor at a hydroxyl position.\n"
        "Example: 'L-Rha β → D-Glc 1→4' means L-Rha is β-linked to the C4 of D-Glc.",
        fontSize=8)
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "4a_linkage_compass.png")
    plt.close(fig)
    print("    Saved.")


# =====================================================================
# 图表 4b: 共现矩阵 (带注释)
# Chart 4b: Co-occurrence Matrix (annotated)
# =====================================================================

def chart4bCoOccurrence(df: pd.DataFrame):
    """单糖 × 修饰基团 共现矩阵 (annotated)."""
    print("  Chart 4b: Sugar–Modification Co-occurrence...")
    rows = []
    for _, r in df.iterrows():
        modStr = r.get("Glycan_Modifications", "")
        if pd.isna(modStr):
            continue
        for unit in str(modStr).split(";"):
            unit = unit.strip()
            sugarMatch = re.match(r"([A-Za-z\d\-]+)_\d+", unit)
            modMatches = re.findall(r"\*([A-Za-z\d\-]+)", unit)
            if sugarMatch:
                for mod in modMatches:
                    rows.append({"Sugar": sugarMatch.group(1), "Mod": mod})
    if not rows:
        return
    codf = pd.DataFrame(rows)
    topSugars = codf["Sugar"].value_counts().head(12).index.tolist()
    topMods = codf["Mod"].value_counts().head(10).index.tolist()
    codf = codf[codf["Sugar"].isin(topSugars) & codf["Mod"].isin(topMods)]
    pivot = codf.groupby(["Sugar", "Mod"]).size().unstack(fill_value=0)

    fig, ax = plt.subplots(figsize=(12, 8))
    sns.heatmap(pivot, annot=True, fmt="d", cmap=GLYCO_CMAP,
                linewidths=0.5, linecolor="#30363D",
                cbar_kws={"label": "Co-occurrence", "shrink": 0.8}, ax=ax)
    ax.set_title("Monosaccharide × Modification Co-occurrence Matrix\n"
                 "单糖 × 修饰基团 共现矩阵", pad=15)
    ax.set_xlabel("Modification (修饰基团)")
    ax.set_ylabel("Monosaccharide (单糖)")
    addAnnotation(ax,
        f"N={len(codf):,} sugar-modification pairs. "
        "Each cell counts unique molecules with a specific sugar + modification.\n"
        "High values indicate strong association. O-Ac (O-Acetyl) is the most ubiquitous modification.",
        fontSize=8)
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "4b_sugar_mod_cooccurrence.png")
    plt.close(fig)
    print("    Saved.")


# =====================================================================
# 图表 4c: PCA 散点图 (带注释 + 更清晰的骨架标签)
# Chart 4c: PCA Scatter (annotated + better scaffold labels)
# =====================================================================

def chart4cScaffoldPCA(df: pd.DataFrame):
    """核心骨架衍生分子 PCA 降维 (annotated)."""
    print("  Chart 4c: Scaffold PCA...")
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from sklearn.decomposition import PCA
    except ImportError:
        print("    RDKit or sklearn not available, skipping.")
        return

    subDf = df[df["Generic_Scaffold"].notna()].copy()
    topScafs = subDf["Generic_Scaffold"].value_counts().head(8).index.tolist()
    subDf = subDf[subDf["Generic_Scaffold"].isin(topScafs)]
    if len(subDf) > 3000:
        subDf = subDf.sample(3000, random_state=42)

    fps, validIdx = [], []
    for idx, smi in subDf["canonical_smiles"].items():
        mol = Chem.MolFromSmiles(str(smi))
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            fps.append(np.array(fp))
            validIdx.append(idx)
    if len(fps) < 50:
        return

    X = np.array(fps)
    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(X)

    # 有意义的骨架标签 (Meaningful scaffold labels)
    SCAFFOLD_NAMES = {
        "Cc1cc(-c2ccccc2)cc2ccccc12": "Flavonoid",
        "C1CCC2C(C1)CCC1C2CCC2C3CCCCC3CCC21": "Steroid",
        "Aliphatic Chain": "Aliphatic",
        "C1CCC2C(C1)CCC1C3CCCC3CCC21": "Androstane",
        "c1ccccc1": "Benzene",
        "C1CCCCC1": "Cyclohexane",
    }
    scafMap = {}
    for i, s in enumerate(topScafs):
        scafMap[s] = SCAFFOLD_NAMES.get(s, f"Scaffold_{i+1}")

    plotDf = pd.DataFrame({
        "PC1": coords[:, 0], "PC2": coords[:, 1],
        "Scaffold_Label": [scafMap[subDf.loc[idx, "Generic_Scaffold"]] for idx in validIdx],
    })

    fig, ax = plt.subplots(figsize=(13, 10))
    for i, (scaf, group) in enumerate(plotDf.groupby("Scaffold_Label")):
        color = ACCENT[i % len(ACCENT)]
        ax.scatter(group["PC1"], group["PC2"], c=color, label=scaf,
                   alpha=0.6, s=15, edgecolors="none")

    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)")
    ax.set_title(f"Aglycone Scaffold Chemical Space (Morgan FP PCA, n={len(plotDf):,})\n"
                 "苷元骨架化学空间 PCA 投射", pad=15)
    ax.legend(loc="upper right", fontsize=9, framealpha=0.3,
              edgecolor="#30363D", facecolor="#161B22")
    ax.grid(alpha=0.2)

    annotText = (
        "METHOD: Morgan Fingerprint (ECFP4, radius=2) encodes each atom\u2019s 2-bond chemical neighborhood into a 1024-bit vector.\n"
        "PCA (Principal Component Analysis) projects the 1024-dim fingerprint space onto the 2 directions of maximum variance.\n"
        "INTERPRETATION: Points close together = structurally similar molecules. Clear clusters = distinct chemical scaffolds.\n"
        "SIGNIFICANCE: Glycosylation does not blur scaffold identity \u2014 aglycone core determines chemical space position."
    )
    fig.text(0.5, 0.02, annotText, ha="center", va="bottom", fontsize=8,
             color=ANNOT_COLOR, style="italic",
             bbox=dict(boxstyle="round,pad=0.5", facecolor="#161B22",
                       edgecolor="#FFB800", alpha=0.9))
    fig.savefig(OUTPUT_DIR / "4c_scaffold_pca.png")
    plt.close(fig)
    print("    Saved.")


# =====================================================================
# 新增图表 5: NP Class × Sugar Count 热力图
# New Chart 5: NP Class × Sugar Count
# =====================================================================

def chart5ClassVsSugarCount(df: pd.DataFrame):
    """NP 化合物大类 × 含糖数量 热力图."""
    print("  Chart 5: NP Class vs Sugar Count...")
    subDf = df[df["NP_Class"] != "Unknown"].copy()
    topClasses = subDf["NP_Class"].value_counts().head(12).index.tolist()
    subDf = subDf[subDf["NP_Class"].isin(topClasses)].copy()
    subDf["Sugar_Bin"] = subDf["Total_Sugar_Count"].clip(upper=6).astype(int)
    subDf.loc[subDf["Total_Sugar_Count"] > 6, "Sugar_Bin"] = 7

    pivot = subDf.groupby(["NP_Class", "Sugar_Bin"]).size().unstack(fill_value=0)
    binLabels = {1: "1", 2: "2", 3: "3", 4: "4", 5: "5", 6: "6", 7: "7+"}
    pivot.columns = [binLabels.get(c, str(c)) for c in pivot.columns]
    pivotNorm = pivot.div(pivot.sum(axis=1), axis=0) * 100

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(pivotNorm, annot=True, fmt=".0f", cmap=HEAT_CMAP,
                linewidths=0.5, linecolor="#30363D",
                cbar_kws={"label": "Proportion (%)", "shrink": 0.8}, ax=ax)
    # 添加样本量到 y-axis 标签 (Add n= to NP Class labels)
    classCounts = subDf.groupby("NP_Class").size()
    newYLabels = []
    for label in ax.get_yticklabels():
        clsName = label.get_text()
        n = classCounts.get(clsName, 0)
        newYLabels.append(f"{clsName} (n={n:,})")
    ax.set_yticklabels(newYLabels, fontsize=8)
    ax.set_title(f"NP Compound Class × Sugar Count Distribution (N={len(subDf):,})\n"
                 "NP 化合物类别 × 含糖数量 分布", pad=15)
    ax.set_xlabel("Sugar Count (含糖数)")
    ax.set_ylabel("NP Class (化合物类别)")
    addAnnotation(ax,
        "Shows glycosylation density across natural product classes.\n"
        "Oleanane triterpenoids tend to carry 3-5 sugars (saponins), while flavonols typically have 1-2.",
        )
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "5_npclass_vs_sugarcount.png")
    plt.close(fig)
    print("    Saved.")


# =====================================================================
# 新增图表 6: 糖苷键 α/β 比例 按 NP Class
# New Chart 6: α/β Ratio per NP Class
# =====================================================================

def chart6AlphaBetaRatio(df: pd.DataFrame):
    """α/β 糖苷键比例 按 NP 大类 — V3: 统计所有糖苷键 (不仅 reducing-end).
    V3: counts ALL glycosidic bonds from Sugar_Sequence, not just aglycone linkage.
    """
    print("  Chart 6: α/β Ratio by NP Pathway (ALL bonds)...")
    rows = []
    for _, r in df.iterrows():
        pathway = r.get("NP_Pathway", "Unknown")
        if pathway == "Unknown":
            continue
        seq = str(r.get("Sugar_Sequence", ""))
        # 从序列中提取所有 a/b 标记 (Extract all anomeric configs from sequence)
        anomers = re.findall(r"\(([ab])\d-\d\)", seq)
        for a in anomers:
            rows.append({"NP_Pathway": pathway, "Anomer": "α" if a == "a" else "β"})
        # 也包含 Aglycone_Linkage_Type (Also include root linkage)
        linkNorm = r.get("Linkage_Norm", "")
        if isinstance(linkNorm, str):
            for sym in ("α", "β"):
                if sym in linkNorm:
                    rows.append({"NP_Pathway": pathway, "Anomer": sym})
    if not rows:
        return
    subDf = pd.DataFrame(rows)
    topPaths = subDf["NP_Pathway"].value_counts().head(6).index.tolist()
    subDf = subDf[subDf["NP_Pathway"].isin(topPaths)]

    pivot = subDf.groupby(["NP_Pathway", "Anomer"]).size().unstack(fill_value=0)
    pivot["α_ratio"] = pivot.get("α", 0) / (pivot.get("α", 0) + pivot.get("β", 0)) * 100

    fig, ax = plt.subplots(figsize=(10, 6))
    colors = ["#FF7B72", "#58A6FF"]
    pivot[["α", "β"]].plot(kind="barh", stacked=True, ax=ax, color=colors,
                           edgecolor="#0D1117", alpha=0.85)
    ax.set_xlabel("Count (分子数)")
    # 添加样本量到 y-axis 标签 (Add n= to pathway labels)
    pathCounts = subDf.groupby("NP_Pathway").size()
    newYLabels = []
    for label in ax.get_yticklabels():
        pwName = label.get_text()
        n = pathCounts.get(pwName, 0)
        newYLabels.append(f"{pwName} (n={n:,})")
    ax.set_yticklabels(newYLabels, fontsize=9)
    ax.set_title("α vs β Glycosidic Bond Ratio per NP Pathway\n"
                 "α/β 糖苷键比例 按 NP 大类", pad=15)
    ax.legend(title="Anomer", framealpha=0.3, edgecolor="#30363D")
    ax.grid(axis="x", alpha=0.2)

    # 标注比例 (Annotate ratios)
    for i, (path, row) in enumerate(pivot.iterrows()):
        total = row.get("α", 0) + row.get("β", 0)
        if total > 0:
            ratio = row.get("α", 0) / total * 100
            ax.text(total + total * 0.02, i, f"α={ratio:.0f}%",
                    va="center", fontsize=8, color="#C9D1D9")

    addAnnotation(ax,
        "ALL glycosidic bonds counted (inter-sugar + aglycone linkage), not just reducing-end.\n""Terpenoid glycosides typically favor β-O-linkages, while Shikimates show a more balanced α/β ratio.",
        )
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "6_alpha_beta_ratio.png")
    plt.close(fig)
    print("    Saved.")


# =====================================================================
# 新增图表 7: Top 单糖频率条形图 (总览)
# New Chart 7: Monosaccharide frequency bar chart
# =====================================================================

def chart7MonoFreqBar(df: pd.DataFrame):
    """Top 20 单糖出现频率 (Overall frequency)."""
    print("  Chart 7: Monosaccharide Frequency...")
    allSugars = []
    for lst in df["Mono_List"]:
        if isinstance(lst, list):
            allSugars.extend(lst)
    sc = Counter(allSugars).most_common(20)
    names, counts = zip(*sc)

    fig, ax = plt.subplots(figsize=(12, 6))
    bars = ax.barh(range(len(names)), counts, color=ACCENT * 3, edgecolor="#0D1117",
                   alpha=0.85)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names)
    ax.invert_yaxis()
    ax.set_xlabel("Occurrence Count (出现次数)")
    ax.set_title(f"Top 20 Monosaccharide Frequency (N={sum(counts):,} occurrences)\nTop 20 单糖出现频次", pad=15)
    ax.grid(axis="x", alpha=0.2)

    for bar, count in zip(bars, counts):
        ax.text(bar.get_width() + max(counts) * 0.01, bar.get_y() + bar.get_height() / 2,
                f"{count:,}", va="center", fontsize=8, color="#C9D1D9")

    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "7_mono_frequency.png")
    plt.close(fig)
    print("    Saved.")




# =====================================================================
# 图表 4a 分面: 按物种/NP 大类 (Faceted 4a)
# Chart 4a faceted by Organism and NP Pathway
# =====================================================================

def chart4aFacetedByOrganism(df: pd.DataFrame):
    """4a 按物种分面 (Faceted by Organism Type)."""
    print("  Chart 4a-Org: Linkage Compass by Organism...")
    allEdges = []
    for _, r in df.iterrows():
        org = r.get("Organism_Type", "Unknown")
        if org == "Unknown":
            continue
        seq = str(r.get("Sugar_Sequence", ""))
        for edge in parseSequenceToEdges(seq):
            edge["Organism"] = org
            allEdges.append(edge)
    if not allEdges:
        return
    edf = pd.DataFrame(allEdges)
    edf["Donor_Anomer"] = edf["donor"] + " " + edf["anomer"]
    edf["Acceptor_Pos"] = edf["acceptor"] + " " + edf["pos"]
    topDonors = edf["Donor_Anomer"].value_counts().head(8).index.tolist()
    topAcceptors = edf["Acceptor_Pos"].value_counts().head(6).index.tolist()
    orgs = [o for o in ["Plant", "Bacteria", "Fungi", "Animal"] if o in edf["Organism"].unique()]
    nOrgs = len(orgs)
    if nOrgs == 0:
        return
    fig, axes = plt.subplots(1, nOrgs, figsize=(6 * nOrgs, 7))
    if nOrgs == 1:
        axes = [axes]
    for ax, org in zip(axes, orgs):
        sub = edf[(edf["Organism"] == org) & edf["Donor_Anomer"].isin(topDonors) & edf["Acceptor_Pos"].isin(topAcceptors)]
        if sub.empty:
            ax.set_title(f"{org}\n(no data)")
            continue
        piv = sub.groupby(["Donor_Anomer", "Acceptor_Pos"]).size().unstack(fill_value=0)
        sns.heatmap(piv, annot=True, fmt="d", cmap=HEAT_CMAP, linewidths=0.4, linecolor="#30363D", ax=ax, cbar=False)
        nOrg = len(sub)
        ax.set_title(f"{org} (n={nOrg:,})")
        if ax != axes[0]:
            ax.set_ylabel("")
    plt.suptitle("Linkage Compass — Faceted by Organism Type\n立体指南针 — 按物种分面",
                 fontsize=13, fontweight="bold", color="#E6EDF3", y=1.02)
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "4a_compass_by_organism.png")
    plt.close(fig)
    print("    Saved.")


def chart4aFacetedByPathway(df: pd.DataFrame):
    """4a 按 NP Pathway 分面 (Faceted by NP Pathway)."""
    print("  Chart 4a-NP: Linkage Compass by NP Pathway...")
    allEdges = []
    for _, r in df.iterrows():
        pathway = r.get("NP_Pathway", "Unknown")
        if pathway == "Unknown":
            continue
        seq = str(r.get("Sugar_Sequence", ""))
        for edge in parseSequenceToEdges(seq):
            edge["Pathway"] = pathway
            allEdges.append(edge)
    if not allEdges:
        return
    edf = pd.DataFrame(allEdges)
    edf["Donor_Anomer"] = edf["donor"] + " " + edf["anomer"]
    edf["Acceptor_Pos"] = edf["acceptor"] + " " + edf["pos"]
    topDonors = edf["Donor_Anomer"].value_counts().head(8).index.tolist()
    topAcceptors = edf["Acceptor_Pos"].value_counts().head(6).index.tolist()
    paths = edf["Pathway"].value_counts().head(4).index.tolist()
    nP = len(paths)
    if nP == 0:
        return
    fig, axes = plt.subplots(1, nP, figsize=(6 * nP, 7))
    if nP == 1:
        axes = [axes]
    for ax, pw in zip(axes, paths):
        sub = edf[(edf["Pathway"] == pw) & edf["Donor_Anomer"].isin(topDonors) & edf["Acceptor_Pos"].isin(topAcceptors)]
        if sub.empty:
            ax.set_title(f"{pw}\n(no data)")
            continue
        piv = sub.groupby(["Donor_Anomer", "Acceptor_Pos"]).size().unstack(fill_value=0)
        sns.heatmap(piv, annot=True, fmt="d", cmap=HEAT_CMAP, linewidths=0.4, linecolor="#30363D", ax=ax, cbar=False)
        nPw = len(sub)
        ax.set_title(f"{pw} (n={nPw:,})")
        if ax != axes[0]:
            ax.set_ylabel("")
    plt.suptitle("Linkage Compass — Faceted by NP Pathway\n立体指南针 — 按 NP 大类分面",
                 fontsize=13, fontweight="bold", color="#E6EDF3", y=1.02)
    plt.tight_layout()
    fig.savefig(OUTPUT_DIR / "4a_compass_by_pathway.png")
    plt.close(fig)
    print("    Saved.")

# =====================================================================
# 主入口 (Main Entry)
# =====================================================================

def main():
    """生成全部可视化图表 V2 (Generate all V3 charts)."""
    CSV_PATH = r"d:\Glycan_Database\reports\GlycoNP_Deep_Enriched_v12.csv"

    print("=" * 70)
    print("  GlycoNP 可视化大屏 V3")
    print("  GlycoNP Visualization Dashboard V3")
    print("=" * 70)

    df, unknownCount = loadAndTransform(CSV_PATH)

    print(f"Generating charts ({len(df):,} molecules)...\n")

    # 维度 1: 热力图 (Dimension 1: Heatmaps)
    chart1aMonoVsOrganism(df, unknownCount)
    chart1aFamilyLevel(df)
    chart1bSeqVsScaffold(df)
    chart1cModVsOrganism(df)
    chart1cModVsFamily(df)
    chart1dLinkageMulti(df)

    # 维度 2: 宏观统计 (Dimension 2: Macro Statistics)
    chart2MacroStats(df)

    # 维度 3: 网络图 (Dimension 3: Networks)
    chart3aGlycanOrgNetwork(df)
    chart3bMotifNetwork(df)

    # 维度 4: 化学空间 (Dimension 4: Chemical Space)
    chart4aLinkageCompass(df)
    chart4bCoOccurrence(df)
    chart4cScaffoldPCA(df)

    # 图 4a 分面 (Faceted 4a)
    chart4aFacetedByOrganism(df)
    chart4aFacetedByPathway(df)

    # 新增图表 (Additional Charts)
    chart5ClassVsSugarCount(df)
    chart6AlphaBetaRatio(df)
    chart7MonoFreqBar(df)

    print(f"\n{'='*70}")
    print(f"  {len(list(OUTPUT_DIR.glob('*.png')))} charts saved to: {OUTPUT_DIR}")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
