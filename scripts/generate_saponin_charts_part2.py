#!/usr/bin/env python3
"""
Part 2: Remaining chart functions (A9-A17, B1-B3, C1-C3, D5-D7) and main().
This file is appended to the main generate_saponin_charts.py at runtime.
"""
# =========================================================================
# A9: Glycan-Organism Network (Interactive HTML via Plotly)
# =========================================================================

def chartA9GlycanOrgNetwork(df: pd.DataFrame):
    """Interactive network: top sugar sequences <-> organism types."""
    print("[A9] Glycan-Organism Network...")
    subDf = df[df["Sugar_Sequence"].notna() & df["Organism_Type"].notna()].copy()
    subDf = subDf[subDf["Organism_Type"] != "Unknown"]
    if subDf.empty:
        print("  [SKIP] No data"); return

    topSeqs = subDf["Sugar_Sequence"].value_counts().head(15).index.tolist()
    subDf = subDf[subDf["Sugar_Sequence"].isin(topSeqs)]

    # Build edge list and node positions using a simple circular layout
    import math
    links = subDf.groupby(["Sugar_Sequence", "Organism_Type"]).size().reset_index(name="Count")
    orgNames = links["Organism_Type"].unique().tolist()
    seqNames = links["Sugar_Sequence"].unique().tolist()
    allNodes = seqNames + orgNames
    nNodes = len(allNodes)
    posMap = {}
    for i, n in enumerate(allNodes):
        angle = 2 * math.pi * i / nNodes
        posMap[n] = (math.cos(angle), math.sin(angle))

    # Build edge traces
    edgeX, edgeY = [], []
    for _, row in links.iterrows():
        x0, y0 = posMap[row["Sugar_Sequence"]]
        x1, y1 = posMap[row["Organism_Type"]]
        edgeX += [x0, x1, None]
        edgeY += [y0, y1, None]

    edgeTrace = go.Scatter(x=edgeX, y=edgeY, mode="lines",
                           line=dict(width=0.5, color="#ccc"),
                           hoverinfo="none")

    # Truncate long labels for readability
    def truncLabel(s, maxLen=25):
        return str(s) if len(str(s)) <= maxLen else str(s)[:maxLen-3] + "..."

    nodeX = [posMap[n][0] for n in allNodes]
    nodeY = [posMap[n][1] for n in allNodes]
    nodeColors = ["#4C78A8"] * len(seqNames) + ["#E45756"] * len(orgNames)
    nodeLabels = [truncLabel(n) for n in allNodes]

    nodeTrace = go.Scatter(x=nodeX, y=nodeY, mode="markers+text",
                           marker=dict(size=14, color=nodeColors, line=dict(width=1, color="white")),
                           text=nodeLabels, textposition="top center",
                           textfont=dict(size=9),
                           hovertext=allNodes, hoverinfo="text")

    fig = go.Figure(data=[edgeTrace, nodeTrace])
    fig.update_layout(
        title="Glycan-Organism Network (Top 15 Sequences)",
        showlegend=False,
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
    )
    saveHtml(fig, "A9_glycan_organism_network.html")


# =========================================================================
# A10: Motif Conservation Network (similar to A9)
# =========================================================================

def chartA10MotifNetwork(df: pd.DataFrame):
    """Network: top disaccharide motifs <-> NP pathways."""
    print("[A10] Motif-Pathway Network...")
    rows = []
    for _, r in df.iterrows():
        pathway = r.get("NP_Pathway", "Unknown")
        if pathway == "Unknown":
            continue
        for motif in extractMotifs(r.get("Sugar_Sequence", ""), windowSize=2):
            rows.append({"Motif": motif, "Pathway": pathway})
    if not rows:
        print("  [SKIP] No motif data"); return

    mdf = pd.DataFrame(rows)
    topMotifs = mdf["Motif"].value_counts().head(15).index.tolist()
    mdf = mdf[mdf["Motif"].isin(topMotifs)]
    links = mdf.groupby(["Motif", "Pathway"]).size().reset_index(name="Count")
    links = links[links["Count"] >= 3]  # Filter weak links

    if links.empty:
        print("  [SKIP] No significant links"); return

    # Simple treemap as alternative visualization
    fig = px.treemap(links, path=["Pathway", "Motif"], values="Count",
                     color="Count", color_continuous_scale=HEATMAP_COLORSCALE)
    fig.update_layout(title="Disaccharide Motif × NP Pathway", height=650)
    saveHtml(fig, "A10_motif_pathway_treemap.html", height=650)


# =========================================================================
# A11/A12/A13: Linkage Compass Heatmaps
# =========================================================================

def chartA11LinkageCompass(df: pd.DataFrame):
    """Heatmap: donor+anomer x acceptor+position — the 'Linkage Compass'."""
    print("[A11-A13] Linkage Compass...")
    allEdges = []
    for _, r in df.iterrows():
        edges = parseSequenceToEdges(r.get("Sugar_Sequence", ""))
        for e in edges:
            e["Organism"] = r.get("Organism_Type", "Unknown")
            e["Pathway"] = r.get("NP_Pathway", "Unknown")
            allEdges.append(e)
    if not allEdges:
        print("  [SKIP] No linkage data"); return

    edf = pd.DataFrame(allEdges)
    edf["Donor_Label"] = edf["donor"] + " " + edf["anomer"]
    edf["Acceptor_Label"] = edf["acceptor"] + " @" + edf["pos"]

    # Overall compass (A11)
    topDonors = edf["Donor_Label"].value_counts().head(12).index.tolist()
    topAcc = edf["Acceptor_Label"].value_counts().head(12).index.tolist()
    sub = edf[edf["Donor_Label"].isin(topDonors) & edf["Acceptor_Label"].isin(topAcc)]
    piv = sub.groupby(["Donor_Label", "Acceptor_Label"]).size().unstack(fill_value=0)

    fig = px.imshow(piv.values, x=piv.columns.tolist(), y=piv.index.tolist(),
                    color_continuous_scale=HEATMAP_COLORSCALE, aspect="auto",
                    labels=dict(color="Count"), text_auto="d")
    fig.update_layout(
        title="Linkage Compass — Donor+Anomer × Acceptor+Position",
        xaxis_title="Acceptor + Position", yaxis_title="Donor + Anomer Type",
        height=650,
    )
    fig.update_xaxes(tickangle=45)
    saveHtml(fig, "A11_linkage_compass.html", height=650)

    # A12: Faceted by Organism type
    for orgType in ["Plant", "Animal"]:
        subOrg = edf[(edf["Organism"] == orgType) &
                     edf["Donor_Label"].isin(topDonors) &
                     edf["Acceptor_Label"].isin(topAcc)]
        if subOrg.empty:
            continue
        piv2 = subOrg.groupby(["Donor_Label", "Acceptor_Label"]).size().unstack(fill_value=0)
        fig2 = px.imshow(piv2.values, x=piv2.columns.tolist(), y=piv2.index.tolist(),
                         color_continuous_scale=HEATMAP_COLORSCALE, aspect="auto",
                         labels=dict(color="Count"), text_auto="d")
        fig2.update_layout(
            title=f"Linkage Compass — {orgType}",
            height=600,
        )
        fig2.update_xaxes(tickangle=45)
        saveHtml(fig2, f"A12_compass_{orgType.lower()}.html", height=600)

    # A13: Faceted by NP Pathway
    for pathway in ["Terpenoids"]:
        subPw = edf[(edf["Pathway"] == pathway) &
                    edf["Donor_Label"].isin(topDonors) &
                    edf["Acceptor_Label"].isin(topAcc)]
        if subPw.empty:
            continue
        piv3 = subPw.groupby(["Donor_Label", "Acceptor_Label"]).size().unstack(fill_value=0)
        fig3 = px.imshow(piv3.values, x=piv3.columns.tolist(), y=piv3.index.tolist(),
                         color_continuous_scale=HEATMAP_COLORSCALE, aspect="auto",
                         labels=dict(color="Count"), text_auto="d")
        fig3.update_layout(title=f"Linkage Compass — {pathway}", height=600)
        fig3.update_xaxes(tickangle=45)
        saveHtml(fig3, f"A13_compass_{pathway.lower()}.html", height=600)


# =========================================================================
# A14: Monosaccharide × Modification Co‑occurrence
# =========================================================================

def chartA14CoOccurrence(df: pd.DataFrame):
    """Heatmap: sugar type vs modification co-occurrence."""
    print("[A14] Sugar × Modification Co-occurrence...")
    rows = []
    for _, r in df.iterrows():
        for sugar in set(r.get("Mono_List", [])):
            for mod in r.get("Mod_Tags", []):
                rows.append({"Sugar": sugar, "Modification": mod})
    if not rows:
        print("  [SKIP] No data"); return
    mdf = pd.DataFrame(rows)
    topSugars = mdf["Sugar"].value_counts().head(12).index.tolist()
    topMods = mdf["Modification"].value_counts().head(10).index.tolist()
    mdf = mdf[mdf["Sugar"].isin(topSugars) & mdf["Modification"].isin(topMods)]
    piv = mdf.groupby(["Sugar", "Modification"]).size().unstack(fill_value=0)

    fig = px.imshow(piv.values, x=piv.columns.tolist(), y=piv.index.tolist(),
                    color_continuous_scale=HEATMAP_COLORSCALE, aspect="auto",
                    labels=dict(color="Count"), text_auto="d")
    fig.update_layout(
        title="Monosaccharide × Modification Co-occurrence",
        xaxis_title="Modification", yaxis_title="Monosaccharide",
    )
    saveHtml(fig, "A14_sugar_mod_cooccurrence.html")


# =========================================================================
# A16: NP Class × Sugar Count Heatmap
# =========================================================================

def chartA16ClassVsSugarCount(df: pd.DataFrame):
    """Heatmap: NP class vs binned sugar count."""
    print("[A16] NP Class × Sugar Count...")
    subDf = df[df["NP_Class"] != "Unknown"].copy()
    topClasses = subDf["NP_Class"].value_counts().head(12).index.tolist()
    subDf = subDf[subDf["NP_Class"].isin(topClasses)]
    subDf["Sugar_Bin"] = subDf["Total_Sugar_Count"].clip(upper=7).astype(int).astype(str)
    piv = subDf.groupby(["NP_Class", "Sugar_Bin"]).size().unstack(fill_value=0)
    piv = piv.reindex(columns=sorted(piv.columns, key=int), fill_value=0)

    # Add N= to row labels so users see sample size per class
    rowTotals = piv.sum(axis=1)
    yLabels = [f"{cls} (N={rowTotals[cls]:,})" for cls in piv.index]

    fig = px.imshow(piv.values, x=piv.columns.tolist(), y=yLabels,
                    color_continuous_scale=HEATMAP_COLORSCALE, aspect="auto",
                    labels=dict(color="Count"), text_auto="d")
    fig.update_layout(
        title="NP Class × Sugar Count Distribution",
        xaxis_title="Total Sugar Count", yaxis_title="NP Class",
    )
    saveHtml(fig, "A16_class_vs_sugar_count.html")


# =========================================================================
# A17: Alpha/Beta Ratio by NP Superclass
# =========================================================================

def chartA17AlphaBetaRatio(df: pd.DataFrame):
    """Grouped bar: alpha vs beta glycosidic bond ratio per NP superclass."""
    print("[A17] Alpha/Beta Ratio...")
    bondRows = []
    for _, r in df.iterrows():
        bondRaw = r.get("Glycan-Aglycone_Bond_Detail", "[]")
        try:
            bonds = json.loads(str(bondRaw)) if pd.notna(bondRaw) and str(bondRaw) not in ("nan","","[]") else []
        except Exception:
            bonds = []
        sc = r.get("NP_Superclass", "Unknown")
        for bd in bonds:
            bond = bd.get("bond", "")
            if "α" in bond or "alpha" in bond.lower():
                bondRows.append({"Superclass": sc, "Anomer": "α (alpha)"})
            elif "β" in bond or "beta" in bond.lower():
                bondRows.append({"Superclass": sc, "Anomer": "β (beta)"})

    if not bondRows:
        print("  [SKIP] No bond data"); return
    bdf = pd.DataFrame(bondRows)
    bdf = bdf[bdf["Superclass"] != "Unknown"]
    # Only saponin-relevant superclasses
    validSc = ["Steroids", "Triterpenoids"]
    bdf = bdf[bdf["Superclass"].isin(validSc)]
    piv = bdf.groupby(["Superclass", "Anomer"]).size().unstack(fill_value=0)

    # Add N= to x-axis labels
    scTotals = bdf["Superclass"].value_counts()
    xLabels = [f"{sc}\n(N={scTotals.get(sc, 0):,})" for sc in piv.index]

    fig = go.Figure()
    for col in piv.columns:
        fig.add_trace(go.Bar(name=col, x=xLabels, y=piv[col].values,
                             text=piv[col].values, textposition="outside"))
    fig.update_layout(
        title="α vs β Glycosidic Bond Count (Steroids vs Triterpenoids)",
        xaxis_title="NP Superclass", yaxis_title="Count",
        barmode="group", height=500,
    )
    saveHtml(fig, "A17_alpha_beta_ratio.html", height=500)


# =========================================================================
# B1: Sankey Diagram — Superclass → Modification → Sugar
# =========================================================================

def chartB1Sankey(df: pd.DataFrame):
    """Sankey flow: NP Superclass -> Modification Tag -> Sugar Name."""
    print("[B1] Sankey Diagram...")
    rows = []
    for _, r in df.iterrows():
        sc = r.get("NP_Superclass", "Unknown")
        if sc == "Unknown":
            continue
        sugars = r.get("Mono_List", [])[:1]  # Take first sugar per molecule
        mods = r.get("Mod_Tags", [])
        if not mods:
            mods = ["None"]
        for mod in mods[:1]:
            for sugar in sugars:
                rows.append({"Superclass": sc, "Mod": mod, "Sugar": sugar})
    if not rows:
        print("  [SKIP] No data"); return

    sdf = pd.DataFrame(rows)
    # Filter to only saponin-relevant superclasses
    validSc = ["Steroids", "Triterpenoids"]
    sdf = sdf[sdf["Superclass"].isin(validSc)]
    topSc = sdf["Superclass"].value_counts().index.tolist()
    topMods = sdf["Mod"].value_counts().head(6).index.tolist()
    topSugars = sdf["Sugar"].value_counts().head(8).index.tolist()
    sdf = sdf[sdf["Superclass"].isin(topSc) & sdf["Mod"].isin(topMods) & sdf["Sugar"].isin(topSugars)]

    # Build labels with N= where useful
    scCounts = sdf["Superclass"].value_counts()
    labelsSc = [f"{s} (N={scCounts.get(s,0):,})" for s in topSc]
    allLabels = labelsSc + topMods + topSugars
    labelIdx = {}
    for i, s in enumerate(topSc):
        labelIdx[s] = i
    for i, m in enumerate(topMods):
        labelIdx[m] = len(topSc) + i
    for i, s in enumerate(topSugars):
        labelIdx[s] = len(topSc) + len(topMods) + i

    # Build links: Superclass -> Mod, Mod -> Sugar
    link1 = sdf.groupby(["Superclass", "Mod"]).size().reset_index(name="Count")
    link2 = sdf.groupby(["Mod", "Sugar"]).size().reset_index(name="Count")

    source = [labelIdx[r["Superclass"]] for _, r in link1.iterrows()] + \
             [labelIdx[r["Mod"]] for _, r in link2.iterrows()]
    target = [labelIdx[r["Mod"]] for _, r in link1.iterrows()] + \
             [labelIdx[r["Sugar"]] for _, r in link2.iterrows()]
    value = link1["Count"].tolist() + link2["Count"].tolist()

    nodeColors = ["#4C78A8"] * len(topSc) + ["#72B7B2"] * len(topMods) + ["#E45756"] * len(topSugars)

    fig = go.Figure(go.Sankey(
        node=dict(pad=20, thickness=20, label=allLabels, color=nodeColors),
        link=dict(source=source, target=target, value=value, color="rgba(200,200,200,0.4)"),
    ))
    fig.update_layout(title="Sankey: NP Superclass → Modification → Sugar", height=650)
    saveHtml(fig, "B1_sankey_flow.html", height=650)


# =========================================================================
# B2: Sugar Sequence × NP Superclass Heatmap (Gold Pair)
# =========================================================================

def chartB2SequenceHeatmap(df: pd.DataFrame):
    """Heatmap: top sugar sequences vs top NP superclasses."""
    print("[B2] Sugar Sequence × Superclass Heatmap...")
    subDf = df[(df["Sugar_Sequence"].notna()) & (df["NP_Superclass"] != "Unknown")]
    # Only saponin-relevant superclasses
    validSc = ["Steroids", "Triterpenoids"]
    subDf = subDf[subDf["NP_Superclass"].isin(validSc)]
    topSeqs = subDf["Sugar_Sequence"].value_counts().head(15).index.tolist()
    subDf = subDf[subDf["Sugar_Sequence"].isin(topSeqs)]

    # Truncate labels for readability
    def trunc(s, n=30):
        return str(s) if len(str(s)) <= n else str(s)[:n-3] + "..."
    subDf = subDf.copy()
    subDf["Seq_Short"] = subDf["Sugar_Sequence"].apply(trunc)
    piv = subDf.groupby(["Seq_Short", "NP_Superclass"]).size().unstack(fill_value=0)
    colTotals = piv.sum(axis=0)
    piv_log = np.log2(piv + 1)
    # Add N= to column labels
    piv_log.columns = [f"{c} (N={int(colTotals[c]):,})" for c in piv_log.columns]

    fig = px.imshow(piv_log.values, x=piv_log.columns.tolist(), y=piv_log.index.tolist(),
                    color_continuous_scale=HEATMAP_COLORSCALE, aspect="auto",
                    labels=dict(color="log₂(Count+1)"), text_auto=".1f")
    fig.update_layout(
        title="Sugar Sequence × NP Superclass (Steroids vs Triterpenoids, log₂ scale)",
        xaxis_title="NP Superclass", yaxis_title="Sugar Sequence",
        height=700,
    )
    saveHtml(fig, "B2_sequence_superclass_heatmap.html", height=700)


# =========================================================================
# B3: Kingdom × Modification Grouped Bar
# =========================================================================

def chartB3KingdomMods(df: pd.DataFrame):
    """Grouped bar: modification rates across biological kingdoms."""
    print("[B3] Kingdom × Modification Fingerprint...")
    rows = []
    for _, r in df.iterrows():
        kingdom = str(r.get("LOTUS_kingdom", ""))
        if not kingdom or kingdom == "nan":
            continue
        mods = r.get("Mod_Tags", [])
        for mod in mods:
            rows.append({"Kingdom": kingdom, "Mod": mod})
    if not rows:
        print("  [SKIP] No data"); return

    mdf = pd.DataFrame(rows)
    topK = mdf["Kingdom"].value_counts().head(4).index.tolist()
    topM = mdf["Mod"].value_counts().head(8).index.tolist()
    mdf = mdf[mdf["Kingdom"].isin(topK) & mdf["Mod"].isin(topM)]
    piv = mdf.groupby(["Kingdom", "Mod"]).size().unstack(fill_value=0)

    fig = go.Figure()
    for mod in piv.columns:
        fig.add_trace(go.Bar(name=mod, x=piv.index.tolist(), y=piv[mod].values))
    fig.update_layout(
        title="Modification Fingerprint Across Biological Kingdoms",
        xaxis_title="Kingdom", yaxis_title="Count",
        barmode="group", height=550,
    )
    saveHtml(fig, "B3_kingdom_mods.html", height=550)


# =========================================================================
# D5: Literature Coverage
# =========================================================================

def chartD5LiteratureCoverage(df: pd.DataFrame):
    """Donut chart and bar showing DOI coverage."""
    print("[D5] Literature Coverage...")
    hasDoi = df["dois"].notna() & (df["dois"] != "")
    counts = hasDoi.value_counts()

    fig = make_subplots(rows=1, cols=2,
        specs=[[{"type": "pie"}, {"type": "bar"}]],
        subplot_titles=["DOI Coverage", "Top Journals"])

    fig.add_trace(go.Pie(
        labels=["Has DOI", "No DOI"],
        values=[counts.get(True, 0), counts.get(False, 0)],
        hole=0.4, marker_colors=["#4C78A8", "#ddd"],
    ), row=1, col=1)

    # Extract journal names from DOIs (approximate)
    allDois = df["dois"].dropna().str.cat(sep=" | ")
    journalPattern = r"10\.\d+/([A-Z]+)"
    journals = re.findall(journalPattern, allDois, re.IGNORECASE)
    journalCounts = Counter(journals).most_common(10)
    if journalCounts:
        jNames, jCounts = zip(*journalCounts)
        fig.add_trace(go.Bar(x=list(jNames), y=list(jCounts), marker_color="#E45756"), row=1, col=2)

    fig.update_layout(title="Literature Coverage Overview", height=450, showlegend=False)
    saveHtml(fig, "D5_literature_coverage.html", height=450)


# =========================================================================
# D6: Bioactivity Overview
# =========================================================================

def chartD6Bioactivity(df: pd.DataFrame):
    """Donut chart and bar for bioactivity data coverage."""
    print("[D6] Bioactivity Overview...")
    hasBio = df["bioactivity_summary"].notna() & (df["bioactivity_summary"] != "")
    counts = hasBio.value_counts()

    fig = make_subplots(rows=1, cols=2,
        specs=[[{"type": "pie"}, {"type": "bar"}]],
        subplot_titles=["Bioactivity Coverage", "Top Targets"])

    fig.add_trace(go.Pie(
        labels=["Has Bioactivity", "No Data"],
        values=[counts.get(True, 0), counts.get(False, 0)],
        hole=0.4, marker_colors=["#72B7B2", "#ddd"],
    ), row=1, col=1)

    hasTargets = df["ChEMBL_Targets"].dropna()
    hasTargets = hasTargets[hasTargets != ""]
    if not hasTargets.empty:
        allTargets = hasTargets.str.cat(sep=" | ").split(" | ")
        tgtCounts = Counter([t.strip() for t in allTargets if t.strip()]).most_common(10)
        if tgtCounts:
            tNames, tCounts = zip(*tgtCounts)
            fig.add_trace(go.Bar(x=list(tNames), y=list(tCounts),
                                 marker_color="#B279A2"), row=1, col=2)

    fig.update_layout(title="Bioactivity Data Overview", height=450, showlegend=False)
    saveHtml(fig, "D6_bioactivity_overview.html", height=450)


# =========================================================================
# D7: TPSA vs MW Scatter (alternative chemical space view)
# =========================================================================

def chartD7TpsaVsMw(df: pd.DataFrame):
    """TPSA vs MW scatter, colored by Total_Sugar_Count."""
    print("[D7] TPSA vs MW Chemical Space...")
    subDf = df[df["topological_polar_surface_area"].notna() &
               df["exact_molecular_weight"].notna()].copy()
    subDf["Sugar_Count_Label"] = subDf["Total_Sugar_Count"].clip(upper=6).astype(int).astype(str)
    subDf.loc[subDf["Total_Sugar_Count"] > 6, "Sugar_Count_Label"] = "7+"

    fig = px.scatter(subDf, x="exact_molecular_weight",
                     y="topological_polar_surface_area",
                     color="Sugar_Count_Label", opacity=0.35,
                     color_discrete_sequence=px.colors.sequential.Viridis,
                     labels={"exact_molecular_weight": "Molecular Weight (Da)",
                             "topological_polar_surface_area": "TPSA (Å²)",
                             "Sugar_Count_Label": "Sugar Count"},
                     hover_data=["name", "Saponin_Type"])
    fig.update_layout(
        title=f"Chemical Space: TPSA vs MW (N={len(subDf):,})",
        height=600,
    )
    saveHtml(fig, "D7_tpsa_vs_mw.html", height=600)


# =========================================================================
# A15: Scaffold PCA (using Plotly scatter)
# =========================================================================

def chartA15ScaffoldPCA(df: pd.DataFrame):
    """PCA of molecular descriptors, colored by Saponin type."""
    print("[A15] Scaffold PCA projection...")
    try:
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
    except ImportError:
        print("  [SKIP] sklearn not available"); return

    featureCols = ["exact_molecular_weight", "alogp", "topological_polar_surface_area",
                   "Total_Sugar_Count", "Max_Chain_Length"]
    subDf = df.dropna(subset=featureCols).copy()
    if len(subDf) < 50:
        print("  [SKIP] Not enough data"); return

    X = subDf[featureCols].values
    scaler = StandardScaler()
    Xs = scaler.fit_transform(X)
    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(Xs)
    subDf["PC1"] = coords[:, 0]
    subDf["PC2"] = coords[:, 1]
    explained = pca.explained_variance_ratio_ * 100

    fig = px.scatter(subDf, x="PC1", y="PC2", color="Saponin_Type", opacity=0.4,
                     color_discrete_map={"Steroidal": "#4C78A8",
                                         "Triterpenoid": "#E45756", "Other": "#72B7B2"},
                     labels={"PC1": f"PC1 ({explained[0]:.1f}%)",
                             "PC2": f"PC2 ({explained[1]:.1f}%)"},
                     hover_data=["name", "Detailed_NP_Class"])
    fig.update_layout(
        title=f"PCA of Physicochemical Descriptors (N={len(subDf):,})",
        height=600,
    )
    saveHtml(fig, "A15_scaffold_pca.html", height=600)


# =========================================================================
# C1: Summary Markdown Report
# =========================================================================

def reportC1Summary(df: pd.DataFrame):
    """Generate a Markdown summary report of glycan frequency stats."""
    print("[C1] Summary Markdown Report...")
    lines = ["# GlycoNP Saponin Database — Statistical Summary\n"]
    lines.append(f"- **Total molecules**: {len(df):,}")
    lines.append(f"- **Steroidal saponins**: {(df['Saponin_Type']=='Steroidal').sum():,}")
    lines.append(f"- **Triterpenoid saponins**: {(df['Saponin_Type']=='Triterpenoid').sum():,}")
    lines.append(f"- **Avg sugar count**: {df['Total_Sugar_Count'].mean():.2f}")
    lines.append(f"- **Avg MW**: {df['exact_molecular_weight'].mean():.1f} Da")
    lines.append(f"- **Max chain length (avg)**: {df['Max_Chain_Length'].mean():.2f}\n")

    # Top sugars
    allSugars = []
    for lst in df["Mono_List"]:
        if isinstance(lst, list):
            allSugars.extend(lst)
    sc = Counter(allSugars).most_common(10)
    lines.append("## Top 10 Monosaccharides\n")
    lines.append("| Rank | Sugar | Count |")
    lines.append("|------|-------|-------|")
    for i, (name, cnt) in enumerate(sc, 1):
        lines.append(f"| {i} | {name} | {cnt:,} |")

    outPath = OUT_DIR / "C1_summary_report.md"
    outPath.write_text("\n".join(lines), encoding="utf-8")
    print(f"  -> Saved: {outPath.name}")


# =========================================================================
# Part 2 函数库 — main() 已移至 generate_saponin_charts.py 主脚本
# Part 2 function library — main() moved to main script
# =========================================================================

