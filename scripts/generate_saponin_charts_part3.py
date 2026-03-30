#!/usr/bin/env python3
"""
==========================================================================
  GlycoNP Saponin Charts — Part 3: Expanded Sugar Chain Analysis
  皂苷图表扩展包: 寡糖序列模式 / 糖链拓扑 / 多糖片段匹配 / 稀有糖 / 共现

  新增 14 张图 (S1-S14), exec 方式并入 generate_saponin_charts.py

  S1  二糖 (Disaccharide) 序列频率 Top-15
  S2  三糖 (Trisaccharide) 序列频率 Top-12
  S3  四糖+ (Tetrasaccharide+) 序列频率 Top-10
  S4  非还原端糖 (Terminal Sugar) 分布
  S5  还原端糖 (Reducing-End Sugar) — 直连苷元的糖
  S6  分支点糖 (Branch-Point Sugar) — 催生支链的糖
  S7  单链 vs 双链皂苷 (Mono- vs Bidesmosidic) 饼图
  S8  糖苷键位置偏好 (Linkage Position Preference) 雷达图
  S9  稀有糖/脱氧糖分布 (Rare & Deoxy Sugar)
  S10 糖组合共现网络 (Sugar Co-occurrence Network)
  S11 多糖片段匹配 (Polysaccharide Motif Matching)
  S12 糖多样性指数 (Sugar Diversity Index) 分布
  S13 糖链长度 × 皂苷细类 热力图
  S14 各 Detailed_NP_Class 的糖指纹雷达图
==========================================================================

[TEST DATA ONLY]
"""

# ── 注意: 本文件通过 exec() 并入 generate_saponin_charts.py ──────────
# 所有全局变量 (df, OUT_DIR, saveHtml, LAYOUT_DEFAULTS, HEATMAP_COLORSCALE,
# QUAL_COLORS, extractMonosaccharideList, parseSequenceToEdges 等)
# 均由主脚本提供。


# ═════════════════════════════════════════════════════════════════════════
# 工具函数 (Utilities for Part 3)
# ═════════════════════════════════════════════════════════════════════════

def _extractChains(seq: str):
    """将序列拆分为独立糖链 (分号分隔)。
    Split sequence into independent sugar chains (semicolon-separated).
    返回: List[str], 每条链为一个子字符串
    """
    if pd.isna(seq) or not seq:
        return []
    return [c.strip() for c in str(seq).split(";") if c.strip()
            and c.strip() not in ("Non_Cyclic_Invalid", "Invalid", "Error")]


def _countSugarsInChain(chain: str) -> int:
    """统计单条链中的糖个数。"""
    pattern = r'(Neu5Ac|Neu5Gc|KDO|Kdn|[DL]-\w+|Hex|dHex|Pen|HexA|HexN|HexNAc|Hept|Oct|Non)'
    return len(re.findall(pattern, chain))


def _extractOrderedSugars(chain: str) -> list:
    """按顺序提取链中的糖名 (保留重复)。
    Extract sugar names in order (keep duplicates).
    E.g. 'D-Glc-(b1-4)-L-Rha-(a1-2)-D-GlcA' -> ['D-Glc', 'L-Rha', 'D-GlcA']
    """
    # 先去掉分支标记 (Remove branch markers)
    clean = re.sub(r'\[', '', re.sub(r'\]-', ' ', str(chain)))
    pattern = r'(Neu5Ac|Neu5Gc|KDO|Kdn|[DL]-(?:d|me|6d)?[A-Z][a-z]+[A-Za-z]*|Hex|dHex|Pen|HexA|HexN|HexNAc|Hept|Oct|Non)'
    return re.findall(pattern, clean)


def _getTerminalSugar(chain: str) -> str:
    """提取非还原端糖 (链的第一个糖)。
    Get non-reducing end sugar (first sugar in chain).
    """
    sugars = _extractOrderedSugars(chain)
    return sugars[0] if sugars else ""


def _getReducingEndSugar(chain: str) -> str:
    """提取还原端糖 (链的最后一个, 直连苷元)。
    Get reducing-end sugar (last sugar, directly attached to aglycone).
    """
    sugars = _extractOrderedSugars(chain)
    return sugars[-1] if sugars else ""


def _getBranchPointSugars(seq: str) -> list:
    """提取分支点糖——即 [X]-Y 结构中 Y 后面的糖。
    Extract branch-point sugars: the sugar that bears a branch.
    Pattern: [branch]-MainSugar -> MainSugar is the branch point.
    """
    if pd.isna(seq) or "[" not in str(seq):
        return []
    # 匹配 ]-Sugar 模式 (Match ]-Sugar pattern)
    pattern = r'\]-(?:\([ab]\d-\d\)-)?((?:Neu5Ac|Neu5Gc|KDO|[DL]-\w+|Hex|dHex|Pen|HexA))'
    return re.findall(pattern, str(seq))


def _buildFullSequenceLabel(chain: str, maxLen: int = 50) -> str:
    """构建带截断的序列标签。"""
    s = str(chain).strip()
    return s if len(s) <= maxLen else s[:maxLen-3] + "..."


# ═════════════════════════════════════════════════════════════════════════
# S1: 二糖序列频率 (Disaccharide Sequence Frequency)
# ═════════════════════════════════════════════════════════════════════════

def chartS1DisaccharideFreq(df: pd.DataFrame):
    """Top-15 most common disaccharide motifs.
    在皂苷数据库中，二糖是最基本的糖链构建单元。
    常见如 D-Glc-(b1-2)-D-GlcA（三萜皂苷的标志性起始二糖）。
    """
    print("[S1] Disaccharide Sequence Frequency...")
    motifs = []
    for seq in df["Sugar_Sequence"].dropna():
        for chain in _extractChains(seq):
            sugars = _extractOrderedSugars(chain)
            if len(sugars) == 2:
                # 提取完整二糖表达式 (Extract full disaccharide expression)
                motifs.append(_buildFullSequenceLabel(chain, 45))
            elif len(sugars) > 2:
                # 滑窗提取所有二糖子序列 (Sliding window: all disaccharide sub-sequences)
                fragments = re.split(r'-\([ab]\d-\d\)-', chain)
                linkages = re.findall(r'\([ab]\d-\d\)', chain)
                for i in range(len(fragments) - 1):
                    if i < len(linkages):
                        motif = f"{fragments[i].strip()}-{linkages[i]}-{fragments[i+1].strip()}"
                        # 清理分支标记 (Clean branch markers)
                        motif = re.sub(r'[\[\]]', '', motif).strip()
                        if motif:
                            motifs.append(_buildFullSequenceLabel(motif, 45))

    if not motifs:
        print("  [SKIP] No disaccharide data"); return

    mc = Counter(motifs).most_common(15)
    names, counts = zip(*mc)
    fig = go.Figure(go.Bar(
        y=list(reversed(names)), x=list(reversed(counts)),
        orientation="h", marker_color="#4C78A8",
        text=list(reversed(counts)), textposition="outside"))
    fig.update_layout(
        title=f"Top 15 Disaccharide Motifs in Saponins (N={sum(counts):,})",
        xaxis_title="Occurrence", yaxis_title="Disaccharide Motif",
        height=600, margin=dict(l=280))
    saveHtml(fig, "S1_disaccharide_frequency.html", height=600, width=1200)


# ═════════════════════════════════════════════════════════════════════════
# S2: 三糖序列频率 (Trisaccharide Sequence Frequency)
# ═════════════════════════════════════════════════════════════════════════

def chartS2TrisaccharideFreq(df: pd.DataFrame):
    """Top-12 most common trisaccharide sequences.
    三糖序列揭示糖基转移酶的连续底物偏好性。
    """
    print("[S2] Trisaccharide Sequence Frequency...")
    seqs = []
    for seq in df["Sugar_Sequence"].dropna():
        for chain in _extractChains(seq):
            nSugars = _countSugarsInChain(chain)
            if nSugars == 3:
                seqs.append(_buildFullSequenceLabel(chain, 55))
    if not seqs:
        print("  [SKIP] No trisaccharide data"); return

    mc = Counter(seqs).most_common(12)
    names, counts = zip(*mc)
    fig = go.Figure(go.Bar(
        y=list(reversed(names)), x=list(reversed(counts)),
        orientation="h", marker_color="#E45756",
        text=list(reversed(counts)), textposition="outside"))
    fig.update_layout(
        title=f"Top 12 Trisaccharide Sequences (N={sum(counts):,})",
        xaxis_title="Occurrence", yaxis_title="Trisaccharide",
        height=550, margin=dict(l=350))
    saveHtml(fig, "S2_trisaccharide_frequency.html", height=550, width=1300)


# ═════════════════════════════════════════════════════════════════════════
# S3: 四糖+序列频率 (Tetrasaccharide+ Frequency)
# ═════════════════════════════════════════════════════════════════════════

def chartS3TetrasaccharideFreq(df: pd.DataFrame):
    """Top-10 most common ≥4-sugar chains.
    四糖及以上序列在皂苷中较罕见，揭示复杂糖基化的生物合成路径。
    """
    print("[S3] Tetrasaccharide+ Sequence Frequency...")
    seqs = []
    for seq in df["Sugar_Sequence"].dropna():
        for chain in _extractChains(seq):
            nSugars = _countSugarsInChain(chain)
            if nSugars >= 4:
                seqs.append(_buildFullSequenceLabel(chain, 70))
    if not seqs:
        print("  [SKIP] No tetrasaccharide+ data"); return

    mc = Counter(seqs).most_common(10)
    names, counts = zip(*mc)
    fig = go.Figure(go.Bar(
        y=list(reversed(names)), x=list(reversed(counts)),
        orientation="h", marker_color="#72B7B2",
        text=list(reversed(counts)), textposition="outside"))
    fig.update_layout(
        title=f"Top 10 Tetrasaccharide+ Sequences (≥4 sugars, N={sum(counts):,})",
        xaxis_title="Occurrence", yaxis_title="Chain Sequence",
        height=500, margin=dict(l=420))
    saveHtml(fig, "S3_tetrasaccharide_frequency.html", height=500, width=1400)


# ═════════════════════════════════════════════════════════════════════════
# S4: 非还原端糖分布 (Terminal / Non-reducing End Sugar)
# ═════════════════════════════════════════════════════════════════════════

def chartS4TerminalSugar(df: pd.DataFrame):
    """Non-reducing end sugar frequency, split by saponin type.
    非还原端糖是糖链延伸的最后一步，反映最后一个糖基转移酶的底物。
    """
    print("[S4] Terminal (Non-reducing End) Sugar Distribution...")
    rows = []
    for _, r in df.iterrows():
        sapType = r.get("Saponin_Type", "Other")
        for chain in _extractChains(r.get("Sugar_Sequence", "")):
            terminal = _getTerminalSugar(chain)
            if terminal:
                rows.append({"Terminal": terminal, "Type": sapType})
    if not rows: print("  [SKIP]"); return
    tdf = pd.DataFrame(rows)
    tdf = tdf[tdf["Type"].isin(["Steroidal", "Triterpenoid"])]
    topT = tdf["Terminal"].value_counts().head(12).index
    tdf = tdf[tdf["Terminal"].isin(topT)]
    piv = tdf.groupby(["Terminal", "Type"]).size().unstack(fill_value=0)

    fig = go.Figure()
    colors = {"Steroidal": "#4C78A8", "Triterpenoid": "#E45756"}
    for col in piv.columns:
        fig.add_trace(go.Bar(name=col, x=piv.index.tolist(), y=piv[col],
                             marker_color=colors.get(col, "#999"),
                             text=piv[col], textposition="outside"))
    fig.update_layout(
        title="Non-reducing End Sugar — Steroidal vs Triterpenoid",
        xaxis_title="Terminal Sugar", yaxis_title="Count",
        barmode="group", height=550)
    saveHtml(fig, "S4_terminal_sugar.html", height=550)


# ═════════════════════════════════════════════════════════════════════════
# S5: 还原端糖 (Reducing-End Sugar) — 直连苷元
# ═════════════════════════════════════════════════════════════════════════

def chartS5ReducingEndSugar(df: pd.DataFrame):
    """Reducing-end sugar (directly attached to aglycone), split by type.
    还原端糖直接与苷元骨架相连，反映初始糖基化酶的底物偏好。
    例如三萜皂苷几乎总是 D-GlcA 连 C-3。
    """
    print("[S5] Reducing-End Sugar Distribution...")
    rows = []
    for _, r in df.iterrows():
        sapType = r.get("Saponin_Type", "Other")
        for chain in _extractChains(r.get("Sugar_Sequence", "")):
            red = _getReducingEndSugar(chain)
            if red:
                rows.append({"ReducEnd": red, "Type": sapType})
    if not rows: print("  [SKIP]"); return
    rdf = pd.DataFrame(rows)
    rdf = rdf[rdf["Type"].isin(["Steroidal", "Triterpenoid"])]
    topR = rdf["ReducEnd"].value_counts().head(12).index
    rdf = rdf[rdf["ReducEnd"].isin(topR)]
    piv = rdf.groupby(["ReducEnd", "Type"]).size().unstack(fill_value=0)

    fig = go.Figure()
    colors = {"Steroidal": "#4C78A8", "Triterpenoid": "#E45756"}
    for col in piv.columns:
        fig.add_trace(go.Bar(name=col, x=piv.index.tolist(), y=piv[col],
                             marker_color=colors.get(col, "#999"),
                             text=piv[col], textposition="outside"))
    fig.update_layout(
        title="Reducing-End Sugar (Directly Attached to Aglycone)<br>"
              "<sub>Reflects initial glycosyltransferase substrate preference</sub>",
        xaxis_title="Reducing-End Sugar", yaxis_title="Count",
        barmode="group", height=550)
    saveHtml(fig, "S5_reducing_end_sugar.html", height=550)


# ═════════════════════════════════════════════════════════════════════════
# S6: 分支点糖 (Branch-Point Sugar)
# ═════════════════════════════════════════════════════════════════════════

def chartS6BranchPointSugar(df: pd.DataFrame):
    """Branch-point sugar frequency — which sugars carry branching.
    分支点糖承载支链——揭示哪些糖具有多连接点糖基化能力。
    """
    print("[S6] Branch-Point Sugar Distribution...")
    allBranch = []
    for seq in df["Sugar_Sequence"].dropna():
        allBranch.extend(_getBranchPointSugars(seq))
    if not allBranch:
        print("  [SKIP] No branched structures"); return

    bc = Counter(allBranch).most_common(10)
    names, counts = zip(*bc)
    fig = go.Figure(go.Bar(
        x=list(names), y=list(counts), marker_color="#B279A2",
        text=list(counts), textposition="outside"))
    fig.update_layout(
        title=f"Branch-Point Sugar Frequency (N={sum(counts):,} branches total)<br>"
              "<sub>Sugars that carry a side-chain attachment in branched saponin glycans</sub>",
        xaxis_title="Branch-Point Sugar", yaxis_title="Count", height=500)
    saveHtml(fig, "S6_branch_point_sugar.html", height=500)


# ═════════════════════════════════════════════════════════════════════════
# S7: 单链 vs 双链皂苷 (Mono- vs Bidesmosidic)
# ═════════════════════════════════════════════════════════════════════════

def chartS7MonoBiDesmosidic(df: pd.DataFrame):
    """Monodesmosidic vs bidesmosidic saponin distribution, by type.
    单链皂苷 (1 条糖链, 通常连 C-3) vs 双链皂苷 (2 条糖链, C-3 + C-28/C-26)。
    """
    print("[S7] Monodesmosidic vs Bidesmosidic...")
    rows = []
    for _, r in df.iterrows():
        chains = _extractChains(r.get("Sugar_Sequence", ""))
        nChains = len(chains)
        sapType = r.get("Saponin_Type", "Other")
        if nChains == 0: continue
        desmo = "Mono" if nChains == 1 else ("Bi" if nChains == 2 else f"Tri+({nChains})")
        rows.append({"Desmosidic": desmo, "Type": sapType})
    if not rows: print("  [SKIP]"); return
    ddf = pd.DataFrame(rows)
    ddf = ddf[ddf["Type"].isin(["Steroidal", "Triterpenoid"])]

    fig = px.histogram(ddf, x="Type", color="Desmosidic", barmode="group",
                       color_discrete_map={"Mono": "#4C78A8", "Bi": "#E45756", "Tri+(3)": "#72B7B2"},
                       text_auto=True)
    fig.update_layout(
        title="Monodesmosidic vs Bidesmosidic Saponins<br>"
              "<sub>Mono=1 sugar chain, Bi=2 chains (e.g. C-3 + C-28), Tri+=3+ chains</sub>",
        xaxis_title="Saponin Type", yaxis_title="Count", height=500)
    saveHtml(fig, "S7_mono_bi_desmosidic.html", height=500)


# ═════════════════════════════════════════════════════════════════════════
# S8: 糖苷键位置偏好雷达图 (Linkage Position Preference Radar)
# ═════════════════════════════════════════════════════════════════════════

def chartS8LinkagePositionRadar(df: pd.DataFrame):
    """Per-sugar linkage position preference radar.
    每种糖偏好哪个位置？L-Rha 多见 (1→2), D-Glc 多见 (1→4)。
    """
    print("[S8] Linkage Position Preference Radar...")
    allEdges = []
    for seq in df["Sugar_Sequence"].dropna():
        allEdges.extend(parseSequenceToEdges(seq))
    if not allEdges: print("  [SKIP]"); return
    edf = pd.DataFrame(allEdges)
    topDonors = edf["donor"].value_counts().head(6).index.tolist()
    edf = edf[edf["donor"].isin(topDonors)]

    # 位置类别 (Position categories)
    positions = ["1→2", "1→3", "1→4", "1→6", "Other"]
    fig = go.Figure()
    for i, sugar in enumerate(topDonors):
        sub = edf[edf["donor"] == sugar]
        posCounts = sub["pos"].value_counts()
        values = [int(posCounts.get(p, 0)) for p in positions[:-1]]
        otherCount = int(sub.shape[0] - sum(values))
        values.append(otherCount)
        values.append(values[0])  # 闭合 (Close)
        fig.add_trace(go.Scatterpolar(
            r=values, theta=positions + [positions[0]],
            fill="toself", name=f"{sugar} (N={sub.shape[0]:,})",
            opacity=0.5, line=dict(color=QUAL_COLORS[i % len(QUAL_COLORS)])))

    fig.update_layout(
        title="Linkage Position Preference by Sugar Type<br>"
              "<sub>Which positions does each sugar prefer as donor?</sub>",
        polar=dict(radialaxis=dict(visible=True)),
        height=650)
    saveHtml(fig, "S8_linkage_position_radar.html", height=650, width=800)


# ═════════════════════════════════════════════════════════════════════════
# S9: 稀有糖/脱氧糖分布 (Rare & Deoxy Sugar Distribution)
# ═════════════════════════════════════════════════════════════════════════

def chartS9RareSugarDist(df: pd.DataFrame):
    """Distribution of rare/unusual sugars across saponin types.
    脱氧糖 (D-Qui, L-Ole, D-Dig, D-Cym, L-Fuc) 是区分
    强心苷 (Cardiac Glycosides) 与普通皂苷的化学标志。
    """
    print("[S9] Rare & Deoxy Sugar Distribution...")
    # 定义稀有/脱氧糖集合 (Define rare/deoxy sugar set)
    rareSugars = {
        "D-Qui", "L-Ole", "D-Dig", "D-Cym", "D-Abe", "L-Fuc",
        "D-GlcN", "D-GalN", "D-GlcNAc", "D-GalNAc",
        "L-Asc", "L-Col", "D-Api", "KDO", "Neu5Ac", "Neu5Gc",
        "D-dGlc", "L-dMan", "D-Aco",
    }
    rows = []
    for _, r in df.iterrows():
        sapType = r.get("Saponin_Type", "Other")
        detailClass = str(r.get("Detailed_NP_Class", ""))
        for sugar in set(r.get("Mono_List", [])):
            if sugar in rareSugars:
                rows.append({"Sugar": sugar, "Type": sapType, "Detail": detailClass})
    if not rows:
        print("  [SKIP] No rare sugars found"); return

    rdf = pd.DataFrame(rows)
    topR = rdf["Sugar"].value_counts().head(12).index
    rdf = rdf[rdf["Sugar"].isin(topR)]

    fig = px.histogram(rdf, x="Sugar", color="Type", barmode="group",
                       color_discrete_map={"Steroidal": "#4C78A8", "Triterpenoid": "#E45756", "Other": "#999"},
                       text_auto=True)
    fig.update_layout(
        title="Rare & Deoxy Sugar Distribution in Saponins<br>"
              "<sub>Deoxy sugars (Qui, Ole, Dig, Cym, Fuc) are biomarkers for cardiac glycosides</sub>",
        xaxis_title="Sugar", yaxis_title="Count", height=550)
    saveHtml(fig, "S9_rare_deoxy_sugar.html", height=550)


# ═════════════════════════════════════════════════════════════════════════
# S10: 糖组合共现网络 (Sugar Co-occurrence Network)
# ═════════════════════════════════════════════════════════════════════════

def chartS10SugarCooccurrence(df: pd.DataFrame):
    """Network of sugar-sugar co-occurrence within same molecule.
    揭示哪些糖经常在同一分子中一起出现。
    例如 Glc+Rha+GlcA 是三萜皂苷的黄金三角。
    """
    print("[S10] Sugar Co-occurrence Network...")
    import math
    pairCounts = Counter()
    for lst in df["Mono_List"]:
        if not isinstance(lst, list) or len(lst) < 2:
            continue
        uniq = sorted(set(lst))
        for i in range(len(uniq)):
            for j in range(i+1, len(uniq)):
                pairCounts[(uniq[i], uniq[j])] += 1

    topPairs = pairCounts.most_common(20)
    if not topPairs:
        print("  [SKIP] No co-occurrence data"); return

    # 构建节点和边 (Build nodes and edges)
    allNodes = set()
    for (a, b), _ in topPairs:
        allNodes.add(a); allNodes.add(b)
    allNodes = sorted(allNodes)
    nNodes = len(allNodes)
    posMap = {}
    for i, n in enumerate(allNodes):
        angle = 2 * math.pi * i / nNodes
        posMap[n] = (math.cos(angle), math.sin(angle))

    edgeX, edgeY, edgeW = [], [], []
    for (a, b), cnt in topPairs:
        x0, y0 = posMap[a]; x1, y1 = posMap[b]
        edgeX += [x0, x1, None]; edgeY += [y0, y1, None]
        edgeW.append(cnt)

    # 节点大小按出现次数 (Node size by total count)
    nodeCounts = Counter()
    for lst in df["Mono_List"]:
        if isinstance(lst, list):
            for s in set(lst):
                if s in allNodes:
                    nodeCounts[s] += 1

    edgeTrace = go.Scatter(x=edgeX, y=edgeY, mode="lines",
                           line=dict(width=1, color="#bbb"), hoverinfo="none")

    nodeX = [posMap[n][0] for n in allNodes]
    nodeY = [posMap[n][1] for n in allNodes]
    nodeSize = [max(12, min(50, nodeCounts.get(n, 1) / 50)) for n in allNodes]
    nodeText = [f"{n} (N={nodeCounts.get(n,0):,})" for n in allNodes]

    nodeTrace = go.Scatter(
        x=nodeX, y=nodeY, mode="markers+text",
        marker=dict(size=nodeSize, color=QUAL_COLORS[:nNodes],
                    line=dict(width=1, color="white")),
        text=allNodes, textposition="top center", textfont=dict(size=10),
        hovertext=nodeText, hoverinfo="text")

    fig = go.Figure(data=[edgeTrace, nodeTrace])
    fig.update_layout(
        title="Sugar Co-occurrence Network<br>"
              "<sub>Edge = sugars appearing in same molecule; node size ∝ frequency</sub>",
        showlegend=False, height=700,
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
    saveHtml(fig, "S10_sugar_cooccurrence_network.html", height=700)


# ═════════════════════════════════════════════════════════════════════════
# S11: 多糖片段匹配 (Polysaccharide Motif Matching)
# ═════════════════════════════════════════════════════════════════════════

def chartS11PolysaccharideMotifMatch(df: pd.DataFrame):
    """Match known polysaccharide structural motifs in saponin sugar chains.
    检查已知多糖的核心二糖重复单元是否出现在皂苷糖链中。
    这揭示了皂苷与多糖生物合成途径的潜在进化关联。
    """
    print("[S11] Polysaccharide Motif Matching...")
    # 已知多糖核心片段 (Known polysaccharide core motifs)
    polysaccharideMotifs = {
        "Cellulose":     ("D-Glc", "b1-4", "D-Glc"),     # β(1→4)-glucan
        "Starch/Amylose":("D-Glc", "a1-4", "D-Glc"),     # α(1→4)-glucan
        "Chitin":        ("D-GlcNAc", "b1-4", "D-GlcNAc"),# β(1→4)-GlcNAc
        "Laminarin":     ("D-Glc", "b1-3", "D-Glc"),     # β(1→3)-glucan
        "Galactan":      ("D-Gal", "b1-4", "D-Gal"),     # β(1→4)-galactan
        "Mannan":        ("D-Man", "a1-6", "D-Man"),      # α(1→6)-mannan
        "Rhamnan":       ("L-Rha", "a1-3", "L-Rha"),     # α(1→3)-rhamnan
        "Xylan":         ("D-Xyl", "b1-4", "D-Xyl"),     # β(1→4)-xylan
        "Pectin_HG":     ("D-GalA", "a1-4", "D-GalA"),   # Homogalacturonan
        "Pectin_RG-I":   ("L-Rha", "a1-4", "D-GalA"),    # Rhamnogalacturonan I
        "Hyaluronan":    ("D-GlcA", "b1-3", "D-GlcNAc"), # Hyaluronic acid
        "Heparin":       ("L-IdoA", "a1-4", "D-GlcN"),   # Heparin/HS
    }

    # 在所有序列中搜索匹配 (Search for matches in all sequences)
    matchCounts = {k: 0 for k in polysaccharideMotifs}
    for seq in df["Sugar_Sequence"].dropna():
        seqStr = str(seq)
        for motifName, (s1, link, s2) in polysaccharideMotifs.items():
            pattern = f"{re.escape(s1)}-({re.escape(link)})-{re.escape(s2)}"
            # 允许括号前缀 (Allow parenthesized linkage format)
            patternAlt = f"{re.escape(s1)}-\\({re.escape(link)}\\)-{re.escape(s2)}"
            if re.search(patternAlt, seqStr):
                matchCounts[motifName] += 1

    # 过滤有匹配的 (Filter to those with matches)
    matched = {k: v for k, v in matchCounts.items() if v > 0}
    if not matched:
        print("  [SKIP] No polysaccharide motifs matched"); return

    # 排序并绘图 (Sort and plot)
    sortedMotifs = sorted(matched.items(), key=lambda x: x[1], reverse=True)
    names = [f"{n}\n({polysaccharideMotifs[n][0]}-{polysaccharideMotifs[n][1]}-{polysaccharideMotifs[n][2]})"
             for n, _ in sortedMotifs]
    counts = [c for _, c in sortedMotifs]

    fig = go.Figure(go.Bar(
        x=names, y=counts, marker_color="#54A24B",
        text=counts, textposition="outside"))
    fig.update_layout(
        title="Polysaccharide Motif Matches in Saponin Glycans<br>"
              "<sub>Core disaccharide repeat units of known polysaccharides found in saponin sugar chains</sub>",
        xaxis_title="Polysaccharide Motif", yaxis_title="Saponin Molecules Containing Motif",
        height=550)
    saveHtml(fig, "S11_polysaccharide_motif_match.html", height=550)


# ═════════════════════════════════════════════════════════════════════════
# S12: 糖多样性指数 (Sugar Diversity Index)
# ═════════════════════════════════════════════════════════════════════════

def chartS12SugarDiversityIndex(df: pd.DataFrame):
    """Shannon diversity index of sugar composition per molecule.
    糖多样性高 = 该分子含有多种不同的糖 (如 Glc+Rha+GlcA+Xyl)。
    糖多样性低 = 该分子只含单一糖 (如纯 Glc 链)。
    """
    print("[S12] Sugar Diversity Index Distribution...")
    import math as _math

    diversities = []
    for _, r in df.iterrows():
        lst = r["Mono_List"]
        if not isinstance(lst, list) or len(lst) == 0:
            continue
        total = len(lst)
        counts = Counter(lst)
        # Shannon diversity index H = -Σ(pi * ln(pi))
        H = 0
        for cnt in counts.values():
            pi = cnt / total
            if pi > 0:
                H -= pi * _math.log(pi)
        diversities.append({
            "H": round(H, 3),
            "Type": r.get("Saponin_Type", "Other"),
            "N_sugars": total,
            "N_unique": len(counts),
        })

    if not diversities: print("  [SKIP]"); return
    hdf = pd.DataFrame(diversities)
    hdf = hdf[hdf["Type"].isin(["Steroidal", "Triterpenoid"])]

    fig = px.histogram(hdf, x="H", color="Type", nbins=30, barmode="overlay",
                       opacity=0.7,
                       color_discrete_map={"Steroidal": "#4C78A8", "Triterpenoid": "#E45756"},
                       labels={"H": "Shannon Diversity Index (H)"})
    fig.update_layout(
        title="Sugar Diversity Index (Shannon H) — Steroidal vs Triterpenoid<br>"
              "<sub>H=0: single sugar type; H>1: highly diverse sugar composition</sub>",
        yaxis_title="Molecule Count", height=500)
    saveHtml(fig, "S12_sugar_diversity_index.html", height=500)


# ═════════════════════════════════════════════════════════════════════════
# S13: 糖链长度 × 皂苷细类 (Chain Length × Detailed NP Class)
# ═════════════════════════════════════════════════════════════════════════

def chartS13ChainLenVsDetailedClass(df: pd.DataFrame):
    """Heatmap: sugar chain length vs top detailed NP classes.
    糖链长度是否与苷元骨架类型相关？例如人参皂苷通常有更长的糖链。
    """
    print("[S13] Chain Length × Detailed NP Class...")
    sub = df[df["Detailed_NP_Class"].notna() & (df["Detailed_NP_Class"] != "")].copy()
    topC = sub["Detailed_NP_Class"].value_counts().head(10).index
    sub = sub[sub["Detailed_NP_Class"].isin(topC)]
    sub["Sugar_Bin"] = sub["Total_Sugar_Count"].clip(upper=7).astype(int).astype(str)

    piv = sub.groupby(["Detailed_NP_Class", "Sugar_Bin"]).size().unstack(fill_value=0)
    piv = piv.reindex(columns=sorted(piv.columns, key=int), fill_value=0)

    # 截断长标签 (Truncate long labels)
    piv.index = [s if len(s) <= 30 else s[:27] + "..." for s in piv.index]
    rowTot = piv.sum(axis=1)
    yLab = [f"{c} (N={rowTot[c]:,})" for c in piv.index]

    fig = px.imshow(piv.values, x=piv.columns.tolist(), y=yLab,
                    color_continuous_scale=HEATMAP_COLORSCALE, aspect="auto",
                    labels=dict(color="Count"), text_auto="d")
    fig.update_layout(
        title="Detailed NP Class × Sugar Count Distribution<br>"
              "<sub>Are certain aglycone types associated with longer sugar chains?</sub>",
        xaxis_title="Total Sugar Count", yaxis_title="Detailed NP Class",
        height=600)
    saveHtml(fig, "S13_chain_vs_detailed_class.html", height=600)


# ═════════════════════════════════════════════════════════════════════════
# S14: 各 NP 细类 糖指纹雷达 (Sugar Fingerprint Radar per Detail Class)
# ═════════════════════════════════════════════════════════════════════════

def chartS14SugarFingerprintRadar(df: pd.DataFrame):
    """Sugar composition fingerprint per detailed NP class (top-6 classes).
    将每个细类的糖组成绘制为雷达图——揭示每种骨架类型的糖偏好特征。
    """
    print("[S14] Sugar Fingerprint Radar by Detailed NP Class...")
    sub = df[df["Detailed_NP_Class"].notna() & (df["Detailed_NP_Class"] != "")]
    topC = sub["Detailed_NP_Class"].value_counts().head(6).index.tolist()

    # 全局 Top-8 糖作为雷达轴 (Global top-8 sugars as axes)
    globalSugars = Counter()
    classData = {}
    for cls in topC:
        subset = sub[sub["Detailed_NP_Class"] == cls]
        sugars = []
        for lst in subset["Mono_List"]:
            if isinstance(lst, list):
                sugars.extend(lst)
        sc = Counter(sugars)
        classData[cls] = sc
        globalSugars.update(sc)

    topSugars = [s for s, _ in globalSugars.most_common(8)]

    fig = go.Figure()
    for i, cls in enumerate(topC):
        sc = classData[cls]
        total = sum(sc.values()) or 1
        values = [sc.get(s, 0) / total * 100 for s in topSugars]
        values.append(values[0])  # 闭合
        shortLabel = cls if len(cls) <= 25 else cls[:22] + "..."
        fig.add_trace(go.Scatterpolar(
            r=values, theta=topSugars + [topSugars[0]],
            fill="toself", name=shortLabel, opacity=0.5,
            line=dict(color=QUAL_COLORS[i % len(QUAL_COLORS)])))

    fig.update_layout(
        title="Sugar Composition Fingerprint — Top 6 Detailed NP Classes<br>"
              "<sub>Each class shows distinct sugar preferences (% of total sugar occurrences)</sub>",
        polar=dict(radialaxis=dict(visible=True)),
        height=700)
    saveHtml(fig, "S14_sugar_fingerprint_radar.html", height=700, width=850)


# =========================================================================
# Part 3 函数库 — main() 在 generate_saponin_charts.py 主脚本统一调用
# Part 3 function library — called from main script via importlib
# =========================================================================

