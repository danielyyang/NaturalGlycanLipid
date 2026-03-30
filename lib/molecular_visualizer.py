"""
Phase 7: 可视化与导出引擎
Phase 7: Visualization & Export Engine

三色高亮体系 (Three-color Highlighting):
  浅红色 = 糖核心骨架 (Sugar core skeleton)
  浅黄色 = 修饰基团原子 (Modification group atoms: O-Ac, NAc, Sulfate, etc.)
  浅蓝色 = 苷元 (Aglycon / non-sugar)
"""
import io
import os
import sys
import tempfile
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, AllChem

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__))))

# 三色定义 (Three-color definitions) — RGB, 0-1 范围
COLOR_SUGAR_CORE = (1.0, 0.70, 0.70)   # 浅红色 — 糖核心骨架 (Light red)
COLOR_MODIFICATION = (1.0, 0.90, 0.50) # 浅黄色 — 修饰基团 (Light yellow)
COLOR_AGLYCON = (0.70, 0.80, 1.0)      # 浅蓝色 — 苷元 (Light blue)

DEFAULT_IMG_SIZE = (450, 300)


# =====================================================================
# 1. 糖/苷元原子高亮 (Glycan/Aglycon Atom Highlighting)
# =====================================================================

def identifySugarAtomZones(
    mol: Chem.Mol,
    sugarUnits: List[Dict],
) -> Tuple[Set[int], Set[int]]:
    """
    基于 Phase 2 糖苷键边界的精确原子归属。
    Precise atom zone assignment based on Phase 2 glycosidic bond boundaries.

    算法 (Algorithm):
      1. 用 findGlycosidicBonds() 识别 sugar_to_aglycon 键
      2. 将这些键的两个端点之间的边标记为"围栏" (fence)
      3. 从所有糖环原子出发做 BFS 洪水填充, 不越过围栏
      4. BFS 可达的所有原子 = 糖区 (包括 C6-OH, NAc, 桥 O, 糖间连接)
      5. 不可达的原子 = 苷元区

    这修正了旧代码的 Bug:
      ✗ 旧: "非环即蓝" → C6-OH 被误标蓝
      ✓ 新: BFS 洪水填充 → C6-OH, NAc 等全部归入糖区

    Args:
        mol: RDKit Mol 对象
        sugarUnits: find_mapped_sugar_units() 返回的糖单元列表

    Returns:
        (glycan_atom_set, aglycon_atom_set)
    """
    from bond_cleavage_engine import findGlycosidicBonds

    allAtoms = set(range(mol.GetNumAtoms()))

    if not sugarUnits:
        return set(), allAtoms

    # =====================================================================
    # v2.1: 分支遍历法 (Branch Traversal Method)
    # 从糖环向外遍历每个分支, 按大小分类:
    #   ≤ MAX_SUBSTITUENT_SIZE → 修饰/取代基 (黄色, 归入 glycan)
    #   > MAX_SUBSTITUENT_SIZE → 苷元 (蓝色)
    # v2.1: traverse branches from ring; classify by size:
    #   ≤ threshold → substituent (yellow, part of glycan)
    #   > threshold → aglycon (blue)
    # =====================================================================
    MAX_SUBSTITUENT_SIZE = 10  # 与 glycan_topology.py 统一 (unified)

    allSugarRingAtoms: Set[int] = set()
    for unit in sugarUnits:
        allSugarRingAtoms.update(unit["ring_atoms"])

    sugarRingAtoms = set(allSugarRingAtoms)
    substituentAtoms: Set[int] = set()
    aglyconAtoms: Set[int] = set()
    glycanAtoms: Set[int] = set(allSugarRingAtoms)

    def _getBranch(startNode: int, originNode: int):
        """BFS 遍历分支, 不进入其他糖环。返回 (分支原子集, 是否连接到其他糖环)"""
        branchNodes = set([startNode])
        queue = [startNode]
        isLinkage = False
        while queue:
            curr = queue.pop(0)
            for nbr in mol.GetAtomWithIdx(curr).GetNeighbors():
                nbrIdx = nbr.GetIdx()
                if nbrIdx == originNode or nbrIdx in branchNodes:
                    continue
                if nbrIdx in allSugarRingAtoms:
                    isLinkage = True
                    continue
                branchNodes.add(nbrIdx)
                queue.append(nbrIdx)
        return branchNodes, isLinkage

    # 遍历每个糖环原子的外环邻居 (Enumerate exocyclic neighbors)
    for rIdx in sorted(allSugarRingAtoms):
        for nbr in mol.GetAtomWithIdx(rIdx).GetNeighbors():
            nbrIdx = nbr.GetIdx()
            if nbrIdx in allSugarRingAtoms or nbrIdx in substituentAtoms or nbrIdx in aglyconAtoms:
                continue

            branch, isLinkage = _getBranch(nbrIdx, rIdx)
            branchSize = len(branch)

            if isLinkage:
                if branchSize > MAX_SUBSTITUENT_SIZE:
                    aglyconAtoms.update(branch)
                else:
                    substituentAtoms.update(branch)
                    glycanAtoms.update(branch)
            else:
                if branchSize <= MAX_SUBSTITUENT_SIZE:
                    substituentAtoms.update(branch)
                    glycanAtoms.update(branch)
                else:
                    aglyconAtoms.update(branch)

    # 补全未遍历原子 (Orphaned atoms → aglycon)
    for idx in allAtoms:
        if idx not in glycanAtoms and idx not in aglyconAtoms:
            aglyconAtoms.add(idx)

    return sugarRingAtoms, substituentAtoms, aglyconAtoms


def cleaveAndCap(
    mol: Chem.Mol,
    glycanAtoms: Set[int],
    aglyconAtoms: Set[int],
    minAglyconHeavyAtoms: int = 3,
    sugarUnits: Optional[List[Dict]] = None,
) -> Tuple[Optional[str], Optional[str]]:
    """
    v3.0: 糖环锚定 + 迭代剥离方案 — 彻底分离 Glycan/Aglycon
    v3.0: Sugar-Ring-Anchor + Iterative Peeling — Clean Glycan/Aglycon Separation.

    核心设计原则 (Core Design Principles):
    1. 糖环原子 (ring_atoms) 是碎片分类的**唯一判据** — 无启发式, 无多数票
       Sugar ring atoms are the SOLE classification criterion — no heuristics
    2. 找到糖区/苷元区边界的 **所有键**, 用 FragmentOnBonds + dummy 切断
       Find ALL boundary bonds, cut with FragmentOnBonds + dummy capping
    3. GetMolFrags + fragsMolAtomMapping 追踪原子回原始分子
       Track atoms back to original molecule via fragsMolAtomMapping
    4. 含糖环原子的碎片 = glycan, 不含 = aglycon
       Fragment with sugar ring atoms = glycan, without = aglycon
    5. 迭代剥离: glycan 碎片中过大的非糖分支被二次剥离归入 aglycon
       Iterative peeling: oversized non-sugar branches in glycan are re-peeled

    Args:
        mol: 完整分子 (RDKit Mol)
        glycanAtoms: 糖区原子索引集合 (来自 identifySugarAtomZones)
        aglyconAtoms: 苷元区原子索引集合 (来自 identifySugarAtomZones)
        minAglyconHeavyAtoms: 苷元碎片的最小重原子数 (默认 3)
        sugarUnits: 糖单元列表 (如果 None 则自动检测)

    Returns:
        (aglyconSmiles, glycanSmiles): 化学上正确的 SMILES, 或 None
    """
    import re

    if not glycanAtoms or not aglyconAtoms:
        if glycanAtoms:
            return None, Chem.MolToSmiles(mol)
        else:
            return Chem.MolToSmiles(mol), None

    # 自动检测糖单元 (Auto-detect if not provided)
    if sugarUnits is None:
        try:
            from glycan_topology import find_mapped_sugar_units
            sugarUnits = find_mapped_sugar_units(mol)
        except Exception:
            sugarUnits = []
    if not sugarUnits:
        return Chem.MolToSmiles(mol), None

    # ===================================================================
    # Step 1: 收集糖环原子集 — 唯一真相源
    # Collect sugar ring atom set — the SINGLE source of truth
    # ===================================================================
    sugarRingAtomSet: Set[int] = set()
    for unit in sugarUnits:
        sugarRingAtomSet.update(unit["ring_atoms"])

    # ===================================================================
    # [已移除] 精确糖苷键切割路径
    # [REMOVED] Precision glycosidic bond cleavage path
    #
    # 原因: cleaveWithConservation 只切 C1-O 糖苷键, 而 identifySugarAtomZones
    # 按分支大小分区 — 两套边界不一致, 导致 glycan 碎片包含苷元残留。
    # 染色正确 = 分区正确, 因此必须始终按分区边界切割。
    # Reason: cleaveWithConservation cuts at C1-O glycosidic bonds, but
    # identifySugarAtomZones classifies by branch size — different boundaries.
    # Coloring is correct = zones are correct, so always cut at zone boundaries.
    # ===================================================================

    # ===================================================================
    # Step 3: 通用分离路径 — FragmentOnBonds + 糖环锚定分类
    # Universal separation path — FragmentOnBonds + sugar ring anchor
    # (用于 C-glycosides, N-glycosides, 以及 cleaveWithConservation 失败的场景)
    # ===================================================================

    # 3a: 找到 glycan/aglycon 边界的所有键
    # Find ALL bonds crossing the glycan/aglycon boundary
    boundaryBondIdxs = []
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if (a1 in glycanAtoms and a2 in aglyconAtoms) or \
           (a1 in aglyconAtoms and a2 in glycanAtoms):
            boundaryBondIdxs.append(bond.GetIdx())

    if not boundaryBondIdxs:
        # 无边界键 → 整个分子是 glycan 或 aglycon
        if sugarRingAtomSet:
            return None, Chem.MolToSmiles(mol)
        else:
            return Chem.MolToSmiles(mol), None

    # 3b: FragmentOnBonds — 切断 + dummy 封端
    fragMol = Chem.FragmentOnBonds(mol, boundaryBondIdxs, addDummies=True)
    origAtomCount = mol.GetNumAtoms()

    # 3c: GetMolFrags + 原子索引追踪
    fragsMolAtomMapping: list = []
    fragMols = Chem.GetMolFrags(
        fragMol, asMols=True, sanitizeFrags=False,
        fragsMolAtomMapping=fragsMolAtomMapping,
    )

    # ===================================================================
    # Step 4: 碎片分类 — 糖环原子是唯一判据
    # Fragment classification — sugar ring atom is the SOLE criterion
    # ===================================================================
    glycanFragSmiles = []
    aglyconFragSmiles = []

    for fragM, atomMapping in zip(fragMols, fragsMolAtomMapping):
        # atomMapping[i] = 该碎片第 i 个原子在 fragMol 中的索引
        # fragMol 中 < origAtomCount 的 = 原始原子, >= origAtomCount = dummy
        hasSugarRingAtom = False
        realAtomCount = 0
        for fragMolIdx in atomMapping:
            if fragMolIdx < origAtomCount:
                realAtomCount += 1
                if fragMolIdx in sugarRingAtomSet:
                    hasSugarRingAtom = True

        fragSmi = Chem.MolToSmiles(fragM)
        if not fragSmi:
            continue

        if hasSugarRingAtom:
            glycanFragSmiles.append(fragSmi)
        else:
            aglyconFragSmiles.append(fragSmi)

    # ===================================================================
    # Step 5: 清理 dummy 标签 + 苷元垃圾回收
    # Clean dummy labels + aglycon garbage collection
    # ===================================================================
    def _cleanDummies(smi: str) -> Optional[str]:
        if not smi:
            return None
        # 去掉 isotope 标签, 保留 [*]
        cleaned = re.sub(r'\[\d+\*\]', '[*]', smi)
        return cleaned

    glycanRaw = '.'.join(glycanFragSmiles) if glycanFragSmiles else None
    aglyconRaw = '.'.join(aglyconFragSmiles) if aglyconFragSmiles else None

    glycanFinal = _cleanDummies(glycanRaw)
    aglyconFinal = _cleanDummies(aglyconRaw)

    # 苷元垃圾回收: 只保留大碎片 (Aglycon GC: filter small debris)
    if aglyconFinal:
        subFrags = aglyconFinal.split('.')
        bigFrags = [sf for sf in subFrags
                    if (m := Chem.MolFromSmiles(sf, sanitize=False))
                    and m.GetNumHeavyAtoms() >= minAglyconHeavyAtoms]
        aglyconFinal = '.'.join(bigFrags) if bigFrags else None

    return aglyconFinal, glycanFinal


def drawHighlightedMolecule(
    smiles: str,
    sugarUnits: Optional[List[Dict]] = None,
    glycanSmiles: Optional[str] = None,
    imgSize: Tuple[int, int] = DEFAULT_IMG_SIZE,
) -> Optional[bytes]:
    """
    绘制三色高亮 2D 分子图 (v2.1 — 基于分支分类)。
    Draw 2D molecule with three-color highlighting (v2.1 — branch-based).

    三色方案 (Three-color scheme):
      浅红 = 糖环骨架原子 (Sugar ring atoms only)
      浅黄 = 糖环上修饰/取代基 (ALL substituents on sugar, ≤10 atoms)
      浅蓝 = 苷元 (Aglycon, >10 atoms or non-sugar)

    v2.1: 不再依赖 SMARTS 匹配来识别修饰; 改用分支大小自动分类。
    v2.1: no more SMARTS for modification detection; branch-size auto-classification.
    """
    if not smiles or smiles in ("NULL", "nan", "", "*"):
        return None

    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None

    # 自动检测糖单元 (Auto-detect if not provided)
    if sugarUnits is None:
        try:
            from glycan_topology import find_mapped_sugar_units
            sugarUnits = find_mapped_sugar_units(mol)
        except Exception:
            sugarUnits = []

    AllChem.Compute2DCoords(mol)

    if not sugarUnits:
        drawer = Draw.rdMolDraw2D.MolDraw2DCairo(imgSize[0], imgSize[1])
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return drawer.GetDrawingText()

    # Step 1: 三区域划分 (Three-zone partitioning)
    sugarRingAtoms, substituentAtoms, aglyconAtoms = identifySugarAtomZones(mol, sugarUnits)

    # Step 2: 直接三色分配 — 基于区域, 不依赖 SMARTS
    # Direct three-color assignment — zone-based, no SMARTS needed
    highlightAtoms = list(range(mol.GetNumAtoms()))
    highlightAtomColors = {}
    for idx in sugarRingAtoms:
        highlightAtomColors[idx] = COLOR_SUGAR_CORE       # 红 (Red)
    for idx in substituentAtoms:
        highlightAtomColors[idx] = COLOR_MODIFICATION      # 黄 (Yellow)
    for idx in aglyconAtoms:
        highlightAtomColors[idx] = COLOR_AGLYCON           # 蓝 (Blue)

    # 键高亮: 按两端原子颜色决定 (Bond color = follow atom colors)
    highlightBonds = []
    highlightBondColors = {}
    for bond in mol.GetBonds():
        bIdx = bond.GetIdx()
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        highlightBonds.append(bIdx)

        c1 = highlightAtomColors.get(a1, COLOR_AGLYCON)
        c2 = highlightAtomColors.get(a2, COLOR_AGLYCON)

        if c1 == c2:
            highlightBondColors[bIdx] = c1
        elif COLOR_MODIFICATION in (c1, c2):
            highlightBondColors[bIdx] = COLOR_MODIFICATION
        elif COLOR_SUGAR_CORE in (c1, c2) and COLOR_AGLYCON in (c1, c2):
            highlightBondColors[bIdx] = COLOR_AGLYCON
        else:
            highlightBondColors[bIdx] = COLOR_SUGAR_CORE

    drawer = Draw.rdMolDraw2D.MolDraw2DCairo(imgSize[0], imgSize[1])
    drawer.drawOptions().addAtomIndices = False
    drawer.DrawMolecule(
        mol,
        highlightAtoms=highlightAtoms,
        highlightAtomColors=highlightAtomColors,
        highlightBonds=highlightBonds,
        highlightBondColors=highlightBondColors,
    )
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


# =====================================================================
# 2. Excel 报告生成 (Excel Report Generation)
# =====================================================================

def generateExcelReport(
    df: pd.DataFrame,
    outputPath: str,
    smilesCol: str = "canonical_smiles",
    inchikeyCol: str = "standard_inchi_key",
    sugarSeqCol: str = "sugar_sequence",
    scaffoldCol: str = "murcko_scaffold",
    classCol: str = "classification",
    glycolipidCol: str = "glycolipid_flag",
    maxRows: int = 100,
    imgSize: Tuple[int, int] = DEFAULT_IMG_SIZE,
) -> str:
    """
    生成带嵌入分子图的 Excel 报告。
    Generate Excel report with embedded molecule images.

    Args:
        df: 包含所有 Phase 产出列的 DataFrame
        outputPath: 输出 .xlsx 路径
        maxRows: 最大行数
        imgSize: 分子图像尺寸

    Returns:
        输出文件路径
    """
    from glycan_topology import find_mapped_sugar_units

    # 限制行数
    reportDf = df.head(maxRows).copy()
    print(f"  Generating Excel report for {len(reportDf)} compounds...")

    # 创建临时目录存放图片 (Temp dir for images)
    tmpDir = tempfile.mkdtemp(prefix="glycan_report_")

    # 准备列 (Prepare columns)
    outputCols = [inchikeyCol, "Structure", sugarSeqCol, scaffoldCol, classCol, glycolipidCol]
    # 确保所需列存在
    for col in [sugarSeqCol, scaffoldCol, classCol, glycolipidCol]:
        if col not in reportDf.columns:
            reportDf[col] = ""

    # 生成图片 (Generate images)
    imgPaths = []
    for i, (_, row) in enumerate(reportDf.iterrows()):
        smiles = str(row.get(smilesCol, ""))
        imgPath = os.path.join(tmpDir, f"mol_{i:04d}.png")

        try:
            mol = Chem.MolFromSmiles(smiles)
            units = find_mapped_sugar_units(mol) if mol else []
            pngData = drawHighlightedMolecule(smiles, units, imgSize)

            if pngData:
                with open(imgPath, "wb") as f:
                    f.write(pngData)
                imgPaths.append(imgPath)
            else:
                imgPaths.append(None)
        except Exception:
            imgPaths.append(None)

        if (i + 1) % 20 == 0:
            print(f"    Image {i+1}/{len(reportDf)} generated")

    # 写入 Excel (Write Excel with xlsxwriter)
    # 每行需要足够高以显示图片
    rowHeight = imgSize[1] * 0.75  # px to ~points
    colWidth = imgSize[0] / 7  # px to ~char width

    with pd.ExcelWriter(outputPath, engine="xlsxwriter") as writer:
        # 写入数据列 (不含图片)
        dataCols = [c for c in outputCols if c != "Structure"]
        exportDf = reportDf[[c for c in dataCols if c in reportDf.columns]].reset_index(drop=True)
        exportDf.to_excel(writer, sheet_name="Report", startcol=0, startrow=0, index=False)

        workbook = writer.book
        worksheet = writer.sheets["Report"]

        # 设置格式 (Set formats)
        headerFmt = workbook.add_format({
            "bold": True, "bg_color": "#2C3E50", "font_color": "white",
            "border": 1, "text_wrap": True, "valign": "vcenter",
        })
        cellFmt = workbook.add_format({
            "text_wrap": True, "valign": "vcenter", "border": 1,
        })

        # 写入表头 (Write headers)
        structureColIdx = 1  # 图片列放在 InChIKey 后面
        # 先移除原有表头，重写包含 Structure 列的表头
        allHeaders = [inchikeyCol, "Structure"] + [c for c in dataCols if c != inchikeyCol]
        for colIdx, header in enumerate(allHeaders):
            worksheet.write(0, colIdx, header, headerFmt)

        # 重写数据列（移位给 Structure 留空间）
        for rowIdx in range(len(exportDf)):
            # InChIKey
            worksheet.write(rowIdx + 1, 0, str(exportDf.iloc[rowIdx].get(inchikeyCol, "")), cellFmt)
            # 其他数据列从 col=2 开始
            otherCols = [c for c in dataCols if c != inchikeyCol]
            for colIdx, col in enumerate(otherCols):
                val = str(exportDf.iloc[rowIdx].get(col, ""))
                worksheet.write(rowIdx + 1, colIdx + 2, val, cellFmt)

        # 插入图片 (Insert images)
        worksheet.set_column(structureColIdx, structureColIdx, colWidth)
        for rowIdx in range(len(reportDf)):
            worksheet.set_row(rowIdx + 1, rowHeight)
            imgPath = imgPaths[rowIdx]
            if imgPath and os.path.exists(imgPath):
                worksheet.insert_image(
                    rowIdx + 1, structureColIdx, imgPath,
                    {"x_scale": 0.8, "y_scale": 0.8, "x_offset": 5, "y_offset": 5},
                )

        # 设置其他列宽 (Set column widths)
        worksheet.set_column(0, 0, 20)  # InChIKey
        for colIdx in range(2, len(allHeaders)):
            worksheet.set_column(colIdx, colIdx, 25)

    # 清理临时文件 (Cleanup temp files)
    for p in imgPaths:
        if p and os.path.exists(p):
            try:
                os.remove(p)
            except Exception:
                pass
    try:
        os.rmdir(tmpDir)
    except Exception:
        pass

    print(f"  Report saved: {outputPath}")
    return outputPath


# =====================================================================
# 3. HTML 报告生成 (HTML Report Generation) — 轻量备选
# =====================================================================

def generateHtmlReport(
    df: pd.DataFrame,
    outputPath: str,
    smilesCol: str = "canonical_smiles",
    inchikeyCol: str = "standard_inchi_key",
    sugarSeqCol: str = "sugar_sequence",
    scaffoldCol: str = "murcko_scaffold",
    classCol: str = "classification",
    glycolipidCol: str = "glycolipid_flag",
    maxRows: int = 100,
    imgSize: Tuple[int, int] = (350, 250),
) -> str:
    """
    生成带内联分子图的 HTML 报告 (使用 Base64 编码图片)。
    Generate HTML report with inline molecule images (Base64 encoded).
    """
    import base64
    from glycan_topology import find_mapped_sugar_units

    reportDf = df.head(maxRows)

    htmlParts = ["""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>GlycoNP Visual Report</title>
<style>
  body { font-family: 'Segoe UI', Arial, sans-serif; background: #f5f6fa; margin: 20px; }
  h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
  .legend { display: flex; gap: 20px; margin: 15px 0; font-size: 14px; }
  .legend span { padding: 4px 12px; border-radius: 4px; }
  .glycan { background: rgba(255,180,180,0.6); }
  .aglycon { background: rgba(180,200,255,0.6); }
  .bridge { background: rgba(255,230,140,0.6); }
  table { border-collapse: collapse; width: 100%; background: white; box-shadow: 0 2px 8px rgba(0,0,0,0.1); }
  th { background: #2c3e50; color: white; padding: 10px; text-align: left; position: sticky; top: 0; }
  td { padding: 8px; border: 1px solid #ddd; vertical-align: middle; }
  tr:nth-child(even) { background: #f8f9fa; }
  tr:hover { background: #eef2f7; }
  img.mol { border: 1px solid #e0e0e0; border-radius: 4px; }
  .tag { display: inline-block; padding: 2px 8px; border-radius: 10px; font-size: 12px; margin: 2px; }
  .tag-glycolipid { background: #ffeaa7; color: #d63031; }
  .tag-class { background: #dfe6e9; color: #2d3436; }
</style>
</head><body>
<h1>GlycoNP Pipeline — Visual Report</h1>
<div class="legend">
  <span class="glycan">■ 糖核心 Sugar Core</span>
  <span class="bridge">■ 修饰基团 Modifications</span>
  <span class="aglycon">■ 苷元 Aglycon</span>
</div>
<table>
<tr><th>#</th><th>Structure</th><th>InChIKey</th><th>Glycan Sequence</th><th>Modifications</th><th>Scaffold</th><th>Classification</th></tr>
"""]

    for i, (_, row) in enumerate(reportDf.iterrows()):
        smiles = str(row.get(smilesCol, ""))
        inchikey = str(row.get(inchikeyCol, ""))[:14]

        # 生成高亮图片 (Generate highlighted image)
        imgTag = ""
        try:
            mol = Chem.MolFromSmiles(smiles)
            units = find_mapped_sugar_units(mol) if mol else []
            pngData = drawHighlightedMolecule(smiles, units, imgSize)
            if pngData:
                b64 = base64.b64encode(pngData).decode("utf-8")
                imgTag = f'<img class="mol" src="data:image/png;base64,{b64}" width="{imgSize[0]}" height="{imgSize[1]}">'
        except Exception:
            imgTag = "<em>Error</em>"

        seq = str(row.get(sugarSeqCol, ""))
        modsStr = str(row.get("Glycan_Modifications", ""))
        scaffold = str(row.get(scaffoldCol, ""))[:35]
        cls = str(row.get(classCol, ""))
        glycoFlag = str(row.get(glycolipidCol, ""))

        clsTags = f'<span class="tag tag-class">{cls}</span>' if cls and cls != "nan" else ""
        if glycoFlag and glycoFlag not in ("", "nan"):
            clsTags += f' <span class="tag tag-glycolipid">{glycoFlag}</span>'

        modsDisplay = modsStr if modsStr and modsStr != "nan" else ""

        htmlParts.append(
            f"<tr><td>{i+1}</td><td>{imgTag}</td><td>{inchikey}</td>"
            f"<td>{seq}</td><td>{modsDisplay}</td><td><code>{scaffold}</code></td><td>{clsTags}</td></tr>\n"
        )

        if (i + 1) % 20 == 0:
            print(f"    HTML row {i+1}/{len(reportDf)}")

    htmlParts.append("</table></body></html>")

    with open(outputPath, "w", encoding="utf-8") as f:
        f.write("".join(htmlParts))

    print(f"  HTML report saved: {outputPath}")
    return outputPath
