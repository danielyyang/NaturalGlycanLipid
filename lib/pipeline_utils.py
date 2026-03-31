"""
[EN] Pipeline Core Utilities (pipeline_utils.py)
     Reusable structural functions shared across all pipeline scripts.
     Contains the data stripping rules and encapsulation logic:
       - Natural Language Processing (NLP) fallback extraction.
       - Rare sugar cross-validation against stereochemical limits.
       - Inter-sugar and Aglycone-sugar Glycosidic bond Directional computations.
       - Universal SVG/PNG 3-color structural rendering.

[CN] GlycoNP 管线核心抽象运算库 (pipeline_utils.py)
     从 V12/V13 主干提炼的通用剥离与合并规则，供所有脚本公共调用。
     核心涵盖：
       - 基于名称的 NLP 糖类救援与罕见糖基纠错。
       - 基于立体化学树的异头碳 C1 寻路与断键向量生成。
       - 识别并切除多修饰基团的尺寸控制逻辑 (SizeThreshold <= 3)。
       - 自动化三色高亮渲染输出：糖(红), 修饰(黄), 苷元(蓝)。
"""
import re
import json
import base64
from typing import Optional

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

from lib.glycan_topology import find_mapped_sugar_units


# =====================================================================
# 罕见糖-名称交叉验证 常量 (Rare Sugar Cross-Validation Constants)
# =====================================================================
RARE_SUGARS = {
    "D-Tal", "L-Tal", "D-All", "L-All", "D-Alt", "L-Alt",
    "D-Gul", "L-Gul", "D-Ido", "L-Ido", "D-TalN", "L-TalN",
}

NAME_SUGAR_MAP = [
    (re.compile(r'glucuronic|glucuronide|glucuronosyl', re.I), 'D-GlcA'),
    (re.compile(r'glucosinolat|glucopyranosid|glucosid|glucosyl|glucofuranos|gluco(?:se)?(?![a-z])', re.I), 'D-Glc'),
    (re.compile(r'galactopyranosid|galactolipid|galactosid|galactosyl|galacto(?:se)?(?![a-z])', re.I), 'D-Gal'),
    (re.compile(r'mannopyranosid|mannoprotein|mannosid|mannosyl|manno(?:se)?(?![a-z])', re.I), 'D-Man'),
    (re.compile(r'xylopyranosid|xylosid|xylosyl|xylan|xylo(?:se)?(?![a-z])', re.I), 'D-Xyl'),
    (re.compile(r'arabino(?:se|syl|furanos|pyranos)?(?![a-z])', re.I), 'L-Ara'),
    (re.compile(r'rhamno(?:se|syl|pyranos)?(?![a-z])', re.I), 'L-Rha'),
    (re.compile(r'fuco(?:se|syl|pyranos)?(?![a-z])', re.I), 'L-Fuc'),
    (re.compile(r'apiose|apiosyl|apifuranos', re.I), 'D-Api'),
    (re.compile(r'quinovose|quinovosyl', re.I), 'D-Qui'),
    (re.compile(r'digitalose', re.I), 'D-Dig'),
    (re.compile(r'digitoxose', re.I), 'D-Dtx'),
    (re.compile(r'oleandrose', re.I), 'L-Ole'),
    (re.compile(r'cymarose', re.I), 'D-Cym'),
    (re.compile(r'thevetose', re.I), 'D-Thv'),
    (re.compile(r'boivinose', re.I), 'D-Boi'),
]

GENERIC_PAT = re.compile(r'\b(Hex|Pen|dHex|HexA|Non|Oct|Hept)\b')


# =====================================================================
# NLP: 名称→糖 提取 (Name-based Sugar Extraction)
# =====================================================================
def extractSugarFromName(
    name: Optional[str],
    iupacName: Optional[str],
    synonyms: Optional[str],
    abstract: Optional[str] = None,
) -> Optional[str]:
    """从名称字段中提取糖线索。
    Extract sugar clue from name fields via regex matching.

    设计意图: 当化学结构匹配返回泛指标签 (Hex/Pen) 时,
    文献名称可以提供具体糖名。这是管线中的"NLP 救援"路径。
    Design: When structural matching returns generic labels,
    literature names can provide specific sugar names.

    Args:
        name: 化合物名称 (Compound name)
        iupacName: IUPAC 名称
        synonyms: 同义词
        abstract: 文献摘要

    Returns:
        str or None: 推断的糖名 (e.g. "D-Glc") 或 None
    """
    textPool = " ".join(filter(None, [
        str(x) if x and str(x) != "nan" else None
        for x in [name, iupacName, synonyms, abstract]
    ]))
    if not textPool.strip():
        return None
    for regex, sugar in NAME_SUGAR_MAP:
        if regex.search(textPool):
            return sugar
    return None


def crossValidateRareSugars(seq: str, nameSugar: Optional[str]) -> str:
    """替换罕见糖为名称推断的常见糖。
    Replace rare sugar with name-inferred common sugar.

    设计意图: 结构匹配可能因手性模糊命中罕见糖 (D-Tal, D-All),
    而文献名称明确指向常见糖 (D-Glc) → 用名称结果替换。
    Design: Structural matching may hit rare sugars due to chirality
    ambiguity. If name clearly indicates a common sugar, replace.

    Args:
        seq: 糖序列字符串 (Sugar sequence string)
        nameSugar: 名称推断的糖 (Name-inferred sugar)

    Returns:
        str: 修正后的序列
    """
    if not nameSugar:
        return seq
    for rare in RARE_SUGARS:
        if rare in seq:
            seq = seq.replace(rare, nameSugar)
    return seq


# =====================================================================
# 异头碳查找 (Anomeric Carbon Finder)
# =====================================================================
def _findAnomericCarbon(mol, ringAtoms):
    """找到糖环的 anomeric carbon (C1): 与环内 O 相邻的环碳。
    Find anomeric carbon (C1): ring carbon adjacent to the ring oxygen.

    如果两个环碳都与环 O 相邻, 返回有环外杂原子 (O/N/S) 邻居的那个。
    If two ring carbons are adjacent to ring O, return the one with an
    exocyclic O/N/S neighbor (the glycosidic bond position).

    Args:
        mol: RDKit molecule
        ringAtoms: list of atom indices in the sugar ring

    Returns:
        int or None: atom index of the anomeric carbon
    """
    ringSet = set(ringAtoms)
    # 找到环内氧 (Find ring oxygen)
    ringOxygenIdx = None
    for idx in ringAtoms:
        if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8:
            ringOxygenIdx = idx
            break
    if ringOxygenIdx is None:
        return None

    # 找与环氧相邻的环碳 (Find ring carbons adjacent to ring oxygen)
    ringO = mol.GetAtomWithIdx(ringOxygenIdx)
    candidateCarbons = []
    for nbr in ringO.GetNeighbors():
        if nbr.GetIdx() in ringSet and nbr.GetAtomicNum() == 6:
            candidateCarbons.append(nbr.GetIdx())

    if not candidateCarbons:
        return None
    if len(candidateCarbons) == 1:
        return candidateCarbons[0]

    # 两个候选: 选择有环外杂原子 (O/N/S) 邻居的那个 → 这是 C1
    # Two candidates: pick the one with exocyclic heteroatom neighbor → C1
    for cIdx in candidateCarbons:
        cAtom = mol.GetAtomWithIdx(cIdx)
        for nbr in cAtom.GetNeighbors():
            nIdx = nbr.GetIdx()
            if nIdx in ringSet:
                continue
            if nbr.GetAtomicNum() in (7, 8, 16):  # O, N, S
                return cIdx
    # 如果都没有环外杂原子, 返回第一个 (Fallback)
    return candidateCarbons[0]


# =====================================================================
# 多糖苷键检测 (Glycosidic Bond Detection)
# =====================================================================
def detectAllGlycosidicBonds(mol, sugarUnits):
    """V3: 仅从 anomeric carbon (C1) 检测 Sugar→Aglycone 糖苷键。
    V3: Only detect glycosidic bonds from anomeric carbon (C1).

    对于每个糖单元, 找到 C1 (与环 O 相邻的碳), 然后检查 C1 的
    环外邻居是否连接 aglycone. 每条独立链产出 1 个 Bond Detail.
    For each sugar unit, find C1, then check if C1's exocyclic neighbor
    connects to aglycone. Each independent chain produces 1 Bond Detail.

    Returns:
        rootBondType: str — 根部糖苷键类型 (e.g. "β-O-linked")
        bondDetailJson: str — JSON 列表, 每条 reducing-end 糖苷键
    """
    if not mol or not sugarUnits:
        return "", "[]"

    try:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    except Exception:
        pass

    allSugarAtoms = set()
    unitIdMap = {}  # atomIdx → unit index
    for uid, u in enumerate(sugarUnits):
        for idx in u.get("ring_atoms", []):
            allSugarAtoms.add(idx)
            unitIdMap[idx] = uid
        for idx in u.get("position_map", {}).keys():
            allSugarAtoms.add(idx)
            if idx not in unitIdMap:
                unitIdMap[idx] = uid
        for idx in u.get("oxygen_map", {}).keys():
            allSugarAtoms.add(idx)
            if idx not in unitIdMap:
                unitIdMap[idx] = uid
        for rIdx in u.get("ring_atoms", []):
            atom = mol.GetAtomWithIdx(rIdx)
            for nbr in atom.GetNeighbors():
                nIdx = nbr.GetIdx()
                if nIdx not in unitIdMap:
                    unitIdMap[nIdx] = uid
                    allSugarAtoms.add(nIdx)

    # 构建糖→糖连接图 (Build sugar-sugar graph)
    sugarGraph = {uid: set() for uid in range(len(sugarUnits))}
    for uid, u in enumerate(sugarUnits):
        ringSet = set(u.get("ring_atoms", []))
        for rIdx in ringSet:
            atom = mol.GetAtomWithIdx(rIdx)
            for nbr in atom.GetNeighbors():
                nIdx = nbr.GetIdx()
                if nIdx in ringSet:
                    continue
                bridge = nbr
                for bNbr in bridge.GetNeighbors():
                    bIdx = bNbr.GetIdx()
                    if bIdx == rIdx:
                        continue
                    if bIdx in unitIdMap and unitIdMap[bIdx] != uid:
                        sugarGraph[uid].add(unitIdMap[bIdx])
                        sugarGraph[unitIdMap[bIdx]].add(uid)

    # Phase 1: Sugar→Aglycone 键检测
    bondDetails = []
    rootBondType = ""
    seenReducingEndUnits = set()

    for uid, u in enumerate(sugarUnits):
        ringAtoms = u.get("ring_atoms", [])
        ringSet = set(ringAtoms)
        sugarName = u.get("name", f"Sugar_{uid}")

        c1Idx = _findAnomericCarbon(mol, ringAtoms)
        if c1Idx is None:
            continue

        c1Atom = mol.GetAtomWithIdx(c1Idx)

        for nbr in c1Atom.GetNeighbors():
            nIdx = nbr.GetIdx()
            if nIdx in ringSet:
                continue

            atomicNum = nbr.GetAtomicNum()
            linkType = ""
            if atomicNum == 8:   linkType = "O"
            elif atomicNum == 16: linkType = "S"
            elif atomicNum == 7:  linkType = "N"
            elif atomicNum == 6:  linkType = "C"
            if not linkType:
                continue

            targetType = "Aglycone"
            for bNbr in nbr.GetNeighbors():
                bIdx = bNbr.GetIdx()
                if bIdx == c1Idx:
                    continue
                if bIdx in unitIdMap and unitIdMap[bIdx] != uid:
                    targetType = "Sugar"
                    break

            if targetType != "Aglycone":
                continue

            # 修饰基团过滤: fragSize ≤ 3 → 修饰, 不是苷元
            # Modification filter: fragSize ≤ 3 → modification, not aglycone
            # 仅用 ring_atoms 作为停止集 (不含 position_map)
            sugarRingOnly = set()
            for su in sugarUnits:
                sugarRingOnly.update(su.get("ring_atoms", []))

            fragSize = 0
            visited = {c1Idx}
            stack = [nIdx]
            while stack:
                cur = stack.pop()
                if cur in visited:
                    continue
                visited.add(cur)
                if cur in sugarRingOnly and cur != nIdx:
                    continue
                fragSize += 1
                for adj in mol.GetAtomWithIdx(cur).GetNeighbors():
                    if adj.GetIdx() not in visited:
                        stack.append(adj.GetIdx())
            # 统一阈值: 所有键类型均为 3
            sizeThreshold = 3
            if fragSize <= sizeThreshold:
                continue

            cipCode = ""
            if c1Atom.HasProp("_CIPCode"):
                cipCode = c1Atom.GetProp("_CIPCode")
            anomer = "?"
            if cipCode == "S": anomer = "α"
            elif cipCode == "R": anomer = "β"

            bondTypeStr = f"{anomer}-{linkType}-linked"

            if uid in seenReducingEndUnits:
                continue
            seenReducingEndUnits.add(uid)

            detail = {
                "sugar": sugarName,
                "sugar_id": uid,
                "target": "Aglycone",
                "target_type": "Aglycone",
                "bond": bondTypeStr,
            }
            bondDetails.append(detail)

            if not rootBondType:
                rootBondType = bondTypeStr
            break

    # Phase 2: Inter-Sugar 键 (仅纯糖/寡糖)
    if not bondDetails:
        for uid, u in enumerate(sugarUnits):
            ringAtoms = u.get("ring_atoms", [])
            ringSet = set(ringAtoms)
            sugarName = u.get("name", f"Sugar_{uid}")

            c1Idx = _findAnomericCarbon(mol, ringAtoms)
            if c1Idx is None:
                continue

            c1Atom = mol.GetAtomWithIdx(c1Idx)

            for nbr in c1Atom.GetNeighbors():
                nIdx = nbr.GetIdx()
                if nIdx in ringSet:
                    continue

                atomicNum = nbr.GetAtomicNum()
                linkType = ""
                if atomicNum == 8:   linkType = "O"
                elif atomicNum == 16: linkType = "S"
                elif atomicNum == 7:  linkType = "N"
                elif atomicNum == 6:  linkType = "C"
                if not linkType:
                    continue

                targetUid = None
                for bNbr in nbr.GetNeighbors():
                    bIdx = bNbr.GetIdx()
                    if bIdx == c1Idx:
                        continue
                    if bIdx in unitIdMap and unitIdMap[bIdx] != uid:
                        targetUid = unitIdMap[bIdx]
                        break

                if targetUid is None:
                    continue

                if uid >= targetUid:
                    continue

                cipCode = ""
                if c1Atom.HasProp("_CIPCode"):
                    cipCode = c1Atom.GetProp("_CIPCode")
                anomer = "?"
                if cipCode == "S": anomer = "α"
                elif cipCode == "R": anomer = "β"

                targetName = sugarUnits[targetUid].get("name", f"Sugar_{targetUid}")
                bondTypeStr = f"{anomer}-{linkType}-linked"

                detail = {
                    "sugar": sugarName,
                    "sugar_id": uid,
                    "target": targetName,
                    "target_type": "Sugar",
                    "bond": bondTypeStr,
                }
                bondDetails.append(detail)
                break

    return rootBondType, json.dumps(bondDetails, ensure_ascii=False)


# =====================================================================
# 物种名 → 生物分类 (Genus-Level Classification)
# =====================================================================
GENUS_TO_TYPE = {}
_PLANT_GENERA = [
    "Acacia","Aconitum","Actaea","Adonis","Aegle","Aesculus","Agave","Ageratum",
    "Allium","Aloe","Alpinia","Ammi","Andrographis","Angelica","Annona","Apium",
    "Arabidopsis","Aralia","Arctium","Aristolochia","Artemisia","Asparagus","Astragalus",
    "Atractylodes","Atropa","Azadirachta","Berberis","Beta","Betula","Brassica","Bryonia",
    "Bupleurum","Calendula","Camellia","Campsis","Cannabis","Capsella","Capsicum","Cassia",
    "Catharanthus","Cephalotaxus","Chrysanthemum","Cinchona","Cinnamomum","Cistanche","Citrus",
    "Clematis","Cnidium","Coffea","Colchicum","Commiphora","Convallaria","Coptis","Coriandrum",
    "Cornus","Crataegus","Croton","Curcuma","Cynara","Daphne","Datura","Dendrobium","Digitalis",
    "Dioscorea","Diospyros","Dracocephalum","Echinacea","Elsholtzia","Ephedra","Epimedium",
    "Erythrina","Eucalyptus","Eucommia","Eupatorium","Euphorbia","Ferula","Ficus","Forsythia",
    "Fraxinus","Galanthus","Garcinia","Gardenia","Gentiana","Ginkgo","Glycine","Glycyrrhiza",
    "Gossypium","Gynostemma","Hedera","Helianthus","Helleborus","Hordeum","Humulus","Hydrangea",
    "Hypericum","Ilex","Illicium","Indigofera","Ipomoea","Iris","Isatis","Juglans","Juncus",
    "Juniperus","Kaempferia","Lactuca","Lamium","Lavandula","Leonurus","Ligusticum","Ligustrum",
    "Lilium","Lindera","Linum","Liquidambar","Lithospermum","Lonicera","Lotus","Lycium","Lycoris",
    "Magnolia","Malus","Mangifera","Matricaria","Medicago","Melissa","Mentha","Mirabilis","Momordica",
    "Morinda","Morus","Murraya","Myristica","Nelumbo","Nicotiana","Nigella","Ocimum","Olea","Oryza",
    "Paeonia","Panax","Papaver","Petroselinum","Peucedanum","Phellodendron","Phyllanthus","Physalis",
    "Pimpinella","Pinus","Piper","Plantago","Platanus","Podophyllum","Polygala","Polygonatum",
    "Polygonum","Prunus","Psoralea","Pueraria","Punica","Quercus","Ranunculus","Rauvolfia",
    "Rehmannia","Rhamnus","Rheum","Rhododendron","Ricinus","Rosa","Rosmarinus","Rubia","Rubus",
    "Ruta","Salix","Salvia","Sambucus","Sanguinaria","Sapindus","Scrophularia","Scutellaria",
    "Senna","Sesamum","Silybum","Sinapis","Siraitia","Smilax","Solanum","Solidago","Sophora",
    "Stevia","Strychnos","Swertia","Syzygium","Tanacetum","Taraxacum","Taxus","Terminalia",
    "Theobroma","Thuja","Thymus","Tilia","Tribulus","Trichosanthes","Trigonella","Triticum",
    "Tussilago","Tylophora","Uncaria","Urtica","Valeriana","Vanilla","Veratrum","Verbena",
    "Vinca","Viola","Viscum","Vitex","Vitis","Withania","Xanthium","Zanthoxylum","Zingiber","Ziziphus",
]
_FUNGI_GENERA = [
    "Acremonium","Alternaria","Amanita","Aspergillus","Aureobasidium","Beauveria",
    "Bipolaris","Botrytis","Candida","Chaetomium","Cladosporium","Claviceps","Colletotrichum",
    "Cordyceps","Cryptococcus","Curvularia","Daldinia","Diaporthe","Emericella","Endothia",
    "Epicoccum","Eupenicillium","Fusarium","Ganoderma","Gliocladium","Hypocrea","Inonotus",
    "Lentinus","Metarhizium","Monascus","Mucor","Myceliophtho","Neosartorya","Nigrospora",
    "Paecilomyces","Penicillium","Pestalotiopsis","Phoma","Phomopsis","Pleurotus","Pycnoporus",
    "Rhizopus","Saccharomyces","Schizophyllum","Talaromyces","Trametes","Trichoderma","Xylaria",
]
_BACTERIA_GENERA = [
    "Actinomadura","Actinoplanes","Amycolatopsis","Bacillus","Brevibacillus","Burkholderia",
    "Clostridium","Corynebacterium","Enterococcus","Erwinia","Escherichia","Frankia",
    "Kitasatospora","Klebsiella","Lactobacillus","Lysobacter","Micromonospora","Mycobacterium",
    "Myxococcus","Nocardia","Nocardiopsis","Paenibacillus","Photorhabdus","Pseudomonas",
    "Rhodococcus","Saccharopolyspora","Salinispora","Serratia","Sorangium","Staphylococcus",
    "Streptococcus","Streptomyces","Streptosporangium","Thermoactinomyces","Vibrio","Xenorhabdus",
]
_MARINE_GENERA = [
    "Aplysina","Axinella","Callyspongia","Cinachyrella","Cliona","Crambe","Discodermia",
    "Dysidea","Halichondria","Haliclona","Hippiospongia","Hymeniacidon","Ircinia","Jaspis",
    "Lissodendoryx","Mycale","Petrosia","Plakortis","Rhopaloeides","Spongia","Stylissa","Theonella",
    "Lyngbya","Moorea","Oscillatoria","Symploca",
]
for g in _PLANT_GENERA: GENUS_TO_TYPE[g] = "Plant"
for g in _FUNGI_GENERA: GENUS_TO_TYPE[g] = "Fungi"
for g in _BACTERIA_GENERA: GENUS_TO_TYPE[g] = "Bacteria"
for g in _MARINE_GENERA: GENUS_TO_TYPE[g] = "Marine"

ORGANISM_TYPE_FALLBACK = [
    (re.compile(r'Plantae|Viridiplantae|Streptophyta|Magnoliopsida|Liliopsida|Pinopsida', re.I), "Plant"),
    (re.compile(r'Fungi|Ascomycota|Basidiomycota|Zygomycota|Chytridiomycota', re.I), "Fungi"),
    (re.compile(r'Bacteria|Actinobacteria|Proteobacteria|Firmicutes|Cyanobacteria|Bacteroidetes', re.I), "Bacteria"),
    (re.compile(r'Animalia|Metazoa|Chordata|Arthropoda|Mollusca|Cnidaria|Porifera|Echinodermata|Annelida', re.I), "Animal"),
]


def inferOrganismType(organisms, family, kingdom, phylum) -> str:
    """物种名→生物分类。
    Genus-level organism type classification.

    Args:
        organisms: 物种名字符串 (可包含多个属名)
        family: 科名
        kingdom: 界名
        phylum: 门名

    Returns:
        str: "Plant", "Fungi", "Bacteria", "Animal", "Marine", or "Unknown"
    """
    orgStr = str(organisms) if organisms and str(organisms) != "nan" else ""
    for word in orgStr.replace(";", " ").replace(",", " ").split():
        cleaned = word.strip().capitalize()
        if cleaned in GENUS_TO_TYPE:
            return GENUS_TO_TYPE[cleaned]
    textPool = " ".join(filter(None, [
        orgStr,
        str(family) if family and str(family) != "nan" else None,
        str(kingdom) if kingdom and str(kingdom) != "nan" else None,
        str(phylum) if phylum and str(phylum) != "nan" else None,
    ]))
    for regex, orgType in ORGANISM_TYPE_FALLBACK:
        if regex.search(textPool):
            return orgType
    return "Unknown"


# =====================================================================
# Bug 5 辅助: 纯糖分子检测 (Pure Sugar Molecule Detection)
# =====================================================================
def _checkPureSugarMolecule(mol) -> bool:
    """判断分子是否为纯糖 (无苷元连接)。
    Check if a molecule is a pure sugar (no aglycone bonds).

    纯糖分子 (如二糖、含 Sulfate 修饰的糖) 没有苷元,
    切割后 glycan == 全分子。
    Pure sugar molecules have no aglycone; cleavage yields glycan == full molecule.
    """
    try:
        units = find_mapped_sugar_units(mol)
        if not units:
            return False
        rootBond, _ = detectAllGlycosidicBonds(mol, units)
        return not bool(rootBond)
    except Exception:
        return False


# =====================================================================
# 分子渲染 (Molecule Rendering)
# =====================================================================
def molToBase64Png(smi, size=(300, 200)) -> str:
    """将 SMILES 渲染为 base64 PNG 字符串。
    Render SMILES to base64-encoded PNG string.

    Args:
        smi: SMILES 字符串
        size: (width, height) 像素

    Returns:
        str: base64 编码的 PNG，空字符串表示失败
    """
    if not smi or str(smi) in ("nan", "None", "", "NULL", "Error"):
        return ""
    if str(smi).startswith("Error"):
        return ""
    try:
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            return ""
        AllChem.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        opts = drawer.drawOptions()
        opts.addStereoAnnotation = True
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return base64.b64encode(drawer.GetDrawingText()).decode("ascii")
    except Exception:
        return ""


def molToHighlightedBase64Png(smi, size=(380, 250)) -> str:
    """三色高亮全分子图 — 糖核心红色, 修饰黄色, 苷元蓝色。
    Three-color highlighted molecule — sugar RED, modifications YELLOW, aglycone BLUE.

    颜色方案 (Color Scheme):
      红 (#FF4444): 糖环原子 + 糖间桥氧 (sugar ring atoms + inter-sugar bridge O)
      黄 (#FFB800): 糖环外直连修饰原子 (exocyclic modification substituents)
      蓝 (#4488FF): 苷元原子 (aglycone = everything else)

    阈值规则 (Threshold):
      ≤10 重原子的连通分量 → 修饰 (黄色)
      >10 重原子 → 苷元 (蓝色)
      与 classify_sugar_parts() 的 MAX_SUBSTITUENT_SIZE=10 保持一致
    """
    if not smi or str(smi) in ("nan", "None", "", "NULL", "Error"):
        return ""
    if str(smi).startswith("Error"):
        return ""
    try:
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            return ""
        AllChem.Compute2DCoords(mol)

        units = find_mapped_sugar_units(mol)
        highlightAtomColors = {}

        RED = (1.0, 0.27, 0.27)       # 糖核心 (sugar core)
        YELLOW = (1.0, 0.72, 0.0)     # 修饰基团 (modifications)
        BLUE = (0.27, 0.53, 1.0)      # 苷元 (aglycone)

        sugarRingAtomSet = set()

        # Pass 1: 标记糖环原子为红色 (Mark sugar ring atoms RED)
        for u in units:
            for idx in u.get("ring_atoms", []):
                sugarRingAtomSet.add(idx)
                highlightAtomColors[idx] = RED

        # Pass 2: 标记糖间桥氧 (Inter-sugar bridge O → RED)
        bridgeOxygens = set()
        for u in units:
            for rIdx in u.get("ring_atoms", []):
                atom = mol.GetAtomWithIdx(rIdx)
                for nbr in atom.GetNeighbors():
                    nIdx = nbr.GetIdx()
                    if nIdx in sugarRingAtomSet:
                        continue
                    if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 2:
                        otherNbrs = [n for n in nbr.GetNeighbors() if n.GetIdx() != rIdx]
                        if otherNbrs and otherNbrs[0].GetIdx() in sugarRingAtomSet:
                            bridgeOxygens.add(nIdx)
                            highlightAtomColors[nIdx] = RED

        # Pass 3: CC 分析分离修饰/苷元 (CC analysis for mod/aglycone separation)
        sugarAndBridge = sugarRingAtomSet | bridgeOxygens
        nonSugarAtoms = set(range(mol.GetNumAtoms())) - sugarAndBridge

        ccVisited = set()
        components = []
        for startIdx in nonSugarAtoms:
            if startIdx in ccVisited:
                continue
            cc = set()
            bfsQueue = [startIdx]
            while bfsQueue:
                cur = bfsQueue.pop(0)
                if cur in ccVisited:
                    continue
                ccVisited.add(cur)
                cc.add(cur)
                curAtom = mol.GetAtomWithIdx(cur)
                for nbr in curAtom.GetNeighbors():
                    nIdx = nbr.GetIdx()
                    if nIdx in nonSugarAtoms and nIdx not in ccVisited:
                        bfsQueue.append(nIdx)
            components.append(cc)

        MAX_HIGHLIGHT_SUBSTITUENT_SIZE = 10
        isPureSugar = _checkPureSugarMolecule(mol)

        if components and not isPureSugar:
            aglyconeAtoms = set()
            modAtoms = set()
            for cc in components:
                if len(cc) <= MAX_HIGHLIGHT_SUBSTITUENT_SIZE:
                    modAtoms |= cc
                else:
                    aglyconeAtoms |= cc
        elif components and isPureSugar:
            aglyconeAtoms = set()
            modAtoms = set()
            for cc in components:
                modAtoms |= cc
        else:
            aglyconeAtoms = set()
            modAtoms = set()

        for idx in modAtoms:
            highlightAtomColors[idx] = YELLOW
        for idx in aglyconeAtoms:
            highlightAtomColors[idx] = BLUE

        if not units:
            for idx in range(mol.GetNumAtoms()):
                highlightAtomColors[idx] = BLUE
        for idx in range(mol.GetNumAtoms()):
            if idx not in highlightAtomColors:
                highlightAtomColors[idx] = BLUE

        highlightAtoms = list(range(mol.GetNumAtoms()))

        highlightBonds = []
        highlightBondColors = {}
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            c1 = highlightAtomColors.get(a1, BLUE)
            c2 = highlightAtomColors.get(a2, BLUE)
            bIdx = bond.GetIdx()
            highlightBonds.append(bIdx)
            if c1 == RED and c2 == RED:
                highlightBondColors[bIdx] = RED
            elif c1 == YELLOW or c2 == YELLOW:
                highlightBondColors[bIdx] = YELLOW
            elif c1 == RED or c2 == RED:
                highlightBondColors[bIdx] = RED
            else:
                highlightBondColors[bIdx] = BLUE

        drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        opts = drawer.drawOptions()
        opts.addStereoAnnotation = True
        drawer.DrawMolecule(
            mol,
            highlightAtoms=highlightAtoms,
            highlightAtomColors=highlightAtomColors,
            highlightBonds=highlightBonds,
            highlightBondColors=highlightBondColors,
        )
        drawer.FinishDrawing()
        return base64.b64encode(drawer.GetDrawingText()).decode("ascii")
    except Exception:
        return molToBase64Png(smi, size)
