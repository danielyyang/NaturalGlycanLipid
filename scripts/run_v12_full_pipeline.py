"""
GlycoNP 全量管线重刷 v12.3 (Full Pipeline Rerun v12.3)
======================================================
v12.3 修复清单 (Round 3 — All Tasks):
  Step 0: 重新切割 Glycan/Aglycon (Re-cleavage for stale data)
  Task 2: 修复空的 Aglycon_SMILES / 错误的 Glycan_SMILES
  Task 4: 多糖苷键逐一检测 (Multi-bond per-sugar detection)
  All v12.2 fixes preserved
[TEST DATA ONLY]
"""
import argparse, os, re, sys, time, json, base64, random
from collections import Counter
from typing import Optional, List, Dict
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
from tqdm import tqdm

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from lib.monosaccharide_identifier import generate_refined_sequence
from lib.glycan_topology import (
    find_mapped_sugar_units,
    get_split_smiles,
    find_glycosidic_linkages,
    get_sugar_units,
)
from lib.feature_extractor import getTopologyScaffoldSmiles

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
REPORT_DIR = os.path.join(BASE_DIR, "reports")
INPUT_CSV = os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched_v12_backup_bond.csv")
OUTPUT_CSV = os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched_v12_Final.csv")
CRITICAL_LOG = os.path.join(REPORT_DIR, "critical_failures.log")

GENERIC_PAT = re.compile(r'\b(Hex|Pen|dHex|HexA|Non|Oct|Hept)\b')

# =====================================================================
# 罕见糖-名称交叉验证 (Rare Sugar Cross-Validation)
# =====================================================================
RARE_SUGARS = {"D-Tal", "L-Tal", "D-All", "L-All", "D-Alt", "L-Alt",
               "D-Gul", "L-Gul", "D-Ido", "L-Ido", "D-TalN", "L-TalN"}

NAME_SUGAR_MAP = [
    (re.compile(r'glucuronic|glucuronide|glucuronosyl', re.I), 'D-GlcA'),
    (re.compile(r'glucosinolat|glucopyranosid|glucosid|glucosyl|glucofuranos|gluco(?:se)?(?![a-z])', re.I), 'D-Glc'),
    (re.compile(r'galactopyranosid|galactolipid|galactosid|galactosyl|galacto(?:se)?(?![a-z])', re.I), 'D-Gal'),
    (re.compile(r'mannopyranosid|mannoprotein|mannosid|mannosyl|manno(?:se)?(?![a-z])', re.I), 'D-Man'),
    (re.compile(r'xylopyranosid|xylosid|xylosyl|xylan|xylo(?:se)?(?![a-z])', re.I), 'D-Xyl'),
    (re.compile(r'arabinopyranosid|arabinosid|arabinosyl|arabinan|arabino(?:se)?(?![a-z])', re.I), 'L-Ara'),
    (re.compile(r'rhamnopyranosid|rhamnosid|rhamnosyl|rhamno(?:se)?(?![a-z])', re.I), 'L-Rha'),
    (re.compile(r'fucopyranosid|fucosid|fucosyl|fucoidan|fuco(?:se)?(?![a-z])', re.I), 'L-Fuc'),
]

def extractSugarFromName(name, iupacName, synonyms, abstract=None):
    """从名称中提取糖线索 / Extract sugar clue from name fields."""
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

def crossValidateRareSugars(seq, nameSugar):
    """替换罕见糖为名称推断的常见糖 / Replace rare sugar with name-inferred common sugar."""
    if not nameSugar or not seq:
        return seq
    for rare in RARE_SUGARS:
        if rare in seq:
            seq = seq.replace(rare, f"{nameSugar}(name_corrected)")
    return seq


# =====================================================================
# Task 4: 糖苷键精准检测 — 仅检测 anomeric carbon (V3 修复)
# Glycosidic Bond Detection — anomeric-carbon-only (V3 Fix)
# =====================================================================
def _findAnomericCarbon(mol, ringAtoms):
    """找到糖环的 anomeric carbon (C1): 与环内 O 相邻的环碳.
    Find anomeric carbon (C1): ring carbon adjacent to the ring oxygen.
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


def detectAllGlycosidicBonds(mol, sugarUnits):
    """V3: 仅从 anomeric carbon (C1) 检测 Sugar→Aglycone 糖苷键.
    V3: Only detect glycosidic bonds from anomeric carbon (C1).

    对于每个糖单元, 找到 C1 (与环 O 相邻的碳), 然后检查 C1 的
    环外邻居是否连接 aglycone. 每条独立链产出 1 个 Bond Detail.

    For each sugar unit, find C1 (carbon adjacent to ring O), then check
    if C1's exocyclic neighbor connects to aglycone. Each independent
    chain produces exactly 1 Bond Detail.

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
    unitIdMap = {}  # atomIdx → unit index (包含环原子 + 环外取代原子)
    for uid, u in enumerate(sugarUnits):
        # 环原子 (Ring atoms)
        for idx in u.get("ring_atoms", []):
            allSugarAtoms.add(idx)
            unitIdMap[idx] = uid
        # position_map 包含 C1-C6 等碳原子 (Exocyclic carbons like C6)
        for idx in u.get("position_map", {}).keys():
            allSugarAtoms.add(idx)
            if idx not in unitIdMap:
                unitIdMap[idx] = uid
        # oxygen_map 包含环外羟基/醚氧 (Exocyclic oxygens)
        for idx in u.get("oxygen_map", {}).keys():
            allSugarAtoms.add(idx)
            if idx not in unitIdMap:
                unitIdMap[idx] = uid
        # 环原子的直接非环邻居也标记为该糖的领地 (防止 C6/OH 漏识别)
        # Tag immediate exocyclic neighbors of ring atoms as belonging to this sugar
        for rIdx in u.get("ring_atoms", []):
            atom = mol.GetAtomWithIdx(rIdx)
            for nbr in atom.GetNeighbors():
                nIdx = nbr.GetIdx()
                if nIdx not in unitIdMap:
                    unitIdMap[nIdx] = uid
                    allSugarAtoms.add(nIdx)

    # 构建糖→糖连接图, 找到每条链的 reducing end
    # Build sugar-sugar graph to find reducing ends per chain
    sugarGraph = {uid: set() for uid in range(len(sugarUnits))}
    for uid, u in enumerate(sugarUnits):
        ringSet = set(u.get("ring_atoms", []))
        for rIdx in ringSet:
            atom = mol.GetAtomWithIdx(rIdx)
            for nbr in atom.GetNeighbors():
                nIdx = nbr.GetIdx()
                if nIdx in ringSet:
                    continue
                # 桥原子 → 看对面是否是另一个糖
                bridge = nbr
                for bNbr in bridge.GetNeighbors():
                    bIdx = bNbr.GetIdx()
                    if bIdx == rIdx:
                        continue
                    if bIdx in unitIdMap and unitIdMap[bIdx] != uid:
                        sugarGraph[uid].add(unitIdMap[bIdx])
                        sugarGraph[unitIdMap[bIdx]].add(uid)

    # 找 reducing-end 糖: 直接连接 Aglycone 的糖
    # Find reducing-end sugars: those with anomeric C connected to aglycone
    bondDetails = []
    rootBondType = ""
    seenReducingEndUnits = set()

    for uid, u in enumerate(sugarUnits):
        ringAtoms = u.get("ring_atoms", [])
        ringSet = set(ringAtoms)
        sugarName = u.get("name", f"Sugar_{uid}")

        # 找到 anomeric carbon (C1)
        c1Idx = _findAnomericCarbon(mol, ringAtoms)
        if c1Idx is None:
            continue

        c1Atom = mol.GetAtomWithIdx(c1Idx)

        # 仅检查 C1 的环外邻居 (Only check exocyclic neighbors of C1)
        for nbr in c1Atom.GetNeighbors():
            nIdx = nbr.GetIdx()
            if nIdx in ringSet:
                continue  # 环内原子, 跳过

            # 确定桥原子类型 (bridge atom type)
            atomicNum = nbr.GetAtomicNum()
            linkType = ""
            if atomicNum == 8:   linkType = "O"
            elif atomicNum == 16: linkType = "S"
            elif atomicNum == 7:  linkType = "N"
            elif atomicNum == 6:  linkType = "C"
            if not linkType:
                continue

            # 桥原子的另一端是否连接 aglycone?
            # Does the other end of the bridge connect to aglycone?
            targetType = "Aglycone"
            for bNbr in nbr.GetNeighbors():
                bIdx = bNbr.GetIdx()
                if bIdx == c1Idx:
                    continue
                if bIdx in unitIdMap and unitIdMap[bIdx] != uid:
                    targetType = "Sugar"
                    break

            if targetType != "Aglycone":
                continue  # 这是 sugar-sugar 键, 跳过

            # ── 修饰基团过滤: 如果"Aglycone"片段很小(≤12 重原子), 则是修饰而非苷元 ──
            # Modification filter: if the "Aglycone" fragment is small (≤12 heavy
            # atoms), it is a modification group (Bz, Ac, sulfate…), not a true
            # aglycone.  True aglycones (terpenoids, flavonoids, steroids) have
            # ≥15 heavy atoms.
            # NOTE: We use ring_atoms only for the stop-set because allSugarAtoms
            # includes exocyclic neighbors (bridge O) which blocks the BFS.
            sugarRingOnly = set()
            for su in sugarUnits:
                sugarRingOnly.update(su.get("ring_atoms", []))
                sugarRingOnly.update(su.get("position_map", {}).keys())

            fragSize = 0
            visited = {c1Idx}  # 不回溯到 C1
            stack = [nIdx]     # 从桥原子开始
            while stack:
                cur = stack.pop()
                if cur in visited:
                    continue
                visited.add(cur)
                if cur in sugarRingOnly and cur != nIdx:
                    continue  # 碰到糖环/骨架原子就停, 但桥原子本身不停
                fragSize += 1
                for adj in mol.GetAtomWithIdx(cur).GetNeighbors():
                    if adj.GetIdx() not in visited:
                        stack.append(adj.GetIdx())
            # Bug 2 修复: S/N/C-linked 糖苷键的苷元通常较小 (如 Sinigrin 仅 ~7 原子),
            # 使用更宽松的阈值以避免误过滤。O-linked 保留原阈值。
            # Bug 2 fix: S/N/C-linked aglycones are often small (e.g. Sinigrin ~7 atoms),
            # use a relaxed threshold to avoid false filtering. O-linked keeps original.
            sizeThreshold = 3 if linkType in ("S", "N", "C") else 12
            if fragSize <= sizeThreshold:
                continue  # 小片段 → 修饰基团, 不是 aglycone

            # 确定 CIP → α/β
            cipCode = ""
            if c1Atom.HasProp("_CIPCode"):
                cipCode = c1Atom.GetProp("_CIPCode")
            anomer = "?"
            if cipCode == "S": anomer = "α"
            elif cipCode == "R": anomer = "β"

            bondTypeStr = f"{anomer}-{linkType}-linked"

            # 避免同一 sugar unit 输出多条 (Avoid duplicate entries per sugar unit)
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
            break  # 每个糖只取 C1 的第一个 aglycone 键

    # ================================================================
    # Phase 2: Inter-Sugar Bond Detection (Bug 4 修复)
    # 仅当 Phase 1 未找到任何 Sugar→Aglycone 键时运行 (纯糖/寡糖)
    # Only runs when Phase 1 found no Sugar→Aglycone bonds (pure sugars)
    # 设计意图: 含苷元分子只需报告 reducing-end 键, 不需要 inter-sugar 键
    # Design: molecules with aglycone only need reducing-end bonds
    # ================================================================
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

                # 桥原子另一端是否是另一个糖? (Does bridge connect to another sugar?)
                targetUid = None
                for bNbr in nbr.GetNeighbors():
                    bIdx = bNbr.GetIdx()
                    if bIdx == c1Idx:
                        continue
                    if bIdx in unitIdMap and unitIdMap[bIdx] != uid:
                        targetUid = unitIdMap[bIdx]
                        break

                if targetUid is None:
                    continue  # 不是 sugar-sugar 键

                # 避免重复: 只记录 uid < targetUid 的方向
                # Avoid duplicates: only record uid < targetUid direction
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
                break  # 每个糖只取 C1 的第一个 inter-sugar 键

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

def inferOrganismType(organisms, family, kingdom, phylum):
    """物种名→生物分类 / Genus-level organism type classification."""
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
    """判断分子是否为纯糖 (检查是否存在糖苷元连接).
    Check if a molecule is a pure sugar (by checking for absence of aglycone bonds).

    纯糖分子 (如二糖、含 Sulfate 修饰的明胶) 没有苷元, 切割后 glycan == 全分子。
    如果有任何 C1 连接的苷元，则返回 False.
    """
    try:
        units = find_mapped_sugar_units(mol)
        if not units:
            return False
        rootBond, _ = detectAllGlycosidicBonds(mol, units)
        # 如果没有确定的糖苷键，则被认为是纯糖（即无苷元）
        return not bool(rootBond)
    except Exception:
        return False


# =====================================================================
# Step 0: 重新切割 Glycan/Aglycon (Re-Cleavage)
# =====================================================================
def step0Recleavage(df):
    """Task 2 修复: 用最新 get_split_smiles() 重新切割所有有糖的分子.
    Task 2 fix: re-run get_split_smiles() with current code for all sugar-containing molecules.
    """
    print("\n" + "=" * 70)
    print("  Step 0: V12.3 Re-cleavage (Glycan/Aglycon SMILES regeneration)")
    print("=" * 70)
    t0 = time.time()

    # 只对含糖行重新切割 (Only re-cleave rows that contain sugars)
    hasSugar = df["contains_sugar"].astype(str).str.lower().isin(["true", "1", "yes"])
    targetIdx = df.index[hasSugar]
    print(f"  Target: {len(targetIdx):,} sugar-containing rows")

    recleavedCount = 0
    errorCount = 0
    noCleavageCount = 0

    for idx in tqdm(targetIdx, desc="  Re-cleave", ncols=80):
        smi = str(df.at[idx, "canonical_smiles"]) if pd.notna(df.at[idx, "canonical_smiles"]) else ""
        if not smi or smi in ("nan", "None", ""):
            continue
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            aglyconSmi, glycanSmi = get_split_smiles(mol)

            # Bug 5: 纯糖分子保留 Glycan_SMILES
            # Bug 5: Pure sugar molecules (disaccharides etc.) keep Glycan_SMILES
            if glycanSmi and glycanSmi == smi:
                if _checkPureSugarMolecule(mol):
                    # 纯糖: 保留全分子作为 Glycan, 无苷元
                    # Pure sugar: keep full molecule as Glycan, no aglycone
                    glycanSmi = smi
                    aglyconSmi = ""
                else:
                    glycanSmi = ""
                    aglyconSmi = ""
                    noCleavageCount += 1

            df.at[idx, "Glycan_SMILES"] = glycanSmi if glycanSmi else ""
            df.at[idx, "Aglycon_SMILES"] = aglyconSmi if aglyconSmi else ""

            # 同时生成拓扑骨架 (Also generate topology scaffold)
            if aglyconSmi:
                df.at[idx, "Generic_Scaffold"] = getTopologyScaffoldSmiles(aglyconSmi)

            recleavedCount += 1
        except Exception:
            errorCount += 1

    elapsed = time.time() - t0
    print(f"\n  Re-cleaved: {recleavedCount:,}")
    print(f"  No-cleavage: {noCleavageCount:,}")
    print(f"  Errors: {errorCount:,}")
    print(f"  Time: {elapsed:.0f}s")
    return df


# =====================================================================
# Step 1: Sugar_Sequence 重匹配 + 交叉验证 + 多糖苷键
# =====================================================================
def step1Rematch(df):
    """V12.3: 重匹配 + 交叉验证 + 多糖苷键逐一检测."""
    print("\n" + "=" * 70)
    print("  Step 1: V12.3 Re-match + Cross-validation + Multi-bond")
    print("=" * 70)
    t0 = time.time()
    seqChanged = 0
    errorCount = 0
    crossCorrected = 0
    hexPreserved = 0
    bondTypeCounter = Counter()

    # 确保 Bond_Detail 列存在 (Ensure Bond_Detail column exists)
    if "Bond_Detail" not in df.columns:
        df["Bond_Detail"] = ""

    for idx in tqdm(df.index, desc="  V12.3 Re-match", ncols=80):
        # 使用全分子做序列生成, 因为 Glycan_SMILES 可能只包含部分糖
        # Use full molecule for sequence generation — Glycan_SMILES may be incomplete
        fullSmi = str(df.at[idx, "canonical_smiles"]) if pd.notna(
            df.at[idx, "canonical_smiles"]) else ""
        if not fullSmi or fullSmi in ("nan", "None", "", "NULL"):
            continue
        try:
            fullMol = Chem.MolFromSmiles(fullSmi)
            if fullMol is None:
                continue
            newSeq, newMods = generate_refined_sequence(fullMol)
            if not newSeq:
                continue

            oldSeq = str(df.at[idx, "Sugar_Sequence"]) if pd.notna(
                df.at[idx, "Sugar_Sequence"]) else ""

            # 防回退 / Anti-regression
            newHasGeneric = bool(GENERIC_PAT.search(newSeq))
            oldHasGeneric = bool(GENERIC_PAT.search(oldSeq)) if oldSeq else True
            if newHasGeneric and not oldHasGeneric and oldSeq:
                hexPreserved += 1
                newSeq = oldSeq

            # 交叉验证 / Cross-validation
            nameVal = df.at[idx, "name"] if pd.notna(df.at[idx, "name"]) else None
            iupacVal = df.at[idx, "iupac_name"] if "iupac_name" in df.columns and pd.notna(df.at[idx, "iupac_name"]) else None
            synVal = df.at[idx, "synonyms"] if "synonyms" in df.columns and pd.notna(df.at[idx, "synonyms"]) else None
            absVal = None
            if "PMC_Abstract" in df.columns and pd.notna(df.at[idx, "PMC_Abstract"]):
                absVal = str(df.at[idx, "PMC_Abstract"])
            if "Crossref_Abstract" in df.columns and pd.notna(df.at[idx, "Crossref_Abstract"]):
                absVal = (absVal or "") + " " + str(df.at[idx, "Crossref_Abstract"])

            nameSugar = extractSugarFromName(nameVal, iupacVal, synVal, absVal)
            origSeq = newSeq
            newSeq = crossValidateRareSugars(newSeq, nameSugar)
            if newSeq != origSeq:
                crossCorrected += 1

            df.at[idx, "Sugar_Sequence"] = newSeq
            if newSeq != oldSeq:
                seqChanged += 1
            if newMods and "Glycan_Modifications" in df.columns:
                df.at[idx, "Glycan_Modifications"] = newMods

            # Task 4: reducing-end 糖苷键检测 (Reducing-end glycosidic bond)
            # 使用同一个 fullMol (全分子)
            if fullMol:
                units = find_mapped_sugar_units(fullMol)
                rootBond, bondJson = detectAllGlycosidicBonds(fullMol, units)
                if rootBond:
                    df.at[idx, "Aglycone_Linkage_Type"] = rootBond
                    bondTypeCounter[rootBond] += 1
                df.at[idx, "Bond_Detail"] = bondJson

        except Exception:
            errorCount += 1

    elapsed = time.time() - t0
    print(f"\n  Results ({elapsed:.0f}s):")
    print(f"    Changed:          {seqChanged:,}")
    print(f"    Cross-corrected:  {crossCorrected:,}")
    print(f"    Hex-preserved:    {hexPreserved:,}")
    print(f"    Errors:           {errorCount:,}")
    print(f"    Bond types:       {dict(bondTypeCounter.most_common(10))}")

    # 物种分类 / Organism Type
    print("\n  Inferring Organism Types (genus lookup)...")
    df["Organism_Type"] = df.apply(
        lambda r: inferOrganismType(
            r.get("organisms"), r.get("LOTUS_Family") or r.get("LOTUS_family"),
            r.get("LOTUS_kingdom"), r.get("LOTUS_phylum")),
        axis=1)
    orgDist = df["Organism_Type"].value_counts()
    for ot, count in orgDist.items():
        print(f"    {ot:<12} {count:>8,}")
    return df


# =====================================================================
# Step 2: NLP Rescue
# =====================================================================
def step2NlpRescue(df):
    """NLP Rescue with enhanced regex."""
    print("\n" + "=" * 70)
    print("  Step 2: NLP Rescue (v12.3)")
    print("=" * 70)
    try:
        from scripts.rescue_generic_sugars import rescueSequence, buildStatisticalPrior
    except ImportError:
        print("  [SKIP] rescue_generic_sugars not importable")
        return df
    t0 = time.time()
    targetMask = df["Sugar_Sequence"].str.contains(GENERIC_PAT, na=False)
    targetCount = targetMask.sum()
    print(f"  Generic rows: {targetCount:,}")
    if targetCount == 0:
        return df
    scaffoldPrior, classPrior = buildStatisticalPrior(df)
    logCounter = {"A": 0, "C": 0, "D": 0, "MISS": 0}
    rescued = 0
    for idx in tqdm(df.index[targetMask], desc="  NLP Rescue", ncols=80):
        oldSeq = str(df.at[idx, "Sugar_Sequence"])
        try:
            newSeq = rescueSequence(
                oldSeq,
                df.at[idx, "name"] if pd.notna(df.at[idx, "name"]) else None,
                df.at[idx, "iupac_name"] if "iupac_name" in df.columns and pd.notna(df.at[idx, "iupac_name"]) else None,
                df.at[idx, "synonyms"] if "synonyms" in df.columns and pd.notna(df.at[idx, "synonyms"]) else None,
                df.at[idx, "Superclass"] if "Superclass" in df.columns and pd.notna(df.at[idx, "Superclass"]) else None,
                df.at[idx, "Murcko_Scaffold"] if "Murcko_Scaffold" in df.columns and pd.notna(df.at[idx, "Murcko_Scaffold"]) else None,
                scaffoldPrior, classPrior, logCounter,
            )
            if newSeq != oldSeq:
                df.at[idx, "Sugar_Sequence"] = newSeq
                rescued += 1
        except Exception:
            pass
    print(f"\n  NLP: {rescued:,} modified / {targetCount:,} target")
    return df


# =====================================================================
# Step 3: Stats
# =====================================================================
def step3Stats(df):
    print("\n" + "=" * 70)
    print("  Step 3: Sugar Distribution")
    print("=" * 70)
    PAT = (r'Neu5Ac|Neu5Gc|KDO|'
           r'[DL]-[A-Z][a-z]+[A-Z]?[a-z]*(?:\([^)]*\))?|'
           r'Hex|dHex|Pen|HexA|HexN|HexNAc|Non|Oct|Hept')
    allTokens = []
    for seq in df["Sugar_Sequence"].dropna():
        allTokens.extend(re.findall(PAT, str(seq)))
    tc = Counter(allTokens)
    total = sum(tc.values())
    for rank, (sugar, count) in enumerate(tc.most_common(20), 1):
        pct = count / total * 100
        print(f"  {rank:2d}. {sugar:<35} {count:>8,} ({pct:5.1f}%)")
    genericCount = sum(v for k, v in tc.items()
                       if k in ("Hex","Pen","dHex","HexA","Non","Oct","Hept"))
    print(f"\n  Generic:        {genericCount:>8,}")
    print(f"  Total:          {total:>8,}")


# =====================================================================
# Step 4: HTML Report (v12.3)
# =====================================================================
def molToBase64Png(smi, size=(300, 200)):
    """从 SMILES 画分子图 / Draw molecule from SMILES as PNG base64.

    处理 * 假原子: 当原始 SMILES 含切割假原子 (*) 无法解析时,
    用 [H] 替换后重试, 确保能渲染出结构图。
    """
    if not smi or str(smi) in ("nan", "None", "", "NULL", "Error"):
        return ""
    if str(smi).startswith("Error"):
        return ""
    try:
        rawSmi = str(smi)
        mol = Chem.MolFromSmiles(rawSmi)
        # 回退: 替换 * 假原子为 [H] (Fallback: replace dummy atoms)
        if mol is None and "*" in rawSmi:
            sanitized = rawSmi.replace("*", "[H]")
            # 移除断裂片段 (Remove disconnected [H] fragments)
            sanitized = re.sub(r"\.\[H\]", "", sanitized)
            sanitized = re.sub(r"^\[H\]\.", "", sanitized)
            mol = Chem.MolFromSmiles(sanitized)
        if mol is None:
            return ""
        drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        opts = drawer.drawOptions()
        opts.addStereoAnnotation = True
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        return base64.b64encode(drawer.GetDrawingText()).decode("ascii")
    except Exception:
        return ""


def molToHighlightedBase64Png(smi, size=(380, 250)):
    """Bug 6: 三色高亮全分子图 — 糖核心红色, 修饰蓝色, 苷元不着色。
    Bug 6: 三色高亮全分子图 — 糖核心红色, 修饰黄色, 苷元蓝色。
    Bug 6: Three-color highlighted molecule — sugar RED, modifications YELLOW, aglycone BLUE.

    颜色方案 (Color Scheme):
      红 (#FF4444): 糖环原子 + 糖间桥氧 (sugar ring atoms + inter-sugar bridge O)
      黄 (#FFB800): 糖环外直连修饰原子 (exocyclic modification substituents)
      蓝 (#4488FF): 苷元原子 (aglycone = everything else)
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
        modAtomSet = set()

        # Pass 1: 标记糖环原子为红色 (Mark sugar ring atoms RED)
        for u in units:
            for idx in u.get("ring_atoms", []):
                sugarRingAtomSet.add(idx)
                highlightAtomColors[idx] = RED

        # Pass 2: BFS 洪泛标记修饰域 (Flood-fill exocyclic modification chains)
        # 策略: 从分子图中去掉糖环原子和桥氧, 在剩余原子中做连通分量分析
        # 最大的连通分量 = 苷元 (蓝色), 其余小分量 = 修饰 (黄色)
        # Strategy: Remove sugar ring atoms + bridge O from graph, do CC analysis
        # Largest CC = aglycone (BLUE), smaller CCs = modifications (YELLOW)

        # 首先标记所有糖间桥氧 (Inter-sugar bridge O → RED)
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

        # 在非糖原子上做 BFS 连通分量 (CC analysis on non-sugar atoms)
        sugarAndBridge = sugarRingAtomSet | bridgeOxygens
        nonSugarAtoms = set(range(mol.GetNumAtoms())) - sugarAndBridge

        # BFS 找所有连通分量 (Find all connected components via BFS)
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

        # 最大连通分量 = 苷元, 其余 = 修饰 (除非纯糖分子)
        # Largest CC = aglycone, rest = mods (unless pure sugar molecule)
        isPureSugar = _checkPureSugarMolecule(mol)

        if components and not isPureSugar:
            # 正常: 有苷元 → 最大分量蓝色, 其余黄色
            components.sort(key=len, reverse=True)
            aglyconeAtoms = components[0]
            modAtoms = set()
            for cc in components[1:]:
                modAtoms |= cc
        elif components and isPureSugar:
            # 纯糖: 所有非糖原子都是修饰 (黄色), 无苷元
            aglyconeAtoms = set()
            modAtoms = set()
            for cc in components:
                modAtoms |= cc
        else:
            aglyconeAtoms = set()
            modAtoms = set()

        # 标记颜色 (Assign colors)
        for idx in modAtoms:
            highlightAtomColors[idx] = YELLOW
        for idx in aglyconeAtoms:
            highlightAtomColors[idx] = BLUE

        # 无糖环 → 整个分子蓝色 (No sugar rings → all blue)
        if not units:
            for idx in range(mol.GetNumAtoms()):
                highlightAtomColors[idx] = BLUE
        # 兜底: 未分类原子 → 蓝色
        for idx in range(mol.GetNumAtoms()):
            if idx not in highlightAtomColors:
                highlightAtomColors[idx] = BLUE

        highlightAtoms = list(range(mol.GetNumAtoms()))

        # 高亮键 (Highlight bonds)
        highlightBonds = []
        highlightBondColors = {}
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            c1 = highlightAtomColors.get(a1, BLUE)
            c2 = highlightAtomColors.get(a2, BLUE)
            bIdx = bond.GetIdx()
            highlightBonds.append(bIdx)
            # 键颜色策略: 两端都红→红; 含黄→黄; 不然→蓝
            if c1 == RED and c2 == RED:
                highlightBondColors[bIdx] = RED
            elif c1 == YELLOW or c2 == YELLOW:
                highlightBondColors[bIdx] = YELLOW
            elif c1 == RED or c2 == RED:
                # 糖环到苷元/修饰过渡键
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
        # 降级到普通画图 (Fallback to plain drawing)
        return molToBase64Png(smi, size)


def step4HtmlReport(df, sampleSize=1000):
    """v12.3: 独立 SMILES 画图 + Bond_Detail + Topology Scaffold."""
    print("\n" + "=" * 70)
    print(f"  Step 4: {sampleSize}-Molecule HTML Report (v12.3)")
    print("=" * 70)

    random.seed(42)
    sampleIdx = random.sample(list(df.index), min(sampleSize, len(df)))
    t0 = time.time()

    html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<title>GlycoNP {sampleSize} Debug Report v12.3</title>
<style>
  * {{ margin:0; padding:0; box-sizing:border-box; }}
  body {{ font-family:'Segoe UI',system-ui,sans-serif; background:#0a0e1a; color:#e0e0e0; padding:8px; }}
  h1 {{ text-align:center; color:#7ec8e3; font-size:1.5em; margin:8px 0; }}
  .meta {{ text-align:center; color:#6b7280; font-size:0.75em; margin-bottom:6px; }}
  table {{ width:100%; border-collapse:collapse; font-size:0.7em; table-layout:fixed; }}
  th {{ background:#1a2332; color:#7ec8e3; padding:5px 3px; border:1px solid #2a3a4a;
    position:sticky; top:0; z-index:10; text-align:center; font-size:0.85em; }}
  td {{ padding:3px 2px; border:1px solid #1a2332; vertical-align:top;
    word-break:break-word; background:#111827; overflow:hidden; }}
  tr:hover td {{ background:#1a2332; }}
  img {{ border-radius:2px; border:1px solid #2a3a4a; background:white;
    max-width:100%; height:auto; display:block; }}
  .seq {{ font-family:'Courier New',monospace; font-size:0.82em; color:#a0d0f0; word-break:break-all; }}
  .name-corrected {{ color:#f97316; font-weight:bold; }}
  .generic {{ color:#f87171; }}
  .bond-tag {{ display:inline-block; padding:1px 4px; border-radius:3px; font-size:0.8em; margin:1px 0; }}
  .bond-O {{ background:#2563eb33; color:#60a5fa; border:1px solid #2563eb; }}
  .bond-S {{ background:#d9770633; color:#fbbf24; border:1px solid #d97706; }}
  .bond-N {{ background:#7c3aed33; color:#c084fc; border:1px solid #7c3aed; }}
  .bond-C {{ background:#05966933; color:#34d399; border:1px solid #059669; }}
  .org-Plant {{ color:#4ade80; }} .org-Fungi {{ color:#c084fc; }}
  .org-Bacteria {{ color:#f87171; }} .org-Marine {{ color:#38bdf8; }}
  .org-Animal {{ color:#fbbf24; }} .org-Unknown {{ color:#6b7280; }}
  .doi {{ color:#60a5fa; font-size:0.78em; word-break:break-all; }}
  .bio {{ color:#34d399; font-size:0.78em; }}
  .no-cleavage {{ color:#f87171; font-style:italic; font-size:0.8em; }}
  .mod-tag {{ display:inline-block; padding:1px 4px; border-radius:3px; font-size:0.78em;
    margin:1px 0; background:#7c3aed22; color:#c084fc; border:1px solid #7c3aed44; }}
  col.c1 {{ width:6%; }} col.c2 {{ width:11%; }} col.c3 {{ width:9%; }}
  col.c4 {{ width:9%; }} col.c5 {{ width:10%; }} col.c6 {{ width:5%; }}
  col.c7 {{ width:5%; }} col.c8 {{ width:7%; }} col.c9 {{ width:6%; }}
  col.c10 {{ width:3%; }} col.c11 {{ width:3%; }} col.c12 {{ width:9%; }}
  col.c13 {{ width:10%; }}
</style>
</head>
<body>
<h1>GlycoNP {sampleSize} Debug Report v12.3</h1>
<p class="meta">Source: {len(df):,} rows | Generated: {time.strftime('%Y-%m-%d %H:%M')}</p>
"""
    html += '<table>\n<colgroup>'
    for i in range(1, 14):
        html += f'<col class="c{i}">'
    html += '</colgroup>\n<thead><tr>'
    headers = ["ID / Name", "Full Molecule", "Glycan Part", "Aglycone Part",
               "Sugar Sequence", "Bond Detail", "Modification", "Organism",
               "Classification", "Sugars", "Chain", "DOI / References", "Bioactivity"]
    for h in headers:
        html += f'<th>{h}</th>'
    html += '</tr></thead>\n<tbody>\n'

    processed = 0
    for idx in tqdm(sampleIdx, desc="  HTML v12.3", ncols=80):
        row = df.loc[idx]

        fullSmi = str(row.get("canonical_smiles", "")) if pd.notna(row.get("canonical_smiles")) else ""
        glycanSmi = str(row.get("Glycan_SMILES", "")) if pd.notna(row.get("Glycan_SMILES")) else ""
        aglyconSmi = str(row.get("Aglycon_SMILES", "")) if pd.notna(row.get("Aglycon_SMILES")) else ""

        # Bug 6: 全分子三色高亮图 (Full molecule three-color highlighting)
        imgFull = molToHighlightedBase64Png(fullSmi, size=(380, 250))
        imgGlycan = molToBase64Png(glycanSmi, size=(280, 200))
        imgAglycon = molToBase64Png(aglyconSmi, size=(280, 200))

        # Sugar Sequence — 优先显示 Consensus (更准确)
        # Prefer Consensus_Sugar_Sequence (LLM/name-validated) over raw Sugar_Sequence
        seq = str(row.get("Sugar_Sequence", "")) if pd.notna(row.get("Sugar_Sequence")) else ""
        consensusSeq = str(row.get("Consensus_Sugar_Sequence", "")) if pd.notna(row.get("Consensus_Sugar_Sequence")) else ""
        # 如果 Consensus 存在且非空, 优先使用 (Prefer consensus if available)
        displaySeq = consensusSeq if consensusSeq and consensusSeq not in ("", "nan", "None") else seq
        seqHtml = displaySeq
        seqHtml = seqHtml.replace("(name_corrected)", '<span class="name-corrected">(名称纠正)</span>')
        for g in ("Hex", "Pen", "dHex", "HexA"):
            if g in seqHtml:
                seqHtml = seqHtml.replace(g, f'<span class="generic">{g}</span>', 1)

        # Bond Detail — 仅显示 reducing-end sugar → Aglycone
        # 从显示序列中提取每条链的末端糖名 (Extract reducing-end sugars from each chain)
        # 策略: 糖链的最后一个糖名是与苷元直连的 reducing-end sugar
        # Strategy: The last sugar name in each semicolon-separated chain is the reducing-end
        import re as _re
        reducingEndSugars = set()
        for chain in displaySeq.split(";"):
            chain = chain.strip()
            if not chain:
                continue
            # 去掉分支标记 [] 内容 (Remove branch markers)
            chainClean = _re.sub(r"\[.*?\]-?", "", chain)
            # 提取所有糖名 (Extract all sugar names)
            sugarNames = _re.findall(r"(?<![A-Za-z])([DL]-(?:d|me|6d)?[A-Z][a-z]+[A-Za-z]*|Hex|Pen|dHex|HexA|Oct)(?![a-z])", chainClean)
            if sugarNames:
                reducingEndSugars.add(sugarNames[-1])  # 最后一个 = reducing end

        bondDetailRaw = str(row.get("Bond_Detail", "")) if pd.notna(row.get("Bond_Detail")) else "[]"
        bondHtml = ""
        try:
            bonds = json.loads(bondDetailRaw) if bondDetailRaw and bondDetailRaw != "nan" else []
            # 两轮过滤: 第一轮尝试 reducing-end 匹配, 失败则回退到全 Aglycone 去重
            # Two-pass: 1st try reducing-end match, fallback to all-Aglycone dedup
            for useReducingFilter in (True, False):
                seenSugars = set()
                candidateHtml = ""
                for bd in bonds:
                    sugar = bd.get("sugar", "?")
                    target = bd.get("target", "?")
                    bond = bd.get("bond", "?")
                    # 过滤: 仅 Aglycone 键 (Filter: Aglycone bonds only)
                    if "Aglycon" not in target and "aglycon" not in target.lower():
                        continue
                    # 第一轮: 过滤 reducing-end (1st pass: reducing-end filter)
                    if useReducingFilter and reducingEndSugars and sugar not in reducingEndSugars:
                        continue
                    # 去重: 每个糖只显示一条键 (Dedup: one bond per sugar)
                    if sugar in seenSugars:
                        continue
                    seenSugars.add(sugar)
                    bondClass = ""
                    if "-O-" in bond: bondClass = "bond-O"
                    elif "-S-" in bond: bondClass = "bond-S"
                    elif "-N-" in bond: bondClass = "bond-N"
                    elif "-C-" in bond: bondClass = "bond-C"
                    candidateHtml += f'<span class="bond-tag {bondClass}">{sugar}→Aglycon: {bond}</span><br>'
                if candidateHtml:
                    # Pass 1: 仅当覆盖所有链时才接受 (Accept only if all chains covered)
                    if useReducingFilter and reducingEndSugars and len(seenSugars) < len(reducingEndSugars):
                        continue  # 部分匹配 → 回退到 Pass 2 (Partial match → fallback)
                    bondHtml = candidateHtml
                    break  # 有结果就用 (Got results, done)
        except Exception:
            pass
        if not bondHtml:
            rootBond = str(row.get("Aglycone_Linkage_Type", "")) if pd.notna(row.get("Aglycone_Linkage_Type")) else ""
            if rootBond:
                bondClass = ""
                if "-O-" in rootBond: bondClass = "bond-O"
                elif "-S-" in rootBond: bondClass = "bond-S"
                elif "-N-" in rootBond: bondClass = "bond-N"
                elif "-C-" in rootBond: bondClass = "bond-C"
                bondHtml = f'<span class="bond-tag {bondClass}">{rootBond}</span>'
            else:
                bondHtml = "—"

        # Organism
        orgType = str(row.get("Organism_Type", "Unknown"))
        organism = str(row.get("organisms", ""))[:35] if pd.notna(row.get("organisms")) else ""
        family = str(row.get("LOTUS_Family", "")) if pd.notna(row.get("LOTUS_Family")) else ""
        if not family:
            family = str(row.get("LOTUS_family", "")) if pd.notna(row.get("LOTUS_family")) else ""
        orgHtml = f'<span class="org-{orgType}"><b>{orgType}</b></span><br>{organism[:25]}<br><small>{family[:20]}</small>'

        # Classification
        superclass = str(row.get("Superclass", "")) if pd.notna(row.get("Superclass")) else ""
        npSuper = str(row.get("np_classifier_superclass", "")) if pd.notna(row.get("np_classifier_superclass")) else ""
        npClass = str(row.get("np_classifier_class", "")) if pd.notna(row.get("np_classifier_class")) else ""
        displayClass = superclass if superclass and superclass != "Unclassified" else npSuper or npClass
        classHtml = f"<b>{displayClass[:22]}</b><br><small>{npClass[:22]}</small>"

        totalSugar = str(row.get("Total_Sugar_Count", "")) if pd.notna(row.get("Total_Sugar_Count")) else ""
        maxChain = str(row.get("Max_Chain_Length", "")) if pd.notna(row.get("Max_Chain_Length")) else ""

        # DOI
        dois = str(row.get("dois", ""))[:80] if pd.notna(row.get("dois")) else ""
        doiHtml = f'<span class="doi">{dois}</span>' if dois else "—"

        # Bioactivity
        bioRaw = str(row.get("bioactivity_summary", "")) if pd.notna(row.get("bioactivity_summary")) else ""
        nlpBio = str(row.get("NLP_Bioactivity_Profile", "")) if pd.notna(row.get("NLP_Bioactivity_Profile")) else ""
        chemblTargets = str(row.get("ChEMBL_Targets", "")) if pd.notna(row.get("ChEMBL_Targets")) else ""
        bioDisplay = bioRaw or nlpBio
        if chemblTargets:
            bioDisplay += f" [{chemblTargets[:30]}]" if bioDisplay else chemblTargets[:30]
        bioHtml = f'<span class="bio">{bioDisplay[:80]}</span>' if bioDisplay else "—"

        # ID
        identifier = str(row.get("identifier", ""))[:15] if pd.notna(row.get("identifier")) else ""
        name = str(row.get("name", ""))[:35] if pd.notna(row.get("name")) else ""
        nameHtml = f"<b>{identifier}</b><br><small>{name}</small>"

        imgFullHtml = f'<img src="data:image/png;base64,{imgFull}">' if imgFull else "—"

        # Glycan/Aglycon: 防御性画图 — 空则显示标记
        if imgGlycan:
            imgGlycanHtml = f'<img src="data:image/png;base64,{imgGlycan}">'
        elif glycanSmi:
            imgGlycanHtml = '<span class="no-cleavage">Parse Error</span>'
        else:
            imgGlycanHtml = '<span class="no-cleavage">No Cleavage</span>'

        if imgAglycon:
            imgAglyconHtml = f'<img src="data:image/png;base64,{imgAglycon}">'
        elif aglyconSmi:
            imgAglyconHtml = '<span class="no-cleavage">Parse Error</span>'
        else:
            imgAglyconHtml = '<span class="no-cleavage">No Aglycone</span>'

        # Bug 6: Modification 列 — 从 Sugar_Sequence 的修饰注释或 Glycan_Modifications 提取
        # Bug 6: Modification column — from per-unit modification data
        modHtml = ""
        glycanMods = str(row.get("Glycan_Modifications", "")) if pd.notna(row.get("Glycan_Modifications")) else ""
        if glycanMods and glycanMods not in ("nan", "None", "", "NULL"):
            # 解析修饰注释格式: "Sugar_1(*O-Ac,*NAc)" → 提取标签
            # Parse modification annotation format
            import re as _re
            modTokens = _re.findall(r'\*([A-Za-z0-9_-]+)', glycanMods)
            seen = set()
            for mt in modTokens:
                if mt not in seen:
                    modHtml += f'<span class="mod-tag">{mt}</span> '
                    seen.add(mt)
        if not modHtml:
            # 从全分子直接扫描修饰
            if fullSmi:
                try:
                    fullMolTmp = Chem.MolFromSmiles(fullSmi)
                    if fullMolTmp:
                        tmpUnits = find_mapped_sugar_units(fullMolTmp)
                        seen = set()
                        for u in tmpUnits:
                            for m in u.get("modifications", []):
                                if m not in seen:
                                    modHtml += f'<span class="mod-tag">{m}</span> '
                                    seen.add(m)
                except Exception:
                    pass
        if not modHtml:
            modHtml = "—"

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
        processed += 1

    html += '</tbody>\n</table>\n</body>\n</html>'
    outPath = os.path.join(REPORT_DIR, "GlycoNP_1000mol_debug_report.html")
    with open(outPath, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"\n  HTML: {processed:,} rows ({time.time()-t0:.0f}s)")
    print(f"  Output: {outPath}")


# =====================================================================
# Main
# =====================================================================
def main():
    print("=" * 70)
    print("  GlycoNP V12.3 Full Pipeline (Tasks 2-4)")
    print("=" * 70)
    t0 = time.time()

    print(f"\n  Loading: {INPUT_CSV}")
    df = pd.read_csv(INPUT_CSV, low_memory=False)
    print(f"  Rows: {len(df):,}, Columns: {len(df.columns)}")

    df = step0Recleavage(df)
    df = step1Rematch(df)
    df = step2NlpRescue(df)
    step3Stats(df)

    print(f"\n  Saving: {OUTPUT_CSV}")
    df.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")
    print(f"  Saved: {len(df):,} rows")

    step4HtmlReport(df, sampleSize=1000)

    print(f"\n  Total: {time.time()-t0:.0f}s ({(time.time()-t0)/60:.1f}min)")
    print("=" * 70)
    print("  V12.4 PIPELINE COMPLETE!")
    print("=" * 70)


if __name__ == "__main__":
    main()
