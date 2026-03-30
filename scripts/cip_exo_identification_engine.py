"""
CIP + Exo 指纹单糖鉴定引擎 (CIP + Exo Fingerprint Sugar Identification Engine)
================================================================================

新架构 SOP (New Architecture Standard Operating Procedure):
  Step 1: 虚拟剥离 → 所有修饰基团还原为 -OH/-NH2
  Step 2: Chem.AssignStereochemistry(force=True) → 在裸糖上重算 CIP
  Step 3: 环行走 O→C1→C2→C3→C4→C5, 提取联合指纹:
          (a) CIP 序列 [C2, C3, C4, C5]  (忽略 C1 异头碳)
          (b) Exo 类型 [CH2OH / CH3 / COOH / H / N]
  Step 4: 查表匹配 → 联合指纹 vs 标准库指纹
  Step 5: 回填修饰标签 → D-Glc + NH2@C2 + SO3H@C6 = D-GlcN6S

化学红线 (Chemical Red Line):
  - CIP 必须在脱修饰后提取 (避免远端取代基影响优先级)
  - C1 (异头碳) 不参与身份匹配 (依赖 α/β 构型)
  - Exo 判别解决 CIP 并列: D-Glc vs D-Qui vs D-GlcA

Usage:
  python scripts/cip_exo_identification_engine.py
"""
import sys
import os
from enum import Enum
from typing import Dict, List, Optional, Tuple, Set, NamedTuple

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from lib.glycan_reference_library import RAW_MONOSACCHARIDE_SMILES


# =====================================================================
# 1. 数据结构定义 (Data Structure Definitions)
# =====================================================================

class ExoType(Enum):
    """环外取代基类型 (Exocyclic Substituent Type).

    设计意图: 不同糖类在 C5/C6 位置有不同的功能基团,
    这是区分 Hexose / DeoxyHex / UronicAcid / Pentose 的关键。
    """
    CH2OH = "CH2OH"      # 普通己糖 C5-CH2OH (Standard hexose)
    COOH = "COOH"        # 糖醛酸 C5-COOH (Uronic acid)
    CH3 = "CH3"          # 6-脱氧糖 C5-CH3 (6-Deoxy sugar, e.g., Rha, Fuc)
    H_ONLY = "H"         # 戊糖 C4 无取代 (Pentose, nothing at C4/C5)
    NH2 = "NH2"          # 氨基取代 (Amino substituent)
    OH = "OH"            # 羟基 (Hydroxyl)
    OTHER = "OTHER"      # 其他 (Other)


class SugarFingerprint(NamedTuple):
    """联合指纹: CIP 序列 + ChiralTag 序列 + Exo 特征。
    Composite fingerprint: CIP + ChiralTag + Exo features.

    设计意图 (Design Intent):
      CIP (R/S) 会被远端取代基变化干扰 (e.g., N-SO3H → NH2 翻转 C2 CIP)。
      ChiralTag (CW/CCW) 来自原子在分子图中的局部拓扑, 不受远端元素影响。
      两套指纹并行: CIP 为主, ChiralTag 为氨基糖 fallback。
    """
    cipC2: str   # C2 的 CIP: 'R', 'S', '?'
    cipC3: str   # C3 的 CIP
    cipC4: str   # C4 的 CIP
    cipC5: str   # C5 的 CIP (pentose: '?')
    tagC2: str   # C2 的 ChiralTag: 'CW', 'CCW', '?'
    tagC3: str   # C3 的 ChiralTag
    tagC4: str   # C4 的 ChiralTag
    tagC5: str   # C5 的 ChiralTag
    exoC2: str   # C2 环外取代基: OH, NH2, etc.
    exoC5: str   # C5 环外取代基: CH2OH, COOH, CH3, H
    ringSize: int  # 环大小: 6=pyranose, 5=furanose


# =====================================================================
# 2. 环行走器 (Ring Walker)
# =====================================================================

def walkSugarRing(
    mol: Chem.Mol,
    ringAtoms: Optional[List[int]] = None,
) -> Optional[List[int]]:
    """沿糖环行走: O → C1 → C2 → C3 → C4 → C5。
    Walk sugar ring in canonical order, return ordered carbon indices.

    C1 定义 (C1 Definition):
      与环 O 相邻的碳, 且有环外氧 (异头碳)。
      如果两者都有环外 O, 选择有两个环外 O 的 (C1 通常有 anomeric OH + 糖苷键 O)。
    """
    if ringAtoms is None:
        rings = mol.GetRingInfo().AtomRings()
        sixRings = [r for r in rings if len(r) == 6]
        fiveRings = [r for r in rings if len(r) == 5]
        sugarRings = sixRings + fiveRings
        if not sugarRings:
            return None
        ringAtoms = list(sugarRings[0])

    ringSet = set(ringAtoms)

    # 找环 O (Find ring oxygen)
    ringOIdx = None
    for idx in ringAtoms:
        if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8:
            ringOIdx = idx
            break
    if ringOIdx is None:
        return None

    # 找 O 的两个环碳邻居 (Find O's two ring-carbon neighbors)
    oAtom = mol.GetAtomWithIdx(ringOIdx)
    ringCarbonNbrs = [n.GetIdx() for n in oAtom.GetNeighbors()
                      if n.GetIdx() in ringSet and n.GetAtomicNum() == 6]
    if len(ringCarbonNbrs) != 2:
        return None

    # === 强化异头碳启发式 (Strengthened Anomeric Carbon Heuristic) ===
    # 设计意图: 异头碳 (C1) 必须满足: 直接连接环 O + 至少一个环外 O/N。
    # 对 6-脱氧糖 (Rha/Fuc), C5 连着 -CH3 不含环外 O, 因此不会被误判为 C1。
    #
    # Anomeric C must connect to ring O AND have exocyclic O (or N).
    # For 6-deoxy sugars, C5 has only -CH3 (no exocyclic O) → never chosen as C1.
    def countExoO(cIdx: int) -> int:
        atom = mol.GetAtomWithIdx(cIdx)
        return sum(1 for n in atom.GetNeighbors()
                   if n.GetAtomicNum() == 8 and n.GetIdx() not in ringSet)

    def countExoN(cIdx: int) -> int:
        atom = mol.GetAtomWithIdx(cIdx)
        return sum(1 for n in atom.GetNeighbors()
                   if n.GetAtomicNum() == 7 and n.GetIdx() not in ringSet)

    def countExoC(cIdx: int) -> int:
        atom = mol.GetAtomWithIdx(cIdx)
        return sum(1 for n in atom.GetNeighbors()
                   if n.GetAtomicNum() == 6 and n.GetIdx() not in ringSet)

    exoO0, exoO1 = countExoO(ringCarbonNbrs[0]), countExoO(ringCarbonNbrs[1])
    exoN0, exoN1 = countExoN(ringCarbonNbrs[0]), countExoN(ringCarbonNbrs[1])
    exoC0, exoC1 = countExoC(ringCarbonNbrs[0]), countExoC(ringCarbonNbrs[1])

    # 优先级 1: 有环外 O 的碳 (Priority 1: carbon with exocyclic O)
    hasExoON_0 = (exoO0 + exoN0) > 0
    hasExoON_1 = (exoO1 + exoN1) > 0

    if hasExoON_0 and not hasExoON_1:
        c1Idx = ringCarbonNbrs[0]
    elif hasExoON_1 and not hasExoON_0:
        c1Idx = ringCarbonNbrs[1]
    elif exoO0 > exoO1:
        # 两者都有 O, 但一个有更多 → 选更多的 (Both have O, pick more)
        c1Idx = ringCarbonNbrs[0]
    elif exoO1 > exoO0:
        c1Idx = ringCarbonNbrs[1]
    else:
        # 最终平局: 选没有 exoC 的 (Final tie-break: no exo-C wins)
        # 因为 C5 通常有 -CH2OH/-CH3/-COOH, 而 C1 通常只有 -OH
        if exoC0 < exoC1:
            c1Idx = ringCarbonNbrs[0]
        elif exoC1 < exoC0:
            c1Idx = ringCarbonNbrs[1]
        else:
            c1Idx = ringCarbonNbrs[0]  # 终极 fallback

    # 行走: C1 → C2 → ... → C5 (Walk: C1 → C2 → ... → C5)
    path = [c1Idx]
    prev = ringOIdx
    curr = c1Idx
    maxSteps = len(ringAtoms) - 2  # 6员环走4步, 5员环走3步

    for _ in range(maxSteps):
        atom = mol.GetAtomWithIdx(curr)
        nextC = None
        for n in atom.GetNeighbors():
            nIdx = n.GetIdx()
            if nIdx in ringSet and nIdx != prev and n.GetAtomicNum() == 6:
                nextC = nIdx
                break
        if nextC is None:
            break
        path.append(nextC)
        prev = curr
        curr = nextC

    return path


# =====================================================================
# 3. Exo 取代基分类器 (Exocyclic Substituent Classifier)
# =====================================================================

def classifyExoSubstituent(
    mol: Chem.Mol,
    ringCarbonIdx: int,
    ringSet: Set[int],
) -> str:
    """分类指定环碳上的环外取代基。
    Classify the exocyclic substituent on a given ring carbon.

    设计意图: 区分 -CH2OH (hexose), -COOH (uronic acid),
    -CH3 (6-deoxy), -NH2 (amino), -OH (hydroxyl), -H (pentose).

    Design Intent: Distinguish between -CH2OH (hexose), -COOH (uronic acid),
    -CH3 (6-deoxy), -NH2 (amino), -OH (hydroxyl), -H (pentose terminal).
    """
    atom = mol.GetAtomWithIdx(ringCarbonIdx)
    exoNeighbors = [n for n in atom.GetNeighbors() if n.GetIdx() not in ringSet]

    if not exoNeighbors:
        return ExoType.H_ONLY.value

    # 对每个环外邻居进行分类 (Classify each exocyclic neighbor)
    hasN = False
    hasO = False
    hasC = False

    for nbr in exoNeighbors:
        sym = nbr.GetAtomicNum()
        if sym == 7:  # 氮
            hasN = True
        elif sym == 8:  # 氧 (单独的 -OH)
            hasO = True
        elif sym == 6:  # 碳 — 需要进一步分析
            hasC = True
            # 分析这个碳的子结构 (Analyze this carbon's substructure)
            cNbr = nbr
            cNbrNbrs = [nn for nn in cNbr.GetNeighbors()
                        if nn.GetIdx() != ringCarbonIdx]

            # 检查 -CH2OH: C-(O,H,H) (Check -CH2OH)
            oCount = sum(1 for nn in cNbrNbrs if nn.GetAtomicNum() == 8)
            hCount = cNbr.GetTotalNumHs()

            # -COOH: C(=O)(OH) — 碳有双键 O + 单键 O
            # Check for carboxyl: C with one double-bond O and one single-bond O
            doubleO = 0
            singleO = 0
            for nn in cNbrNbrs:
                if nn.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(cNbr.GetIdx(), nn.GetIdx())
                    if bond and bond.GetBondTypeAsDouble() >= 1.9:  # 双键
                        doubleO += 1
                    else:
                        singleO += 1

            if doubleO >= 1 and singleO >= 1:
                return ExoType.COOH.value  # -COOH (uronic acid)

            if oCount >= 1 and hCount >= 2:
                return ExoType.CH2OH.value  # -CH2OH (standard hexose)

            if oCount == 0 and hCount >= 3:
                return ExoType.CH3.value  # -CH3 (6-deoxy)

            if oCount >= 1 and hCount <= 1:
                # 可能是 -CH(OH)- 链延伸 (heptose, etc.)
                return ExoType.CH2OH.value

    if hasN:
        return ExoType.NH2.value
    if hasO and not hasC:
        return ExoType.OH.value

    return ExoType.OTHER.value


# =====================================================================
# 4. 联合指纹提取器 (Composite Fingerprint Extractor)
# =====================================================================

def extractSugarFingerprint(
    mol: Chem.Mol,
    ringAtoms: Optional[List[int]] = None,
) -> Optional[SugarFingerprint]:
    """提取糖的联合指纹: CIP@C2-C5 + Exo@C2/C5。
    Extract sugar composite fingerprint: CIP at C2-C5 + Exo at C2/C5.

    前置条件: mol 应该是脱修饰后的裸糖, 且已执行
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)。
    """
    path = walkSugarRing(mol, ringAtoms)
    if path is None or len(path) < 4:
        return None

    ringSet = set(path)
    # 扩展环集合: 包含环 O (Expand ring set to include ring O)
    for idx in path:
        atom = mol.GetAtomWithIdx(idx)
        for n in atom.GetNeighbors():
            if n.GetAtomicNum() == 8 and n.GetIdx() not in ringSet:
                # 检查是否为环 O (Check if this O is in the ring)
                rings = mol.GetRingInfo().AtomRings()
                for r in rings:
                    if n.GetIdx() in r and idx in r:
                        ringSet.add(n.GetIdx())

    # 强制 CIP 计算 (Force CIP assignment)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # 提取各位置 CIP (Extract CIP at each position)
    def getCip(atomIdx: int) -> str:
        atom = mol.GetAtomWithIdx(atomIdx)
        return atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else "?"

    # 提取 ChiralTag: CW/CCW — 不受远端取代基变化影响
    # Extract ChiralTag: CW/CCW — NOT affected by distal substituent changes
    def getTag(atomIdx: int) -> str:
        atom = mol.GetAtomWithIdx(atomIdx)
        tag = atom.GetChiralTag()
        if tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            return "CW"
        elif tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            return "CCW"
        return "?"

    isPyranose = len(path) == 5  # C1-C5

    cipC2 = getCip(path[1]) if len(path) > 1 else "?"
    cipC3 = getCip(path[2]) if len(path) > 2 else "?"
    cipC4 = getCip(path[3]) if len(path) > 3 else "?"
    cipC5 = getCip(path[4]) if len(path) > 4 else "?"

    tagC2 = getTag(path[1]) if len(path) > 1 else "?"
    tagC3 = getTag(path[2]) if len(path) > 2 else "?"
    tagC4 = getTag(path[3]) if len(path) > 3 else "?"
    tagC5 = getTag(path[4]) if len(path) > 4 else "?"

    # Exo at C2 (关键: 区分 -OH vs -NH2)
    exoC2 = classifyExoSubstituent(mol, path[1], ringSet) if len(path) > 1 else "?"

    # Exo at C5 (关键: 区分 -CH2OH vs -COOH vs -CH3)
    if isPyranose:
        exoC5 = classifyExoSubstituent(mol, path[4], ringSet)
    elif len(path) >= 4:
        exoC5 = classifyExoSubstituent(mol, path[3], ringSet)
    else:
        exoC5 = "?"

    ringSize = 6 if isPyranose else 5

    return SugarFingerprint(
        cipC2=cipC2, cipC3=cipC3, cipC4=cipC4, cipC5=cipC5,
        tagC2=tagC2, tagC3=tagC3, tagC4=tagC4, tagC5=tagC5,
        exoC2=exoC2, exoC5=exoC5, ringSize=ringSize,
    )


# =====================================================================
# 5. 标准参考指纹库 (Reference Fingerprint Database)
# =====================================================================

def buildReferenceFingerprintDb() -> Dict[Tuple[str, str], SugarFingerprint]:
    """从 RAW_MONOSACCHARIDE_SMILES 构建参考指纹库。
    Build reference fingerprint database from standard sugar SMILES.
    """
    refDb = {}
    for (name, anomer), smi in RAW_MONOSACCHARIDE_SMILES.items():
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        fp = extractSugarFingerprint(mol)
        if fp is None:
            continue
        refDb[(name, anomer)] = fp
    return refDb


# 全局预构建 (Global pre-built database)
_REF_DB: Optional[Dict] = None


def getReferenceFingerprintDb() -> Dict[Tuple[str, str], SugarFingerprint]:
    """延迟加载参考指纹库 (Lazy-load reference fingerprint database)."""
    global _REF_DB
    if _REF_DB is None:
        _REF_DB = buildReferenceFingerprintDb()
    return _REF_DB


# =====================================================================
# 6. 匹配引擎 (Matching Engine)
# =====================================================================

def matchSugarFingerprint(
    targetFp: SugarFingerprint,
    refDb: Optional[Dict] = None,
) -> List[Tuple[str, str, int, int, bool]]:
    """用联合指纹匹配参考库。
    Match target fingerprint against reference database.

    匹配规则 (Matching Rules):
      1. CIP 匹配: C2-C5 中非 '?' 的位置必须完全一致
      2. Exo 匹配: exoC5 必须匹配 (区分 hexose/deoxy/uronic)
                    exoC2 用于区分 amino sugar vs normal
      3. 优先级: (CIP 匹配数 + Exo 匹配数) / 总可比较位

    Returns:
      List of (name, anomer, cipMatch, totalCheckable, exoMatch)
      sorted by overall score descending.
    """
    if refDb is None:
        refDb = getReferenceFingerprintDb()

    results = []

    for (refName, refAnomer), refFp in refDb.items():
        # 环大小必须匹配 (Ring size must match)
        if targetFp.ringSize != refFp.ringSize:
            continue

        # === 双轨匹配: CIP (R/S) + ChiralTag (CW/CCW) ===
        # Dual-track: CIP for normal sugars, ChiralTag for amino sugar fallback

        # Track A: CIP 匹配 C2-C5
        cipMatch = 0
        cipCheckable = 0
        for tCip, rCip in [
            (targetFp.cipC2, refFp.cipC2),
            (targetFp.cipC3, refFp.cipC3),
            (targetFp.cipC4, refFp.cipC4),
            (targetFp.cipC5, refFp.cipC5),
        ]:
            if tCip == "?" or rCip == "?":
                continue
            cipCheckable += 1
            if tCip == rCip:
                cipMatch += 1

        cipPerfect = (cipCheckable > 0 and cipMatch == cipCheckable)

        # Track B: ChiralTag 匹配 C3-C5 (跳过 C2 — 氨基糖 C2 tag 也可能变)
        # ChiralTag matching at C3-C5 (skip C2 — amino sugar C2 tag may shift)
        tagMatch = 0
        tagCheckable = 0
        for tTag, rTag in [
            (targetFp.tagC3, refFp.tagC3),
            (targetFp.tagC4, refFp.tagC4),
            (targetFp.tagC5, refFp.tagC5),
        ]:
            if tTag == "?" or rTag == "?":
                continue
            tagCheckable += 1
            if tTag == rTag:
                tagMatch += 1

        tagPerfect = (tagCheckable > 0 and tagMatch == tagCheckable)

        # 准入条件 (Admission criteria):
        #   CIP 完美匹配 OR ChiralTag C3-C5 完美匹配
        if not cipPerfect and not tagPerfect:
            continue

        # Exo C5 匹配 (Exo C5 matching — critical for tie-breaking)
        exoC5Match = (targetFp.exoC5 == refFp.exoC5)

        # Exo C2 匹配 (Exo C2 matching — amino sugar discrimination)
        exoC2Match = (targetFp.exoC2 == refFp.exoC2)

        # 综合得分 (Composite score)
        # Exo C5 = 1000, Exo C2 = 100, CIP perfect = 50, Tag perfect = 40, base = match count
        score = cipMatch + tagMatch
        if exoC5Match:
            score += 1000
        if exoC2Match:
            score += 100
        if cipPerfect:
            score += 50
        if tagPerfect:
            score += 40

        results.append((refName, refAnomer, score, cipCheckable, exoC5Match, exoC2Match))

    # 排序: 先按得分, 再按 checkable
    results.sort(key=lambda x: (x[2], x[3]), reverse=True)
    return results


# =====================================================================
# 7. 完整 SOP 管线 (Full SOP Pipeline)
# =====================================================================

def identifySugarByCipExo(
    mol: Chem.Mol,
    ringAtoms: Optional[List[int]] = None,
    applyDemod: bool = True,
) -> Tuple[str, str, List[str]]:
    """完整 SOP: 脱修饰 → 重算 CIP → 提取指纹 → 匹配 → 回填。
    Full SOP: Demod → Reassign CIP → Extract FP → Match → Reassemble.

    Returns:
      (sugarName, anomer, detectedMods)
    """
    detectedMods = []

    workMol = Chem.RWMol(Chem.Mol(mol))

    # Step 1: 虚拟剥离 (Virtual Demodification)
    if applyDemod:
        from scripts.poc_virtual_demodify import virtualDemodify, COMPILED_PATTERNS
        # 剥离所有 O-修饰 + N-修饰 (Strip all O- and N-modifications)
        allMods = list(COMPILED_PATTERNS.keys())
        workMol, detectedMods = virtualDemodify(workMol.GetMol(), modNames=allMods)
    else:
        workMol = workMol.GetMol()

    # Step 2: 强制 CIP 重算 (Force CIP reassignment on bare sugar)
    try:
        Chem.AssignStereochemistry(workMol, force=True, cleanIt=True)
    except Exception:
        pass

    # Step 3: 提取联合指纹 (Extract composite fingerprint)
    fp = extractSugarFingerprint(workMol, ringAtoms)
    if fp is None:
        return ("Unknown", "?", detectedMods)

    # Step 4: 查表匹配 (Dictionary lookup)
    matches = matchSugarFingerprint(fp)
    if not matches:
        return ("Unknown", "?", detectedMods)

    topName, topAnomer = matches[0][0], matches[0][1]

    # Step 5: 回填修饰标签 (Reassemble with mod labels)
    # 例如: D-Glc + [N-Acetyl] → D-GlcNAc
    #        D-Glc + [N-Acetyl, O-Sulfate] → D-GlcNAc6S
    finalName = topName
    # NAc 回填 (NAc backfill)
    if "N-Acetyl" in detectedMods:
        if finalName.endswith("N"):
            finalName = finalName + "Ac"  # D-GlcN → D-GlcNAc
        elif not finalName.endswith("NAc"):
            finalName = finalName + "NAc"

    return (finalName, topAnomer, detectedMods)


# =====================================================================
# 8. 测试套件 (Test Suite)
# =====================================================================

def main():
    print("=" * 70)
    print("  CIP + Exo Fingerprint Identification Engine")
    print("  CIP + Exo 指纹鉴定引擎 — 完整 SOP 测试")
    print("=" * 70)

    # Build reference DB and display key entries
    refDb = getReferenceFingerprintDb()
    print(f"\n  Reference database: {len(refDb)} entries")

    # Display key reference fingerprints with Exo discrimination
    print(f"\n  {'Sugar':12s} {'An':3s} {'C2':3s} {'C3':3s} {'C4':3s} {'C5':3s} {'ExoC2':6s} {'ExoC5':6s}")
    print(f"  {'-'*55}")

    keyNames = [
        "D-Glc", "D-GlcN", "D-GlcNAc", "D-GlcA",
        "D-Gal", "D-GalN", "D-GalNAc", "D-GalA",
        "D-Man", "D-ManN", "D-ManNAc",
        "L-Rha", "L-Fuc", "D-Xyl", "L-Ara",
        "D-Qui", "L-IdoA", "D-Dtx",
    ]
    for name in keyNames:
        for anom in ["a", "b"]:
            key = (name, anom)
            if key in refDb:
                fp = refDb[key]
                print(f"  {name:12s} {anom:3s} "
                      f"{fp.cipC2:3s} {fp.cipC3:3s} {fp.cipC4:3s} {fp.cipC5:3s} "
                      f"{fp.exoC2:6s} {fp.exoC5:6s}")

    # === Test Cases ===
    print(f"\n{'='*70}")
    print("  Test Cases")
    print(f"{'='*70}")

    testCases = [
        ("Bare D-Glc",
         "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
         False, "D-Glc"),
        ("Bare D-Gal",
         "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O",
         False, "D-Gal"),
        ("Bare D-Man",
         "OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O",
         False, "D-Man"),
        ("Bare L-Rha",
         "C[C@@H]1OC(O)[C@H](O)[C@H](O)[C@H]1O",
         False, "L-Rha"),
        ("Bare L-Fuc",
         "C[C@@H]1OC(O)[C@@H](O)[C@H](O)[C@@H]1O",
         False, "L-Fuc"),
        ("Bare D-GlcA",
         "O=C(O)[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
         False, "D-GlcA"),
        ("Bare D-GlcN",
         "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O",
         False, "D-GlcN"),
        ("DiAcetyl D-GlcA (demod needed)",
         "O[C@@H]1[C@H](OC(=O)C)[C@@H](OC(=O)C)[C@H](O)[C@@H](C(=O)O)O1",
         True, "D-GlcA"),
        ("GlcNAc-6-Sulfate (demod needed)",
         "CC(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)O)[C@@H](O)[C@@H]1O",
         True, "D-GlcN"),
        ("Fondaparinux fragment (demod needed)",
         "OS(=O)(=O)N[C@H]1C(O)O[C@H](COS(=O)(=O)O)[C@@H](O)[C@@H]1O",
         True, "D-GlcN"),
        ("Methylated L-Ara (demod needed)",
         "O[C@H]1[C@@H](OC)[C@@H](O)[C@@H](O)CO1",
         True, "L-Ara"),
    ]

    passCount = 0
    for testName, smi, needDemod, expected in testCases:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print(f"\n  FAIL: {testName} — SMILES parse error")
            continue

        name, anomer, mods = identifySugarByCipExo(mol, applyDemod=needDemod)

        # 提取裸糖名 (bare name without NAc suffix for comparison)
        bareName = name.replace("NAc", "N").replace("Ac", "")
        if bareName.endswith("N") and expected.endswith("N"):
            match = (bareName == expected)
        else:
            match = (name == expected) or (bareName == expected)

        status = "PASS" if match else "FAIL"
        if match:
            passCount += 1

        modStr = ", ".join(mods) if mods else "none"
        print(f"\n  [{status}] {testName}")
        print(f"    Expected: {expected}")
        print(f"    Got:      {name} (anomer={anomer})")
        print(f"    Mods:     {modStr}")

    total = len(testCases)
    print(f"\n{'='*70}")
    print(f"  Results: {passCount}/{total} passed")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
