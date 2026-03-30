"""
GlycoNP 泛指糖精确修复管线 (Generic Sugar Rescue Pipeline)

将 Hex/Pen/dHex/Hept/Oct/Non 等泛化标签修复为精确单糖名称。
Uses 4 strategies in priority order:
  A. Text-Mining from compound name/iupac_name/synonyms
  C. Rule-Based scaffold-specific mapping
  D. Statistical conditional probability from same-scaffold prior
  (B. API retrieval skipped — no PubChem CID column)

使用方法 (Usage):
  python scripts/rescue_generic_sugars.py
"""
import os
import re
import sys
import time
from collections import Counter, defaultdict
from typing import Optional

import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# =====================================================================
# 策略 A: 化合物名称文本挖掘字典 (Text-Mining Dictionary)
# 从化合物名称中提取明确的糖名称线索
# =====================================================================
# 优先级排序: 最长匹配优先, 避免 "glucoside" 误匹配 "galactoside"
TEXT_MINING_RULES = [
    # 精确多字匹配 (Multi-word exact, highest priority)
    # Highest priority: uronic acids & amino sugars
    (r'glucuronic\s*acid|glucuronide|glucuronosyl', 'D-GlcA'),
    (r'galacturonic\s*acid|galacturonide', 'D-GalA'),
    (r'N-acetylglucosamin|GlcNAc', 'D-GlcNAc'),
    (r'N-acetylgalactosamin|GalNAc', 'D-GalNAc'),
    (r'neuraminic|sialic|Neu5Ac|NeuAc', 'Neu5Ac'),

    # 精确单糖名 + 复合词根匹配 (Exact names + compound word roots)
    # v12 增强: 移除 \b 限制, 改用包含匹配, 确保 glucosinolate 等复合词可识别
    # v12 enhanced: removed \b constraints, use containment matching
    # to handle compounds like glucosinolate, galactolipid, mannoprotein
    (r'glucopyranosid|glucosinolat|glucosid|glucosyl|glucofuranos|gluco(?:se)?(?![a-z])', 'D-Glc'),
    (r'galactopyranosid|galactolipid|galactosid|galactosyl|galacto(?:se)?(?![a-z])', 'D-Gal'),
    (r'mannopyranosid|mannoprotein|mannosid|mannosyl|manno(?:se)?(?![a-z])', 'D-Man'),
    (r'xylopyranosid|xylosid|xylosyl|xylan|xylo(?:se)?(?![a-z])', 'D-Xyl'),
    (r'arabinopyranosid|arabinosid|arabinosyl|arabinan|arabino(?:se)?(?![a-z])', 'L-Ara'),
    (r'rhamnopyranosid|rhamnosid|rhamnosyl|rhamno(?:se)?(?![a-z])', 'L-Rha'),
    (r'fucopyranosid|fucosid|fucosyl|fucoidan|fuco(?:se)?(?![a-z])', 'L-Fuc'),
    (r'ribopyranosid|ribosid|ribosyl|ribo(?:se)?(?![a-z])', 'D-Rib'),
    (r'fructopyranosid|fructosid|fructosyl|fructan|fructo(?:se)?(?![a-z])', 'D-Fru'),

    # 脱氧糖特异名 (Deoxy sugar specific names)
    (r'digitalosid|digitalose', 'D-Dig'),
    (r'oleandrosid|oleandrose', 'L-Ole'),
    (r'cymarosid|cymarose', 'D-Cym'),
    (r'thevetosid|thevetose', 'D-The'),
    (r'boivinosid|boivinose', 'D-Boi'),
    (r'quinovo(?:se|sid)', 'D-Qui'),
]
# 预编译正则 (Pre-compile regex)
COMPILED_TEXT_RULES = [
    (re.compile(pattern, re.IGNORECASE), sugar)
    for pattern, sugar in TEXT_MINING_RULES
]


def strategyA_textMining(
    name: Optional[str],
    iupacName: Optional[str],
    synonyms: Optional[str],
) -> Optional[str]:
    """
    策略 A: 从化合物名称/IUPAC名/同义词中文本挖掘精确糖名。
    Strategy A: Text-mine exact sugar name from compound name fields.

    Returns:
        精确糖名 or None (if no match found)
    """
    # 合并所有可用文本线索
    textPool = " ".join(filter(None, [
        str(name) if name and str(name) != "nan" else None,
        str(iupacName) if iupacName and str(iupacName) != "nan" else None,
        str(synonyms) if synonyms and str(synonyms) != "nan" else None,
    ]))

    if not textPool.strip():
        return None

    for regex, sugar in COMPILED_TEXT_RULES:
        if regex.search(textPool):
            return sugar

    return None


# =====================================================================
# 策略 C: 骨架专一性映射字典 (Rule-Based Scaffold Mapping)
# 特定化合物大类天然产物中, 糖链选择高度专一
# =====================================================================
SCAFFOLD_RULES = {
    # 大环内酯类 — 庚糖分支专一
    "Macrolides": {
        "Hept": "L-Mycarose",  # 大环内酯 100% → 庚糖类(红霉素系)
        "Hex": "D-Desosamine",  # 脱氧氨基糖
    },
    # 强心苷类甾体 — 脱氧糖主导
    "Steroids": {
        "Hex": "D-Glc",
    },
}


def strategyC_ruleBasedMapping(
    superclass: str,
    genericSugar: str,
) -> Optional[str]:
    """
    策略 C: 基于骨架大类的专一性字典映射。
    Strategy C: Rule-based mapping from scaffold class.

    Returns:
        精确糖名 or None
    """
    scClass = str(superclass).strip()
    if "(Tanimoto=" in scClass:
        scClass = scClass[:scClass.index("(Tanimoto=")].strip()

    if scClass in SCAFFOLD_RULES:
        mapping = SCAFFOLD_RULES[scClass]
        if genericSugar in mapping:
            return mapping[genericSugar]

    return None


# =====================================================================
# 策略 D: 统计学条件概率 (Statistical Conditional Probability)
# P(精确糖 | 骨架, 泛化标签) — 从已精确匹配的同类数据中推断
# =====================================================================
MIN_EVIDENCE = 5        # 最低证据数, 低于此不敢推断
MIN_CONFIDENCE = 0.50   # 最低置信度阈值


def buildStatisticalPrior(df: pd.DataFrame) -> dict:
    """
    构建条件概率表: P(Sugar_exact | Scaffold, Generic_class)

    从已精确识别的样本中, 按 Murcko_Scaffold 分组统计糖分布。
    Build conditional probability table from precisely matched data.
    """
    print("  [Strategy D] Building statistical prior...")

    # 只取精确匹配行 (没有泛化标签)
    genericPattern = r'\bHex\b|\bPen\b|\bdHex\b|\bHexA\b|\bNon\b|\bOct\b|\bHept\b'
    preciseMask = (
        df["Sugar_Sequence"].notna()
        & (df["Sugar_Sequence"] != "")
        & (~df["Sugar_Sequence"].str.contains(genericPattern, na=False, regex=True))
    )
    preciseDf = df.loc[preciseMask].copy()
    print(f"    Precise samples for training: {len(preciseDf):,}")

    # 提取第一个糖名 (作为最主要的糖)
    def extractFirstPreciseSugar(seq: str) -> Optional[str]:
        tokens = re.findall(
            r'Neu5Ac|Neu5Gc|KDO|'
            r'[DL]-[A-Z][a-z]+[A-Z]?[a-z]*',
            str(seq))
        return tokens[0] if tokens else None

    preciseDf["_primary_sugar"] = preciseDf["Sugar_Sequence"].apply(
        extractFirstPreciseSugar)
    preciseDf = preciseDf[preciseDf["_primary_sugar"].notna()]

    # 按骨架分组统计 P(sugar | scaffold)
    priorTable = {}  # {scaffold: {sugar: count}}
    for scaffold, group in preciseDf.groupby("Murcko_Scaffold"):
        sugarDist = group["_primary_sugar"].value_counts().to_dict()
        priorTable[scaffold] = sugarDist

    # 同时按 Superclass 分组 (更粗粒度的后备)
    classPrior = {}
    for sc, group in preciseDf.groupby("Superclass"):
        scClean = str(sc).strip()
        if "(Tanimoto=" in scClean:
            scClean = scClean[:scClean.index("(Tanimoto=")].strip()
        sugarDist = group["_primary_sugar"].value_counts().to_dict()
        classPrior[scClean] = sugarDist

    print(f"    Scaffold priors: {len(priorTable):,}")
    print(f"    Class priors: {len(classPrior):,}")

    return priorTable, classPrior


def strategyD_statisticalInference(
    scaffold: Optional[str],
    superclass: Optional[str],
    genericSugar: str,
    scaffoldPrior: dict,
    classPrior: dict,
) -> Optional[str]:
    """
    策略 D: 从同骨架/同大类中最高频的精确糖推断。
    Strategy D: Infer from most frequent precise sugar in same scaffold/class.

    Returns:
        (精确糖名, 置信度) or (None, 0)
    """
    # 尝试 1: 骨架级别先验 (最精确)
    if scaffold and str(scaffold) != "nan":
        dist = scaffoldPrior.get(str(scaffold), {})
        total = sum(dist.values())
        if total >= MIN_EVIDENCE:
            topSugar = max(dist, key=dist.get)
            confidence = dist[topSugar] / total
            if confidence >= MIN_CONFIDENCE:
                return f"{topSugar}_predicted"

    # 尝试 2: 大类级别先验 (更粗)
    scClean = str(superclass).strip() if superclass else ""
    if "(Tanimoto=" in scClean:
        scClean = scClean[:scClean.index("(Tanimoto=")].strip()
    if scClean and scClean != "nan":
        dist = classPrior.get(scClean, {})
        total = sum(dist.values())
        if total >= MIN_EVIDENCE:
            topSugar = max(dist, key=dist.get)
            confidence = dist[topSugar] / total
            if confidence >= MIN_CONFIDENCE:
                return f"{topSugar}_predicted"

    return None


# =====================================================================
# 主管线: 逐行修复 (Main Pipeline)
# =====================================================================
def rescueSingleToken(
    token: str,
    name: Optional[str],
    iupacName: Optional[str],
    synonyms: Optional[str],
    superclass: Optional[str],
    scaffold: Optional[str],
    scaffoldPrior: dict,
    classPrior: dict,
    logCounter: dict,
) -> str:
    """
    对单个泛化 token 执行修复, 按优先级 A → C → D 尝试。
    Rescue single generic token via A → C → D cascade.
    """
    # 不是泛化标签, 直接返回
    if token not in ("Hex", "Pen", "dHex", "HexA", "Non", "Oct", "Hept"):
        return token

    # 策略 A: 文本匹配
    resultA = strategyA_textMining(name, iupacName, synonyms)
    if resultA:
        logCounter["A"] += 1
        return resultA

    # 策略 C: 骨架规则
    resultC = strategyC_ruleBasedMapping(str(superclass), token)
    if resultC:
        logCounter["C"] += 1
        return resultC

    # 策略 D: 统计推断
    resultD = strategyD_statisticalInference(
        scaffold, superclass, token, scaffoldPrior, classPrior)
    if resultD:
        logCounter["D"] += 1
        return resultD

    # 全部失败, 保留原标签
    logCounter["MISS"] += 1
    return token


def rescueSequence(
    seq: str,
    name: Optional[str],
    iupacName: Optional[str],
    synonyms: Optional[str],
    superclass: Optional[str],
    scaffold: Optional[str],
    scaffoldPrior: dict,
    classPrior: dict,
    logCounter: dict,
) -> str:
    """
    修复整条 Sugar_Sequence 中的泛化 token。
    Rescue all generic tokens in a Sugar_Sequence string.
    """
    genericPattern = re.compile(r'\b(Hex|Pen|dHex|HexA|Non|Oct|Hept)\b')

    def replaceMatch(m):
        token = m.group(0)
        return rescueSingleToken(
            token, name, iupacName, synonyms,
            superclass, scaffold,
            scaffoldPrior, classPrior, logCounter)

    return genericPattern.sub(replaceMatch, seq)


def main():
    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")
    inputPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full_Cleaned.csv")

    print("=" * 70)
    print("  GlycoNP Generic Sugar Rescue Pipeline")
    print("  泛指糖精确修复管线")
    print("=" * 70)

    t0 = time.time()
    df = pd.read_csv(inputPath, low_memory=False, encoding="utf-8-sig")
    total = len(df)
    print(f"  Loaded: {total:,} rows ({time.time()-t0:.1f}s)")

    # ===== Step 1: 审计 =====
    genericPattern = r'\bHex\b|\bPen\b|\bdHex\b|\bHexA\b|\bNon\b|\bOct\b|\bHept\b'
    targetMask = df["Sugar_Sequence"].str.contains(
        genericPattern, na=False, regex=True)
    targetCount = targetMask.sum()
    print(f"\n  [AUDIT] Rows with generic sugars: {targetCount:,}")

    # 旧分布 (Old distribution)
    oldTokens = []
    for seq in df.loc[targetMask, "Sugar_Sequence"].dropna():
        oldTokens.extend(re.findall(genericPattern, str(seq)))
    oldCounter = Counter(oldTokens)
    print(f"  Generic tokens breakdown:")
    for k, v in oldCounter.most_common():
        print(f"    {k:10s} {v:>6,}")

    # ===== Step 2: 策略评估 =====
    print(f"\n  [EVALUATE] Strategy Assessment:")
    nameAvail = df.loc[targetMask, "name"].notna().sum()
    iupacAvail = df.loc[targetMask, "iupac_name"].notna().sum()
    synoAvail = df.loc[targetMask, "synonyms"].notna().sum()
    scAvail = df.loc[targetMask, "Superclass"].notna().sum()
    scaffAvail = df.loc[targetMask, "Murcko_Scaffold"].notna().sum()

    print(f"    Strategy A (Text-Mining):  {nameAvail:,} names, "
          f"{iupacAvail:,} IUPAC, {synoAvail:,} synonyms → FEASIBLE")
    print(f"    Strategy B (API):          No PubChem_CID column → SKIPPED")
    print(f"    Strategy C (Rule-Based):   {scAvail:,} Superclass → FEASIBLE")
    print(f"    Strategy D (Statistical):  {scaffAvail:,} scaffolds → FEASIBLE")

    # ===== Step 3: 构建统计先验 =====
    scaffoldPrior, classPrior = buildStatisticalPrior(df)

    # ===== Step 4: 执行修复 =====
    print(f"\n  [EXECUTE] Rescuing {targetCount:,} rows...")
    logCounter = {"A": 0, "C": 0, "D": 0, "MISS": 0}

    df["Sugar_Sequence_Before"] = df["Sugar_Sequence"].copy()

    rescued = 0
    for idx in df.index[targetMask]:
        oldSeq = str(df.at[idx, "Sugar_Sequence"])
        newSeq = rescueSequence(
            oldSeq,
            df.at[idx, "name"] if pd.notna(df.at[idx, "name"]) else None,
            df.at[idx, "iupac_name"] if pd.notna(df.at[idx, "iupac_name"]) else None,
            df.at[idx, "synonyms"] if pd.notna(df.at[idx, "synonyms"]) else None,
            df.at[idx, "Superclass"] if pd.notna(df.at[idx, "Superclass"]) else None,
            df.at[idx, "Murcko_Scaffold"] if pd.notna(df.at[idx, "Murcko_Scaffold"]) else None,
            scaffoldPrior, classPrior, logCounter,
        )
        if newSeq != oldSeq:
            df.at[idx, "Sugar_Sequence"] = newSeq
            rescued += 1

    print(f"\n  [RESULT] Rescue Summary:")
    print(f"    Total rows modified: {rescued:,} / {targetCount:,}")
    print(f"    Strategy A (Text-Mining):   {logCounter['A']:>6,} tokens")
    print(f"    Strategy C (Rule-Based):    {logCounter['C']:>6,} tokens")
    print(f"    Strategy D (Statistical):   {logCounter['D']:>6,} tokens")
    print(f"    Unrescued (MISS):           {logCounter['MISS']:>6,} tokens")

    totalTokensRescued = logCounter["A"] + logCounter["C"] + logCounter["D"]
    totalTokens = totalTokensRescued + logCounter["MISS"]
    rescueRate = totalTokensRescued / totalTokens * 100 if totalTokens > 0 else 0
    print(f"    Rescue rate: {rescueRate:.1f}%")

    # ===== 新分布 =====
    newTokens = []
    for seq in df["Sugar_Sequence"].dropna():
        newTokens.extend(re.findall(genericPattern, str(seq)))
    newCounter = Counter(newTokens)

    print(f"\n  [COMPARISON] Before → After:")
    print(f"  {'Token':<12s} {'Before':>8s} {'After':>8s} {'Delta':>8s}")
    print(f"  {'-'*38}")
    for key in ["Hex", "Pen", "dHex", "HexA", "Non", "Oct", "Hept"]:
        old = oldCounter.get(key, 0)
        new = newCounter.get(key, 0)
        delta = new - old
        if old > 0 or new > 0:
            sign = "+" if delta > 0 else ""
            print(f"  {key:<12s} {old:>8,} {new:>8,} {sign}{delta:>7,}")

    # 新增精确糖统计 (Newly assigned precise sugars)
    print(f"\n  [NEW ASSIGNMENTS] Top precise sugars assigned:")
    allNewTokens = []
    for seq in df.loc[targetMask, "Sugar_Sequence"].dropna():
        tokens = re.findall(
            r'Neu5Ac|Neu5Gc|KDO|'
            r'[DL]-[A-Z][a-z]+[A-Z]?[a-z]*(?:_predicted)?',
            str(seq))
        allNewTokens.extend(tokens)
    assignedCounter = Counter(allNewTokens)
    for k, v in assignedCounter.most_common(15):
        print(f"    {k:25s} {v:>6,}")

    # ===== 保存 =====
    df.drop(columns=["Sugar_Sequence_Before"], inplace=True)
    outPath = os.path.join(reportDir, "GlycoNP_Pipeline_Full_Cleaned.csv")
    df.to_csv(outPath, index=False, encoding="utf-8-sig")
    print(f"\n  Updated: {outPath}")
    print(f"  Time: {time.time()-t0:.0f}s")
    print(f"  {'='*70}")


if __name__ == "__main__":
    main()
