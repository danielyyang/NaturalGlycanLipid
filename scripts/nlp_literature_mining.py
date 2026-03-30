"""
Local NLP Literature Mining Engine
===================================

纯本地 NLP 文本挖掘引擎: 无需外部 API,
利用增强版正则字典从 PMC Title + Abstract 中提取:
  - 精确糖名 (Sugar Name Resolution)
  - 生物活性数据 (Bioactivity Extraction)
  - 靶点/疾病 (Target/Disease Extraction)

然后执行零污染安全合并 (Pandas 防火墙) 和 HTML 报告生成。

Usage:
  python scripts/nlp_literature_mining.py
  python scripts/nlp_literature_mining.py --full   # 全量 PMC 抓取 + 挖掘
"""
import argparse
import html
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from typing import Dict, List, Optional, Tuple

import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
REPORT_DIR = os.path.join(BASE_DIR, "reports")
INPUT_CSV = os.path.join(REPORT_DIR, "GlycoNP_LLM_Rescued.csv")
OUTPUT_CSV = os.path.join(REPORT_DIR, "GlycoNP_NLP_Enriched.csv")

GENERIC_PATTERN = re.compile(r'\b(Hex|Pen|dHex|HexA|Non|Oct|Hept)\b')


# =====================================================================
# 增强版糖名正则字典 (Enhanced Sugar Name Dictionary)
# =====================================================================
# 设计意图: 比 rescue_generic_sugars.py 更深度地挖掘文献文本,
# 包含 IUPAC 系统名、通俗名、拉丁形态等全方位匹配

SUGAR_TEXT_RULES = [
    # === 醛糖酸 / 氨基糖 (Uronic acids / Amino sugars, highest priority) ===
    (re.compile(r'glucuronic\s*acid|glucuronide|glucuronosyl|GlcUA|GlcA',
                re.I), 'D-GlcA'),
    (re.compile(r'galacturonic\s*acid|galacturonide|GalUA|GalA', re.I),
     'D-GalA'),
    (re.compile(r'N-acetyl[-\s]*glucosamin|GlcNAc|acetylglucosamine', re.I),
     'D-GlcNAc'),
    (re.compile(r'N-acetyl[-\s]*galactosamin|GalNAc|acetylgalactosamine',
                re.I), 'D-GalNAc'),
    (re.compile(r'neuraminic|sialic\s*acid|Neu5Ac|NeuAc|NANA', re.I),
     'Neu5Ac'),

    # === 六碳糖 (Hexoses) ===
    (re.compile(r'glucopyranosid|glucosid[ey]|glucosyl|\bgluco(?:se)?\b|'
                r'\bGlc[^N]|\bGlc$|\b[Bb]eta-[Dd]-glucosid|\b[Dd]-glucos',
                re.I), 'D-Glc'),
    (re.compile(r'galactopyranosid|galactosid[ey]|galactosyl|'
                r'\bgalacto(?:se)?\b|\bGal[^ANnU]|\bGal$', re.I), 'D-Gal'),
    (re.compile(r'mannopyranosid|mannosid[ey]|mannosyl|\bmanno(?:se)?\b',
                re.I), 'D-Man'),

    # === 五碳糖 (Pentoses) ===
    (re.compile(r'xylopyranosid|xylosid[ey]|xylosyl|\bxylo(?:se)?\b|'
                r'\bXyl\b', re.I), 'D-Xyl'),
    (re.compile(r'arabinopyranosid|arabinosid[ey]|arabinosyl|arabinoside|'
                r'\barabino(?:se)?\b|\bAra\b', re.I), 'L-Ara'),
    (re.compile(r'ribopyranosid|ribosid[ey]|\bribo(?:se)?\b|\bRib\b',
                re.I), 'D-Rib'),
    (re.compile(r'apiopyranosid|apiosid[ey]|apiosyl|\bapio(?:se)?\b|'
                r'\bApi\b', re.I), 'D-Api'),

    # === 脱氧糖 (Deoxysugars) ===
    (re.compile(r'rhamnopyranosid|rhamnosid[ey]|rhamnosyl|\brhamno(?:se)?\b|'
                r'\bRha\b', re.I), 'L-Rha'),
    (re.compile(r'fucopyranosid|fucosid[ey]|fucosyl|\bfuco(?:se)?\b|'
                r'\bFuc\b', re.I), 'L-Fuc'),
    (re.compile(r'quinovopyranosid|quinovosid[ey]|\bquinovo(?:se)?\b',
                re.I), 'D-Qui'),
    (re.compile(r'digitalosid[ey]|\bdigitalo(?:se)?\b', re.I),
     'D-Dig'),
    (re.compile(r'oleandrosid[ey]|\boleandro(?:se)?\b', re.I),
     'L-Ole'),
    (re.compile(r'cymarosid[ey]|\bcymaro(?:se)?\b', re.I), 'D-Cym'),
    (re.compile(r'thevetos(?:id[ey]|e)', re.I), 'D-The'),
    (re.compile(r'boivinosid[ey]|\bboivino(?:se)?\b', re.I),
     'D-Boi'),

    # === 酮糖 (Ketoses) ===
    (re.compile(r'fructopyranosid|fructosid[ey]|fructosyl|fructofuranosid|'
                r'\bfructo(?:se)?\b|\bFru\b', re.I), 'D-Fru'),

    # === 氨基酮糖/特殊 (Amino/Special sugars) ===
    (re.compile(r'glucosamin[ey]|GlcN\b', re.I), 'D-GlcN'),
    (re.compile(r'galactosamin[ey]|GalN\b', re.I), 'D-GalN'),
]


# =====================================================================
# 生物活性关键词字典 (Bioactivity Keyword Dictionary)
# =====================================================================

BIOACTIVITY_RULES = [
    # 细胞毒 / 抗肿瘤 (Cytotoxicity / Antitumor)
    (re.compile(r'cytotoxi[cn]|anti[-\s]?(?:tumor|tumour|cancer|proliferati)|'
                r'antineo(?:plastic)?|tumor\s*inhibit', re.I),
     'Cytotoxic'),
    (re.compile(r'apoptos[ie]s|pro[-\s]?apoptotic|caspase', re.I),
     'Apoptosis-inducing'),

    # 抗菌 / 抗真菌 (Antimicrobial)
    (re.compile(r'anti[-\s]?(?:bacteri|microbi)|bactericid|bacteriostat',
                re.I), 'Antibacterial'),
    (re.compile(r'anti[-\s]?fung|fungicid', re.I), 'Antifungal'),
    (re.compile(r'anti[-\s]?viral|virucid', re.I), 'Antiviral'),
    (re.compile(r'anti[-\s]?(?:malari|plasmod)', re.I), 'Antimalarial'),

    # 抗炎 / 免疫 (Anti-inflammatory / Immune)
    (re.compile(r'anti[-\s]?inflammat|COX[-\s]?[12]|NF[-\s]?kB|TNF|'
                r'interleukin|IL[-\s]?\d', re.I), 'Anti-inflammatory'),
    (re.compile(r'immuno[-\s]?modulat|immuno[-\s]?stimulat|immuno[-\s]?'
                r'suppress', re.I), 'Immunomodulatory'),

    # 抗氧化 (Antioxidant)
    (re.compile(r'anti[-\s]?oxid|radical\s*scaveng|DPPH|ABTS|ORAC|ROS',
                re.I), 'Antioxidant'),

    # 酶抑制 (Enzyme inhibition)
    (re.compile(r'inhibit\w*\s+(?:acetylcholinestera|AChE)', re.I),
     'AChE-inhibitory'),
    (re.compile(r'inhibit\w*\s+(?:alpha[-\s]?glucosidas|amylas)', re.I),
     'alpha-Glucosidase-inh'),
    (re.compile(r'inhibit\w*\s+(?:tyrosinas|PTP)', re.I),
     'Tyrosinase-inh'),

    # 心血管 (Cardiovascular)
    (re.compile(r'cardiotonic|inotropic|cardiac\s*glycosid', re.I),
     'Cardiotonic'),
    (re.compile(r'anti[-\s]?hypertens|vasorelax|ACE\s*inhibit', re.I),
     'Antihypertensive'),

    # 疗效指标 (Quantitative bioactivity markers)
    (re.compile(r'IC\s*50|IC50', re.I), 'IC50-reported'),
    (re.compile(r'MIC\b|minimum\s*inhibitory', re.I), 'MIC-reported'),
    (re.compile(r'EC\s*50|EC50', re.I), 'EC50-reported'),
    (re.compile(r'LD\s*50|LD50', re.I), 'LD50-reported'),
    (re.compile(r'Ki\s*=|Ki\s*value', re.I), 'Ki-reported'),

    # 其他 (Others)
    (re.compile(r'anti[-\s]?diabet|hypoglycemi', re.I), 'Antidiabetic'),
    (re.compile(r'neuro[-\s]?protect|neuroprotect', re.I), 'Neuroprotective'),
    (re.compile(r'hepato[-\s]?protect', re.I), 'Hepatoprotective'),
    (re.compile(r'anti[-\s]?leishman', re.I), 'Antileishmanial'),
    (re.compile(r'anti[-\s]?trypanosom', re.I), 'Antitrypanosomal'),
    (re.compile(r'anti[-\s]?helminth|anthelmint', re.I), 'Anthelmintic'),
    (re.compile(r'wound\s*heal', re.I), 'Wound-healing'),
]

# 靶点/细胞系/病原体字典 (Target/Cell line/Pathogen dictionary)
TARGET_RULES = [
    # 肿瘤细胞系 (Cancer cell lines)
    (re.compile(r'\b(?:MCF[-\s]?7|MDA[-\s]?MB[-\s]?231)\b'), 'MCF-7/MDA-MB'),
    (re.compile(r'\bA549\b'), 'A549'),
    (re.compile(r'\bHeLa\b', re.I), 'HeLa'),
    (re.compile(r'\bHepG2\b'), 'HepG2'),
    (re.compile(r'\bHL[-\s]?60\b'), 'HL-60'),
    (re.compile(r'\bHCT[-\s]?116\b'), 'HCT-116'),
    (re.compile(r'\bHT[-\s]?29\b'), 'HT-29'),
    (re.compile(r'\bK562\b'), 'K562'),
    (re.compile(r'\bPC[-\s]?3\b'), 'PC-3'),
    (re.compile(r'\bSK[-\s]?BR\b'), 'SK-BR'),

    # 病原体 (Pathogens)
    (re.compile(r'Staphylococcus\s*aureus|S\.\s*aureus|MRSA', re.I),
     'S.aureus'),
    (re.compile(r'Escherichia\s*coli|E\.\s*coli', re.I), 'E.coli'),
    (re.compile(r'Plasmodium\s*falciparum|P\.\s*falciparum', re.I),
     'P.falciparum'),
    (re.compile(r'Mycobacterium\s*tuberculosis|M\.\s*tuberculosis', re.I),
     'M.tuberculosis'),
    (re.compile(r'Candida\s*albicans|C\.\s*albicans', re.I), 'C.albicans'),
    (re.compile(r'Trypanosoma|T\.\s*(?:brucei|cruzi)', re.I), 'Trypanosoma'),
    (re.compile(r'Leishmania', re.I), 'Leishmania'),

    # 蛋白靶点 (Protein targets)
    (re.compile(r'\bCOX[-\s]?[12]\b'), 'COX-1/2'),
    (re.compile(r'\bNF[-\s]?[kK]B\b'), 'NF-kB'),
    (re.compile(r'\bTNF[-\s]?alpha\b', re.I), 'TNF-alpha'),
    (re.compile(r'\bHIV\b'), 'HIV'),
]


# =====================================================================
# NLP 文本挖掘引擎 (NLP Text Mining Engine)
# =====================================================================

def mineText(
    title: str,
    abstract: str,
) -> Dict[str, any]:
    """
    从 Title+Abstract 中提取糖名、活性、靶点。
    Mine sugar names, bioactivities, and targets from Title + Abstract.

    返回:
      {
        "sugars": ["D-Glc", "L-Rha", ...],  去重有序
        "bioactivities": ["Cytotoxic", ...],
        "targets": ["A549", ...],
      }
    """
    text = f"{title} {abstract}".strip()
    if not text:
        return {"sugars": [], "bioactivities": [], "targets": []}

    # 糖名提取 (Sugar extraction)
    foundSugars = []
    seen = set()
    for regex, sugar in SUGAR_TEXT_RULES:
        if regex.search(text) and sugar not in seen:
            foundSugars.append(sugar)
            seen.add(sugar)

    # 活性提取 (Bioactivity extraction)
    foundBioact = []
    seen = set()
    for regex, label in BIOACTIVITY_RULES:
        if regex.search(text) and label not in seen:
            foundBioact.append(label)
            seen.add(label)

    # 靶点提取 (Target extraction)
    foundTargets = []
    seen = set()
    for regex, label in TARGET_RULES:
        if regex.search(text) and label not in seen:
            foundTargets.append(label)
            seen.add(label)

    return {
        "sugars": foundSugars,
        "bioactivities": foundBioact,
        "targets": foundTargets,
    }


# =====================================================================
# Task 1: Europe PMC 增量抓取 (Incremental PMC Fetch)
# =====================================================================

def fetchEuropePmc(doi: str) -> Dict[str, str]:
    """Europe PMC REST API fetch."""
    encodedDoi = urllib.parse.quote(doi, safe="")
    url = (f"https://www.ebi.ac.uk/europepmc/webservices/rest/search"
           f"?query=DOI:{encodedDoi}&resultType=core&format=json")
    try:
        req = urllib.request.Request(url)
        req.add_header("User-Agent",
                       "GlycoNP/1.0 (mailto:glyconp@research.dev)")
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            results = data.get("resultList", {}).get("result", [])
            if not results:
                return {"title": "", "abstract": ""}
            entry = results[0]
            title = re.sub(r'<[^>]+>', '', str(entry.get("title", "")))
            abstract = re.sub(r'<[^>]+>', '',
                              str(entry.get("abstractText", "")))
            return {"title": title.strip(), "abstract": abstract.strip()}
    except Exception:
        return {"title": "", "abstract": ""}


def runPmcFetch(df: pd.DataFrame, targetIdxs: list) -> pd.DataFrame:
    """增量抓取 PMC 数据 (仅处理未抓取过的)。"""
    print("\n" + "=" * 70)
    print("  PMC Incremental Fetch")
    print("=" * 70)
    t0 = time.time()

    for col in ["PMC_Title", "PMC_Abstract"]:
        if col not in df.columns:
            df[col] = ""

    needFetch = []
    for idx in targetIdxs:
        existing = str(df.at[idx, "PMC_Title"]) if pd.notna(
            df.at[idx, "PMC_Title"]) else ""
        if existing in ("", "nan", "None"):
            needFetch.append(idx)

    print(f"  Already fetched: {len(targetIdxs) - len(needFetch):,}")
    print(f"  Need fetch: {len(needFetch):,}")

    fetched = 0
    for i, idx in enumerate(needFetch):
        doiStr = str(df.at[idx, "dois"]).strip()
        firstDoi = doiStr.split("|")[0].strip().split(";")[0].strip()
        if not firstDoi or firstDoi in ("nan", "None", ""):
            continue
        time.sleep(0.2)
        result = fetchEuropePmc(firstDoi)
        if result["title"]:
            df.at[idx, "PMC_Title"] = result["title"]
            fetched += 1
        if result["abstract"]:
            df.at[idx, "PMC_Abstract"] = result["abstract"][:3000]

        if (i + 1) % 100 == 0:
            print(f"    [{i+1}/{len(needFetch)}] fetched={fetched} "
                  f"({time.time()-t0:.0f}s)")

    print(f"  New titles fetched: {fetched:,}")
    print(f"  Time: {time.time()-t0:.0f}s")
    return df


# =====================================================================
# Task 2: NLP 批量挖掘 (NLP Batch Mining)
# =====================================================================

def runNlpMining(df: pd.DataFrame, targetIdxs: list) -> pd.DataFrame:
    """
    批量 NLP 文本挖掘: 从 PMC Title+Abstract 提取糖名+活性+靶点。
    Batch NLP text mining from PMC Title + Abstract.
    """
    print("\n" + "=" * 70)
    print("  NLP Text Mining Engine")
    print("  NLP 文本挖掘引擎")
    print("=" * 70)
    t0 = time.time()

    # 初始化列 (Initialize columns)
    for col in ["NLP_Inferred_Sugars", "NLP_Bioactivity_Profile",
                 "NLP_Targets"]:
        if col not in df.columns:
            df[col] = ""

    mined = 0
    sugarHits = 0
    bioactHits = 0
    targetHits = 0
    sugarCounter = Counter()  # 跟踪提取到的糖名分布

    for idx in targetIdxs:
        title = str(df.at[idx, "PMC_Title"]) if pd.notna(
            df.at[idx, "PMC_Title"]) else ""
        abstract = str(df.at[idx, "PMC_Abstract"]) if pd.notna(
            df.at[idx, "PMC_Abstract"]) else ""
        if title in ("nan", "None"):
            title = ""
        if abstract in ("nan", "None"):
            abstract = ""

        # 也扫描已有的 name/iupac_name 字段
        name = str(df.at[idx, "name"]) if pd.notna(
            df.at[idx, "name"]) else ""
        if name in ("nan", "None"):
            name = ""
        fullText = f"{title} {abstract} {name}"

        if not fullText.strip():
            continue

        result = mineText(title, abstract + " " + name)

        if result["sugars"]:
            df.at[idx, "NLP_Inferred_Sugars"] = "; ".join(result["sugars"])
            sugarHits += 1
            sugarCounter.update(result["sugars"])

        if result["bioactivities"]:
            df.at[idx, "NLP_Bioactivity_Profile"] = "; ".join(
                result["bioactivities"])
            bioactHits += 1

        if result["targets"]:
            df.at[idx, "NLP_Targets"] = "; ".join(result["targets"])
            targetHits += 1

        mined += 1

    elapsed = time.time() - t0
    print(f"\n  NLP Mining Results:")
    print(f"    Total processed: {mined:,}")
    print(f"    Sugar resolved: {sugarHits:,}")
    print(f"    Bioactivity found: {bioactHits:,}")
    print(f"    Targets found: {targetHits:,}")
    print(f"\n  Extracted sugar distribution:")
    for k, v in sugarCounter.most_common(15):
        print(f"    {k:12s}: {v:,}")
    print(f"  Time: {elapsed:.0f}s")

    return df


# =====================================================================
# 碳数一致性字典 (Carbon-Class Consistency Map)
# =====================================================================
# 设计意图: 每种泛指标签只允许被替换为化学上同碳类的精确糖。
# 违反此规则的替换被视为"跨类污染" — NLP 抓到了文献中其他分子的糖。

SUGAR_CLASS_MAP = {
    'Hex':  ['D-Glc', 'D-Gal', 'D-Man', 'L-Glc', 'L-Gal', 'D-Tal',
             'D-All', 'D-Gul', 'L-Ido', 'D-Fru'],
    'Pen':  ['D-Xyl', 'L-Ara', 'D-Rib', 'D-Api', 'D-Lyx', 'L-Lyx'],
    'dHex': ['L-Rha', 'L-Fuc', 'D-Fuc', 'L-Ole', 'D-Cym', 'D-The',
             'D-Dig', 'D-Boi', 'D-Qui'],
    'HexN': ['D-GlcN', 'D-GalN', 'D-ManN'],
    'HexA': ['D-GlcA', 'D-GalA', 'D-ManA'],
    'Non':  ['Neu5Ac', 'Neu5Gc', 'KDO'],
    'Oct':  [],   # 八碳糖极罕见, 不允许自动替换
    'Hept': [],   # 七碳糖极罕见, 不允许自动替换
}

# 反向映射: 精确糖 → 其所属化学大类 (Reverse: precise sugar → class)
SUGAR_TO_CLASS = {}
for cls, sugars in SUGAR_CLASS_MAP.items():
    for s in sugars:
        SUGAR_TO_CLASS[s] = cls


def _isClassCompatible(genericLabel: str, preciseSugar: str) -> bool:
    """
    检查精确糖是否与泛指标签的碳类一致。
    Check if a precise sugar is carbon-class compatible with a generic label.

    例: Hex + D-Glc = True (同为六碳)
        Pen + L-Rha = False (五碳 vs 六碳脱氧, 碳类冲突)
        dHex + L-Rha = True (脱氧六碳一致)
    """
    allowed = SUGAR_CLASS_MAP.get(genericLabel, [])
    return preciseSugar in allowed


def runSafeMerge(df: pd.DataFrame) -> pd.DataFrame:
    """
    Pandas 防火墙 v3: 碳数一致性 + 多重候选保护锁。

    条件 A (防污染): 原序列已精确 -> 保持原样
    条件 B (精准反哺, 严格化学审核):
      B1: 碳类一致性 — NLP 糖必须属于泛指标签的同碳类
      B2: 多糖长度保护
      B3: 逐 token 精准对齐
      B4: 多重候选保护锁 — 同类多候选+多token时拒绝替换
    """
    print("\n" + "=" * 70)
    print("  Zero-Pollution Safe Merge v3 (Multi-Candidate Backoff)")
    print("  零污染安全合并 v3 (多重候选保护锁)")
    print("=" * 70)

    seqCol = "Sugar_Sequence"
    nlpCol = "NLP_Inferred_Sugars"
    consensusCol = "Consensus_Sugar_Sequence"

    df[consensusCol] = df[seqCol].copy()

    # 统计计数器 (Counters)
    protectedCount = 0       # 条件 A: 原序列已精确
    mergedCount = 0          # 条件 B: 成功合并
    noNlpCount = 0           # 无 NLP 数据
    carbonClashCount = 0     # 碳类冲突拦截
    multiCandidateBlock = 0  # 多候选拦截 (NEW in v3)
    noCompatibleCount = 0    # NLP 有糖但全部不兼容

    # 详细日志 (Detailed logs)
    acceptedLogs = []
    rejectedLogs = []

    for idx in df.index:
        origSeq = str(df.at[idx, seqCol]) if pd.notna(
            df.at[idx, seqCol]) else ""
        nlpSugars = str(df.at[idx, nlpCol]) if pd.notna(
            df.at[idx, nlpCol]) else ""

        if nlpSugars in ("", "nan", "None"):
            noNlpCount += 1
            continue

        # 条件 A: 原序列是否已精确? (Is original already precise?)
        hasGeneric = bool(GENERIC_PATTERN.search(origSeq))
        if not hasGeneric:
            protectedCount += 1
            continue

        # 解析 NLP 糖列表 (Parse NLP sugar list)
        nlpSugarList = [s.strip() for s in nlpSugars.split(";")
                        if s.strip()]
        if not nlpSugarList:
            noNlpCount += 1
            continue

        # 提取原序列中的所有泛指 token
        genericTokens = GENERIC_PATTERN.findall(origSeq)

        # ============================================================
        # 规则 B4: 多重候选保护锁 (Multiple-Candidate Backoff)
        # ============================================================
        # 对每种泛指标签, 统计:
        #   - 序列中该标签出现次数 (tokenCount)
        #   - NLP 中属于该碳类的候选糖去重数 (candidateCount)
        #
        # 唯一候选放行: candidateCount == 1 -> 大概率均多糖, 允许替换
        # 多重候选拦截: tokenCount >= 2 AND candidateCount >= 2
        #   -> 位置对应关系未知, 拒绝替换该标签类别

        from collections import Counter as _C
        tokenCounts = _C(genericTokens)  # e.g., {"Hex": 3, "Pen": 1}

        # 按碳类分组 NLP 候选 (Group NLP candidates by carbon class)
        classToNlpSugars = {}
        for nlpSugar in nlpSugarList:
            sugarClass = SUGAR_TO_CLASS.get(nlpSugar)
            if sugarClass:
                classToNlpSugars.setdefault(sugarClass, set()).add(nlpSugar)

        # 判断每种泛指标签是否被阻止
        blockedClasses = set()  # 被多候选保护锁阻止的碳类
        for genericLabel, count in tokenCounts.items():
            candidates = classToNlpSugars.get(genericLabel, set())
            if count >= 2 and len(candidates) >= 2:
                # 多 token + 多候选 = 拒绝!
                blockedClasses.add(genericLabel)
                multiCandidateBlock += 1

        # === 逐 token 替换 (带 B1+B4 双重检查) ===
        newSeq = origSeq
        anyReplaced = False
        anyRejected = False
        rejectionReason = ""

        def replaceTokenV3(m):
            nonlocal anyReplaced, anyRejected, rejectionReason
            nonlocal carbonClashCount

            token = m.group(0)  # e.g., "Hex", "Pen"
            labels = {"Hex", "Pen", "dHex", "HexA", "Non", "Oct", "Hept"}
            if token not in labels:
                return token

            # B4: 多重候选保护锁
            if token in blockedClasses:
                anyRejected = True
                candidates = classToNlpSugars.get(token, set())
                rejectionReason = (
                    f"{token}x{tokenCounts[token]} vs "
                    f"NLP[{','.join(sorted(candidates))}] "
                    f"(MULTI-CANDIDATE BACKOFF)")
                return token  # 保留原泛指标签

            # B1: 碳类一致性 — 寻找唯一兼容候选
            candidates = classToNlpSugars.get(token, set())
            if len(candidates) == 1:
                sugar = next(iter(candidates))
                anyReplaced = True
                return f"{sugar}(NLP)"

            # 回退: 在完整 NLP 列表中找第一个兼容的
            for nlpSugar in nlpSugarList:
                if _isClassCompatible(token, nlpSugar):
                    anyReplaced = True
                    return f"{nlpSugar}(NLP)"

            # 没找到兼容的 — 碳类冲突
            anyRejected = True
            rejectionReason = (
                f"{token} vs NLP[{','.join(nlpSugarList)}] "
                f"(class mismatch)")
            carbonClashCount += 1
            return token

        newSeq = GENERIC_PATTERN.sub(replaceTokenV3, origSeq)

        if anyReplaced and newSeq != origSeq:
            df.at[idx, consensusCol] = newSeq
            mergedCount += 1
            acceptedLogs.append({
                "idx": int(idx),
                "original": origSeq[:70],
                "nlp": nlpSugars[:50],
                "consensus": newSeq[:70],
            })
        elif anyRejected:
            rejectedLogs.append({
                "idx": int(idx),
                "original": origSeq[:70],
                "nlp": nlpSugars[:50],
                "reason": rejectionReason,
            })
            if not anyReplaced:
                noCompatibleCount += 1

    # 统计输出 (Stats output)
    hasGenericOrig = df[seqCol].str.contains(
        r'\bHex\b|\bPen\b|\bdHex\b|\bNon\b|\bOct\b|\bHept\b',
        na=False, regex=True).sum()
    hasGenericConsensus = df[consensusCol].str.contains(
        r'\bHex\b|\bPen\b|\bdHex\b|\bNon\b|\bOct\b|\bHept\b',
        na=False, regex=True).sum()

    print(f"\n  Firewall v3 Stats:")
    print(f"    Condition A (protected, precise): {protectedCount:,}")
    print(f"    Condition B (NLP merged):         {mergedCount:,}")
    print(f"    Carbon-class BLOCKED:             {carbonClashCount:,}")
    print(f"    Multi-candidate BLOCKED:          {multiCandidateBlock:,}")
    print(f"    No compatible sugar found:        {noCompatibleCount:,}")
    print(f"    No NLP data:                      {noNlpCount:,}")
    print(f"\n  Sugar_Sequence with generic: {hasGenericOrig:,}")
    print(f"  Consensus with generic:      {hasGenericConsensus:,}")
    print(f"  Reduction:                   -{hasGenericOrig - hasGenericConsensus:,}")

    # 多候选拦截日志 (Multi-candidate rejection log)
    multiBlocked = [r for r in rejectedLogs
                    if "MULTI-CANDIDATE" in r.get("reason", "")]
    if multiBlocked:
        print(f"\n  Multi-candidate BLOCKED samples (first 5):")
        for r in multiBlocked[:5]:
            print(f"    Row {r['idx']}: {r['original']}")
            print(f"      NLP: {r['nlp']}")
            print(f"      Reason: {r['reason']}")

    # 碳类拦截日志 (Carbon-class rejection log)
    carbonBlocked = [r for r in rejectedLogs
                     if "class mismatch" in r.get("reason", "")]
    if carbonBlocked:
        print(f"\n  Carbon-class REJECTED samples (first 5):")
        for r in carbonBlocked[:5]:
            print(f"    Row {r['idx']}: {r['original']}")
            print(f"      NLP: {r['nlp']}")
            print(f"      Reason: {r['reason']}")

    # 成功样本 (Accepted samples)
    if acceptedLogs:
        print(f"\n  Chemistry-APPROVED merge samples (first 5):")
        for a in acceptedLogs[:5]:
            print(f"    Row {a['idx']}: {a['original']}")
            print(f"      -> {a['consensus']}")

    return df


# =====================================================================
# Task 4: HTML Debug Report
# =====================================================================

def generateHtmlReport(
    df: pd.DataFrame,
    targetIdxs: list,
    outputPath: str,
) -> None:
    """生成 HTML 审核报告。"""
    print("\n" + "=" * 70)
    print("  HTML Debug Report Generation")
    print("=" * 70)

    subset = df.loc[targetIdxs].copy()

    # 统计 (Stats)
    hasNlpSugar = (subset.get("NLP_Inferred_Sugars", pd.Series(dtype=str))
                   .notna() & ~subset.get(
                       "NLP_Inferred_Sugars", pd.Series(dtype=str)
                   ).astype(str).isin(["", "nan", "None"])).sum()
    hasBioact = (subset.get("NLP_Bioactivity_Profile", pd.Series(dtype=str))
                 .notna() & ~subset.get(
                     "NLP_Bioactivity_Profile", pd.Series(dtype=str)
                 ).astype(str).isin(["", "nan", "None"])).sum()
    hasConsensus = subset.get(
        "Consensus_Sugar_Sequence", pd.Series(dtype=str)
    ).astype(str).str.contains("NLP", na=False).sum()

    htmlParts = [f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>GlycoNP NLP Literature Rescue Report</title>
<style>
body {{ font-family:'Segoe UI',Arial,sans-serif; margin:20px; background:#0d1117; color:#c9d1d9; }}
h1 {{ color:#58a6ff; border-bottom:2px solid #30363d; padding-bottom:10px; }}
table {{ border-collapse:collapse; width:100%; margin:10px 0; background:#161b22; }}
th {{ background:#21262d; color:#58a6ff; padding:10px; text-align:left; font-size:11px; border:1px solid #30363d; }}
td {{ padding:8px; border:1px solid #30363d; font-size:11px; max-width:280px; word-wrap:break-word; }}
tr:hover {{ background:#1c2333; }}
.fuzzy {{ background:#634e10; color:#f0c000; padding:2px 5px; border-radius:3px; font-weight:bold; }}
.precise {{ background:#0d3117; color:#3fb950; padding:2px 5px; border-radius:3px; }}
.nlp {{ background:#0c2d6b; color:#58a6ff; padding:2px 5px; border-radius:3px; }}
.bioact {{ background:#3d1f00; color:#f78166; padding:2px 5px; border-radius:3px; font-size:10px; }}
.target {{ background:#2d1530; color:#bc8cff; padding:2px 5px; border-radius:3px; font-size:10px; }}
.consensus {{ background:#1b3a28; color:#7ee787; padding:2px 5px; border-radius:3px; font-weight:bold; }}
.stat-row {{ display:flex; gap:15px; margin:15px 0; flex-wrap:wrap; }}
.stat-box {{ background:#161b22; border:1px solid #30363d; padding:15px 25px; border-radius:10px; text-align:center; }}
.stat-num {{ font-size:28px; font-weight:bold; color:#58a6ff; }}
.stat-label {{ font-size:11px; color:#8b949e; margin-top:5px; }}
.title-cell {{ color:#8b949e; font-size:10px; }}
</style></head><body>
<h1>GlycoNP NLP Literature Rescue - Debug Report</h1>
<div class="stat-row">
  <div class="stat-box"><div class="stat-num">{len(targetIdxs)}</div><div class="stat-label">Total Records</div></div>
  <div class="stat-box"><div class="stat-num">{hasNlpSugar}</div><div class="stat-label">Sugar Resolved</div></div>
  <div class="stat-box"><div class="stat-num">{hasBioact}</div><div class="stat-label">Bioactivity</div></div>
  <div class="stat-box"><div class="stat-num">{hasConsensus}</div><div class="stat-label">NLP Merged</div></div>
</div>
<table><tr>
  <th>#</th><th>Name</th><th>Sugar_Sequence</th><th>NLP_Sugars</th>
  <th>Consensus</th><th>Bioactivity</th><th>Targets</th>
  <th>PMC Title</th><th>DOI</th>
</tr>"""]

    for i, idx in enumerate(targetIdxs):
        row = df.loc[idx]
        name = html.escape(str(row.get("name", ""))[:50])
        origSeq = str(row.get("Sugar_Sequence", ""))
        nlpSug = str(row.get("NLP_Inferred_Sugars", ""))
        consensus = str(row.get("Consensus_Sugar_Sequence", ""))
        bioact = str(row.get("NLP_Bioactivity_Profile", ""))
        targets = str(row.get("NLP_Targets", ""))
        pmcTitle = str(row.get("PMC_Title", ""))[:100]
        doi = str(row.get("dois", ""))[:35]

        origHtml = GENERIC_PATTERN.sub(
            r'<span class="fuzzy">\1</span>', html.escape(origSeq))

        nlpHtml = "-"
        if nlpSug and nlpSug not in ("", "nan", "None"):
            nlpHtml = f'<span class="nlp">{html.escape(nlpSug)}</span>'

        consHtml = html.escape(consensus)
        if "(NLP)" in consensus:
            consHtml = f'<span class="consensus">{consHtml}</span>'
        elif not GENERIC_PATTERN.search(consensus):
            consHtml = f'<span class="precise">{consHtml}</span>'

        bioHtml = "-"
        if bioact and bioact not in ("", "nan", "None"):
            bioHtml = f'<span class="bioact">{html.escape(bioact[:80])}</span>'

        targHtml = "-"
        if targets and targets not in ("", "nan", "None"):
            targHtml = f'<span class="target">{html.escape(targets[:60])}</span>'

        titleHtml = f'<span class="title-cell">{html.escape(pmcTitle)}</span>'

        htmlParts.append(f"""<tr>
  <td>{i+1}</td><td>{name}</td><td>{origHtml}</td><td>{nlpHtml}</td>
  <td>{consHtml}</td><td>{bioHtml}</td><td>{targHtml}</td>
  <td>{titleHtml}</td><td style="font-size:9px">{html.escape(doi)}</td>
</tr>""")

    htmlParts.append("</table></body></html>")
    with open(outputPath, "w", encoding="utf-8") as f:
        f.write("\n".join(htmlParts))
    print(f"  HTML report: {outputPath}")


# =====================================================================
# Main
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        description="NLP Literature Mining Engine")
    parser.add_argument("--full", action="store_true",
                        help="Full mode: fetch all records with DOI")
    parser.add_argument("--limit", type=int, default=500)
    parser.add_argument("--skip-fetch", action="store_true",
                        help="Skip PMC fetch, use existing data")
    args = parser.parse_args()

    print("=" * 70)
    print("  GlycoNP Local NLP Literature Mining Engine")
    print("  本地 NLP 文献挖掘引擎 (无需 API)")
    print("=" * 70)

    inputCsv = INPUT_CSV
    if not os.path.exists(inputCsv):
        # 回退到其他输入文件
        for alt in ["GlycoNP_Fully_Enriched.csv",
                     "GlycoNP_Imputed_Merged.csv"]:
            altPath = os.path.join(REPORT_DIR, alt)
            if os.path.exists(altPath):
                print(f"  Fallback input: {alt}")
                inputCsv = altPath
                break

    df = pd.read_csv(inputCsv, dtype=str, low_memory=False,
                     encoding="utf-8-sig")
    total = len(df)
    print(f"  Loaded: {total:,} rows from {INPUT_CSV}")

    # 选择目标集 (Select targets)
    hasDoi = (df["dois"].notna() &
              ~df["dois"].astype(str).isin(["", "nan", "None"]))

    if args.full:
        # 全量: 所有有 DOI 的记录
        targetIdxs = df.index[hasDoi].tolist()
    else:
        # 试运行: 混合采样
        hasGeneric = df["Sugar_Sequence"].str.contains(
            r'\b(?:Hex|Pen|dHex|HexA|Non|Oct|Hept)\b',
            na=False, regex=True)
        fuzzyDoi = df.index[hasDoi & hasGeneric].tolist()
        preciseDoi = df.index[hasDoi & ~hasGeneric].tolist()

        import random
        random.seed(42)
        nFuzzy = min(int(args.limit * 0.8), len(fuzzyDoi))
        nPrecise = min(args.limit - nFuzzy, len(preciseDoi))
        targetIdxs = (random.sample(fuzzyDoi, nFuzzy) +
                      random.sample(preciseDoi, nPrecise))
        random.shuffle(targetIdxs)

    print(f"  Target records: {len(targetIdxs):,}")

    # ============ PMC Fetch ============
    if not args.skip_fetch:
        df = runPmcFetch(df, targetIdxs)
    else:
        print("\n  [SKIP] PMC fetch (--skip-fetch)")

    # ============ NLP Mining ============
    df = runNlpMining(df, targetIdxs)

    # ============ Safe Merge ============
    df = runSafeMerge(df)

    # ============ HTML Report ============
    htmlPath = os.path.join(REPORT_DIR,
                            "debug_sample_500_nlp_rescued.html")
    generateHtmlReport(df, targetIdxs, htmlPath)

    # ============ Pandas Merge Comparison Log ============
    print("\n" + "=" * 70)
    print("  Pandas Firewall Merge Comparison (Sample 10 rows)")
    print("=" * 70)
    mergedIdxs = [
        idx for idx in targetIdxs
        if "(NLP)" in str(df.at[idx, "Consensus_Sugar_Sequence"])
    ]
    for idx in mergedIdxs[:10]:
        print(f"\n  Row {idx}:")
        print(f"    Original:  {str(df.at[idx, 'Sugar_Sequence'])[:70]}")
        print(f"    NLP:       {str(df.at[idx, 'NLP_Inferred_Sugars'])[:70]}")
        print(f"    Consensus: {str(df.at[idx, 'Consensus_Sugar_Sequence'])[:70]}")
        bioact = str(df.at[idx, 'NLP_Bioactivity_Profile'])
        if bioact not in ("", "nan", "None"):
            print(f"    Bioact:    {bioact[:70]}")

    # ============ Save ============
    df.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")
    print(f"\n  Output: {OUTPUT_CSV}")
    print("=" * 70)


if __name__ == "__main__":
    main()
