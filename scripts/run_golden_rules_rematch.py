"""
贝叶斯先验 + 全量竞争 v6.0 终极重刷 + NLP 否决
Bayesian Prior v6.0 + Full Competition + NLP Veto

Tier 1-3 集成在 identify_monosaccharide_v2 (全量竞争 + 贝叶斯裁决)
NLP 否决权在此脚本中后处理
"""
import sys, os, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))

import pandas as pd
from rdkit import Chem
from tqdm import tqdm
from collections import Counter

from lib.monosaccharide_identifier import analyze_glycan

# ---- 配置 ----
INPUT_CSV = "reports/GlycoNP_Deep_Enriched_v7.csv"
OUTPUT_CSV = "reports/GlycoNP_Deep_Enriched_v8.csv"

# 碳数查表 (Carbon count lookup for NLP rescue)
SUGAR_CARBON_COUNT = {
    "D-Glc": 6, "L-Glc": 6, "D-Gal": 6, "L-Gal": 6, "D-Man": 6, "L-Man": 6,
    "D-All": 6, "D-Alt": 6, "D-Gul": 6, "D-Ido": 6, "D-Tal": 6,
    "D-GlcA": 6, "D-GalA": 6, "D-ManA": 6, "D-GulA": 6, "D-IdoA": 6,
    "D-GlcNAc": 8, "D-GalNAc": 8, "D-ManNAc": 8,
    "D-GlcN": 6, "D-GalN": 6, "D-ManN": 6,
    "L-Rha": 6, "D-Rha": 6, "L-Fuc": 6, "D-Fuc": 6, "D-Qui": 6,
    "D-Xyl": 5, "L-Ara": 5, "D-Ara": 5, "D-Rib": 5, "D-Lyx": 5, "L-Lyx": 5,
    "D-Api": 5, "D-Fru": 6, "Neu5Ac": 9, "KDO": 8,
}

print("=" * 70)
print("贝叶斯先验 v6.0 + 全量竞争 + NLP 否决 / Bayesian Prior v6.0")
print("=" * 70)

# 加载数据
df = pd.read_csv(INPUT_CSV, low_memory=False)
print(f"Loaded {len(df)} rows from {INPUT_CSV}")

glycanCol = "Glycan_SMILES"
nlpCol = "NLP_Inferred_Sugars"
hasGlycan = df[glycanCol].notna() & (df[glycanCol] != "") & (df[glycanCol] != "nan")
print(f"Rows with Glycan_SMILES: {hasGlycan.sum()}")
nlpAvail = df[nlpCol].notna().sum() if nlpCol in df.columns else 0
print(f"Rows with NLP_Inferred_Sugars: {nlpAvail}")

# ================================================================
# Step 1: 运行 v6.0 全量竞争引擎
# ================================================================
print("\nStep 1: Running v6.0 engine (Bayesian + Full Competition)...")
newSeqs = []
newMods = []
for idx, row in tqdm(df.iterrows(), total=len(df), desc="v6.0 Match"):
    smiles = row.get(glycanCol, "")
    if pd.isna(smiles) or str(smiles).strip() in ("", "nan", "NULL"):
        newSeqs.append("")
        newMods.append("")
    else:
        try:
            seq, mod = analyze_glycan(str(smiles))
        except Exception:
            seq, mod = "", ""
        newSeqs.append(seq)
        newMods.append(mod)

df["Sugar_Sequence"] = newSeqs
df["Sugar_Mods_Detail"] = newMods

# ================================================================
# Step 2: 法则 4 — NLP 文献交叉验证
# ================================================================
print("\nStep 2: NLP Literature Rescue (Rule 4)...")

nlpRescueCount = 0
nlpRescueDetails = Counter()

for idx in tqdm(range(len(df)), desc="NLP Rescue"):
    seq = str(df.at[idx, "Sugar_Sequence"]).strip()
    if not seq or seq in ("", "nan"):
        continue

    nlpRaw = df.at[idx, nlpCol] if nlpCol in df.columns else None
    if pd.isna(nlpRaw) or not str(nlpRaw).strip():
        continue

    # 解析 NLP 候选糖列表
    nlpSugars = [s.strip() for s in str(nlpRaw).split(";") if s.strip()]
    if not nlpSugars:
        continue

    # 检查序列中是否有需要救援的标签
    newSeq = seq
    needSave = False
    nlpSet = set(nlpSugars)

    # 策略 V: NLP 一票否决权 (NLP Veto)
    # 罕见糖 (D-Tal, D-All 等) 如果 NLP 明确指出常见替代 → 强制覆写
    RARE_VETO_TARGETS = {
        "D-Tal": ["D-Gal", "D-Man"],  # D-Tal 是 D-Gal 的 C2 差向异构体
        "D-All": ["D-Glc"],
        "D-Alt": ["D-Man"],
        "D-Gul": ["D-Gal"],
        "D-Ido": ["D-Gal"],
        "L-Tal": ["L-Gal"],
    }
    for rareSugar, commonAlts in RARE_VETO_TARGETS.items():
        if rareSugar in newSeq:
            for alt in commonAlts:
                if alt in nlpSet:
                    newSeq = newSeq.replace(rareSugar, f"{alt}(Veto)")
                    needSave = True
                    nlpRescueDetails[f"{rareSugar}→{alt}(Veto)"] += seq.count(rareSugar)
                    break

    # 策略 A: Hex → NLP 唯一 6C 候选
    if "Hex" in newSeq:
        nlp6C = [s for s in nlpSugars if SUGAR_CARBON_COUNT.get(s, 0) == 6]
        if len(set(nlp6C)) <= 2 and len(nlp6C) > 0:
            bestNlpHex = Counter(nlp6C).most_common(1)[0][0]
            newSeq = newSeq.replace("Hex", f"{bestNlpHex}(NLP)")
            needSave = True
            nlpRescueDetails[f"Hex→{bestNlpHex}(NLP)"] += seq.count("Hex")

    # 策略 B: D-Rib → NLP 明确的其他 5C 糖
    if "D-Rib" in newSeq:
        nlp5C = [s for s in nlpSugars if SUGAR_CARBON_COUNT.get(s, 0) == 5 and s != "D-Rib"]
        if nlp5C:
            bestNlp5C = Counter(nlp5C).most_common(1)[0][0]
            newSeq = newSeq.replace("D-Rib", f"{bestNlp5C}(NLP)")
            needSave = True
            nlpRescueDetails[f"D-Rib→{bestNlp5C}(NLP)"] += seq.count("D-Rib")

    # 策略 C: Pen → NLP 具体 5C 糖
    if "Pen" in newSeq:
        nlp5C = [s for s in nlpSugars if SUGAR_CARBON_COUNT.get(s, 0) == 5]
        if len(set(nlp5C)) <= 2 and len(nlp5C) > 0:
            bestNlpPen = Counter(nlp5C).most_common(1)[0][0]
            newSeq = newSeq.replace("Pen", f"{bestNlpPen}(NLP)")
            needSave = True
            nlpRescueDetails[f"Pen→{bestNlpPen}(NLP)"] += seq.count("Pen")

    if needSave:
        df.at[idx, "Sugar_Sequence"] = newSeq
        nlpRescueCount += 1

print(f"\nNLP Rescue applied to {nlpRescueCount} rows")
print("NLP Rescue details:")
for detail, count in nlpRescueDetails.most_common(20):
    print(f"  {detail}: {count}")

# 保存
df.to_csv(OUTPUT_CSV, index=False)
print(f"\nSaved to {OUTPUT_CSV}")

# ================================================================
# Step 3: 最终统计
# ================================================================
print("\n" + "=" * 70)
print("Final Top 20 Sugar Distribution (v4.0 Golden Rules + NLP)")
print("=" * 70)

allNames = []
for seq in df["Sugar_Sequence"]:
    seq = str(seq).strip()
    if not seq or seq in ("", "nan", "Invalid", "Error"):
        continue
    for chain in seq.split(" ; "):
        names = re.findall(
            r'([DL]-[A-Za-z0-9]+(?:\((?:NLP|Rescued|Veto)\))?|Hex[N]?|Pen[N]?|dHex|Oct|Hept|Neu5Ac|Neu5Gc|KDO|Kdn)',
            chain
        )
        allNames.extend(names)

counter = Counter(allNames)
total = len(allNames)
print(f"Total sugar instances: {total}")
print(f"Unique names: {len(counter)}")
for rank, (name, count) in enumerate(counter.most_common(20), 1):
    pct = count / total * 100 if total > 0 else 0
    print(f"  {rank:>2}. {name:<30} {count:>8}   ({pct:.1f}%)")

dRib = counter.get("D-Rib", 0) + counter.get("D-Ribf", 0)
hexG = counter.get("Hex", 0)
dTal = counter.get("D-Tal", 0)
print(f"\n=== KEY METRICS ===")
print(f"D-Rib (total):    {dRib}")
print(f"Hex (generic):    {hexG}")
print(f"D-Tal:            {dTal}")
print(f"D-Rib + Hex:      {dRib + hexG}")
