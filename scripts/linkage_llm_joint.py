"""
Aglycone Linkage Probe + Local LLM Rescue Pipeline
====================================================

四大任务综合脚本:

Task 1: 苷元-糖链初始连接键探测 (Root Glycosidic Bond Probe)
  - O-苷 / C-苷 / N-苷 类型判定
  - α/β 构型推断 (CIP R/S → D/L config → alpha/beta)

Task 2: 本地 LLM 挖掘 (Local Ollama/llama-cpp Bridge)
  - 调用本地 Ollama 运行的模型
  - JSON Prompt 结构化提取

Task 3: v3 碳数防火墙 (Firewall Merge)
  - SUGAR_CLASS_MAP 碳类校验
  - Multi-Candidate Backoff

Task 4: 100-Sample Joint Trial
  - RDKit linkage probe + LLM + firewall 联合试跑

Usage:
  python scripts/linkage_llm_joint.py
"""
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

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from lib.bond_cleavage_engine import findAnomericCarbons, findGlycosidicBonds
from lib.glycan_topology import find_mapped_sugar_units

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
REPORT_DIR = os.path.join(BASE_DIR, "reports")

GENERIC_PATTERN = re.compile(r'\b(Hex|Pen|dHex|HexA|Non|Oct|Hept)\b')


# =====================================================================
# 碳类映射 (Carbon-Class Map) — 复用 v3 防火墙
# =====================================================================

SUGAR_CLASS_MAP = {
    'Hex':  ['D-Glc', 'D-Gal', 'D-Man', 'L-Glc', 'L-Gal', 'D-Tal',
             'D-All', 'D-Gul', 'L-Ido', 'D-Fru'],
    'Pen':  ['D-Xyl', 'L-Ara', 'D-Rib', 'D-Api', 'D-Lyx', 'L-Lyx'],
    'dHex': ['L-Rha', 'L-Fuc', 'D-Fuc', 'L-Ole', 'D-Cym', 'D-The',
             'D-Dig', 'D-Boi', 'D-Qui'],
    'HexN': ['D-GlcN', 'D-GalN', 'D-ManN'],
    'HexA': ['D-GlcA', 'D-GalA', 'D-ManA'],
    'Non':  ['Neu5Ac', 'Neu5Gc', 'KDO'],
    'Oct':  [],
    'Hept': [],
}
SUGAR_TO_CLASS = {}
for cls, sugars in SUGAR_CLASS_MAP.items():
    for s in sugars:
        SUGAR_TO_CLASS[s] = cls


# =====================================================================
# Task 1: 苷元-糖链初始连接键探测 (Root Glycosidic Bond Probe)
# =====================================================================

# D-糖的 C1 立体映射: CIP R→α, S→β (对 D-构型六碳糖)
# 这是简化规则, 不完美但覆盖大部分主流天然产物
# For D-sugars (most common): R at anomeric → α, S at anomeric → β

def probeRootLinkage(smiles: str) -> Dict:
    """
    探测苷元-糖链的初始连接键类型。
    Probe the root glycosidic bond connecting aglycone to first sugar.

    返回:
      {
        "linkage_type": "beta-O" / "alpha-O" / "R-O" / "?-C" / "N-glyc" / ...
        "bridge_element": "O" / "C" / "N" / "S"
        "cip_code": "R" / "S" / None
        "stereo_label": "alpha" / "beta" / "?" / None
        "bonds_found": int
      }
    """
    result = {
        "linkage_type": "",
        "bridge_element": "",
        "cip_code": "",
        "stereo_label": "",
        "bonds_found": 0,
    }

    if not smiles or str(smiles) in ("nan", "", "None", "NULL"):
        return result

    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return result

        # 分配立体化学 (Assign stereochemistry for CIP codes)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

        # 寻找糖环 (Find sugar rings)
        sugarUnits = find_mapped_sugar_units(mol)
        if not sugarUnits:
            return result

        # 寻找糖苷键 (Find glycosidic bonds)
        glycBonds = findGlycosidicBonds(mol, sugarUnits)
        result["bonds_found"] = len(glycBonds)

        # 仅关注 sugar_to_aglycon 键 (苷元连接)
        rootBonds = [b for b in glycBonds
                     if b["bond_type"] == "sugar_to_aglycon"]

        if not rootBonds:
            return result

        # 取第一条根键 (Take first root bond)
        rootBond = rootBonds[0]
        bridgeElem = rootBond["bridge_element"]
        anomericCIdx = rootBond["anomeric_carbon_idx"]

        result["bridge_element"] = bridgeElem

        # 桥原子类型 → 苷类型
        if bridgeElem == "O":
            glycType = "O"
        elif bridgeElem == "N":
            glycType = "N"
        elif bridgeElem == "S":
            glycType = "S"
        else:
            glycType = "C"

        # 提取异头碳 CIP 编码 (Extract CIP code of anomeric carbon)
        anomericAtom = mol.GetAtomWithIdx(anomericCIdx)
        cipCode = ""

        # 方法 1: 尝试 _CIPCode 属性
        if anomericAtom.HasProp("_CIPCode"):
            cipCode = anomericAtom.GetProp("_CIPCode")

        # 方法 2: 尝试 ChiralTag
        if not cipCode:
            chiralTag = anomericAtom.GetChiralTag()
            if chiralTag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
                cipCode = "CW"
            elif chiralTag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                cipCode = "CCW"

        result["cip_code"] = cipCode

        # α/β 推断 (Infer alpha/beta)
        # 简化规则: 对常见的 D-六碳糖 (大部分天然产物):
        #   CIP Code S at anomeric → beta (equatorial glycosidic bond)
        #   CIP Code R at anomeric → alpha (axial glycosidic bond)
        # 对 L-糖 (如 L-Rha): 反转
        stereoLabel = "?"
        if cipCode in ("R", "S"):
            # 默认 D-config (大部分天然产物糖)
            if cipCode == "R":
                stereoLabel = "alpha"
            elif cipCode == "S":
                stereoLabel = "beta"
        elif cipCode in ("CW", "CCW"):
            # 有手性标签但无 CIP → 保留方向信息
            stereoLabel = cipCode

        result["stereo_label"] = stereoLabel

        # 合成最终标签 (Compose final label)
        if stereoLabel in ("alpha", "beta"):
            result["linkage_type"] = f"{stereoLabel}-{glycType}"
        elif cipCode:
            result["linkage_type"] = f"{cipCode}-{glycType}"
        else:
            result["linkage_type"] = f"?-{glycType}"

    except Exception as e:
        result["linkage_type"] = f"ERROR:{str(e)[:30]}"

    return result


# =====================================================================
# Task 2: 本地 LLM 挖掘引擎 (Local Ollama Bridge)
# =====================================================================

# JSON 提取 Prompt — 复用之前为 Claude 设计的防污染模板
LLM_SYSTEM_PROMPT = """You are a glycochemistry expert. Extract sugar sequence and bioactivity information from scientific paper titles and abstracts about natural product glycosides.

Output ONLY valid JSON with these fields:
{
  "compound_name": "name of the glycosylated compound",
  "llm_sugar_sequence": ["D-Glc", "L-Rha", ...],
  "bioactivities": ["cytotoxic", "antibacterial", ...],
  "targets_and_diseases": ["A549", "S. aureus", ...]
}

Rules:
- Use standard sugar abbreviations: D-Glc, D-Gal, D-Man, L-Rha, L-Fuc, D-Xyl, L-Ara, D-GlcA, D-GalA, D-GlcNAc, Neu5Ac
- Only include sugars explicitly mentioned in the text
- If no sugar information found, return empty arrays
- Output ONLY the JSON, no other text"""


def callOllamaLocal(
    title: str,
    abstract: str,
    model: str = "llama3.2",
    temperature: float = 0.1,
    timeout: int = 30,
) -> Optional[Dict]:
    """
    调用本地 Ollama API 进行结构化提取。
    Call local Ollama API for structured extraction.
    """
    userPrompt = f"""Title: {title}

Abstract: {abstract[:1500] if abstract else 'N/A'}

Extract the sugar names and bioactivity from this paper."""

    payload = {
        "model": model,
        "messages": [
            {"role": "system", "content": LLM_SYSTEM_PROMPT},
            {"role": "user", "content": userPrompt},
        ],
        "stream": False,
        "options": {
            "temperature": temperature,
            "num_predict": 300,
        },
        "format": "json",
    }

    url = "http://localhost:11434/api/chat"

    try:
        data = json.dumps(payload).encode("utf-8")
        req = urllib.request.Request(
            url, data=data,
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            respData = json.loads(resp.read().decode("utf-8"))
            content = respData.get("message", {}).get("content", "")
            return json.loads(content)
    except urllib.error.URLError:
        return None  # Ollama not running
    except json.JSONDecodeError:
        return None
    except Exception:
        return None


def checkOllamaAvailable() -> bool:
    """检查 Ollama 本地服务是否可用。"""
    try:
        url = "http://localhost:11434/api/tags"
        req = urllib.request.Request(url, method="GET")
        with urllib.request.urlopen(req, timeout=3) as resp:
            data = json.loads(resp.read().decode("utf-8"))
            models = [m.get("name", "") for m in data.get("models", [])]
            if models:
                print(f"  Ollama available! Models: {', '.join(models)}")
                return True
    except Exception:
        pass
    return False


# =====================================================================
# Task 3: v3 碳数防火墙 (Firewall Merge for LLM Output)
# =====================================================================

def firewallMergeLlm(
    origSeq: str,
    llmSugars: List[str],
) -> Tuple[str, str]:
    """
    对 LLM 输出执行 v3 碳数一致性防火墙合并。
    Apply v3 chemistry-aware firewall to LLM-extracted sugars.

    Returns:
      (consensus_sequence, status) where status is:
        "MERGED", "BLOCKED:carbon_mismatch", "BLOCKED:multi_candidate",
        "NO_GENERIC", "NO_LLM"
    """
    if not llmSugars:
        return origSeq, "NO_LLM"

    if not GENERIC_PATTERN.search(origSeq):
        return origSeq, "NO_GENERIC"

    genericTokens = GENERIC_PATTERN.findall(origSeq)
    tokenCounts = Counter(genericTokens)

    # 按碳类分组 LLM 候选
    classToSugars = {}
    for sugar in llmSugars:
        sugarClass = SUGAR_TO_CLASS.get(sugar)
        if sugarClass:
            classToSugars.setdefault(sugarClass, set()).add(sugar)

    # 多重候选保护锁
    blockedClasses = set()
    for label, count in tokenCounts.items():
        candidates = classToSugars.get(label, set())
        if count >= 2 and len(candidates) >= 2:
            blockedClasses.add(label)

    anyReplaced = False
    anyBlocked = False
    blockReason = ""

    def replaceToken(m):
        nonlocal anyReplaced, anyBlocked, blockReason
        token = m.group(0)
        labels = {"Hex", "Pen", "dHex", "HexA", "Non", "Oct", "Hept"}
        if token not in labels:
            return token

        if token in blockedClasses:
            anyBlocked = True
            blockReason = "multi_candidate"
            return token

        candidates = classToSugars.get(token, set())
        if len(candidates) == 1:
            anyReplaced = True
            return f"{next(iter(candidates))}(LLM)"

        for sugar in llmSugars:
            allowed = SUGAR_CLASS_MAP.get(token, [])
            if sugar in allowed:
                anyReplaced = True
                return f"{sugar}(LLM)"

        anyBlocked = True
        blockReason = "carbon_mismatch"
        return token

    newSeq = GENERIC_PATTERN.sub(replaceToken, origSeq)

    if anyReplaced and newSeq != origSeq:
        return newSeq, "MERGED"
    elif anyBlocked:
        return origSeq, f"BLOCKED:{blockReason}"
    else:
        return origSeq, "NO_COMPATIBLE"


# =====================================================================
# Task 4: 100-Sample Joint Trial
# =====================================================================

def main():
    print("=" * 70)
    print("  Aglycone Linkage Probe + Local LLM Joint Trial")
    print("  苷元连接探测 + 本地 LLM 联合试跑")
    print("=" * 70)

    # 读取数据 (Load data)
    candidates = [
        "GlycoNP_Deep_Enriched.csv",
        "GlycoNP_NLP_Enriched.csv",
        "GlycoNP_LLM_Rescued.csv",
    ]
    inputCsv = None
    for c in candidates:
        p = os.path.join(REPORT_DIR, c)
        if os.path.exists(p):
            inputCsv = p
            break

    df = pd.read_csv(inputCsv, dtype=str, low_memory=False,
                     encoding="utf-8-sig")
    total = len(df)
    print(f"  Loaded: {total:,} rows from {os.path.basename(inputCsv)}")

    # 选择100个复杂样本 (Select 100 complex samples)
    # 要求: 有 SMILES + 有模糊标签 + 有 PMC Title
    import random
    random.seed(42)

    seqCol = "Consensus_Sugar_Sequence"
    if seqCol not in df.columns:
        seqCol = "Sugar_Sequence"

    hasFuzzy = df[seqCol].str.contains(
        r'\b(?:Hex|Pen|dHex|HexA|Non)\b', na=False, regex=True)
    hasSmiles = (df["canonical_smiles"].notna() &
                 (df["canonical_smiles"].str.len() > 10))
    hasPmc = (df.get("PMC_Title", pd.Series(dtype=str)).notna() &
              ~df.get("PMC_Title", pd.Series(dtype=str)).astype(str).isin(
                  ["", "nan", "None"]))

    complexMask = hasFuzzy & hasSmiles
    complexIdxs = df.index[complexMask].tolist()
    pmcComplexIdxs = df.index[complexMask & hasPmc].tolist()

    # 优先选带 PMC 的, 补足 100 个
    if len(pmcComplexIdxs) >= 100:
        sampleIdxs = random.sample(pmcComplexIdxs, 100)
    else:
        sampleIdxs = pmcComplexIdxs.copy()
        remaining = [i for i in complexIdxs if i not in set(sampleIdxs)]
        nMore = min(100 - len(sampleIdxs), len(remaining))
        sampleIdxs += random.sample(remaining, nMore)

    print(f"  Trial samples: {len(sampleIdxs)} "
          f"(with PMC: {len(pmcComplexIdxs)}, total complex: {len(complexIdxs):,})")

    # ============================
    # Task 1: Root Linkage Probe
    # ============================
    print("\n" + "=" * 70)
    print("  Task 1: Root Glycosidic Bond Probe")
    print("  苷元-糖链初始连接键探测")
    print("=" * 70)

    linkageResults = []
    for i, idx in enumerate(sampleIdxs):
        smi = str(df.at[idx, "canonical_smiles"])
        result = probeRootLinkage(smi)
        linkageResults.append(result)
        if (i + 1) % 25 == 0:
            print(f"    [{i+1}/100] processed")

    # 写入新列 (Write to new column)
    for col in ["Aglycone_Linkage_Type", "Root_Bridge_Element",
                 "Root_CIP_Code"]:
        if col not in df.columns:
            df[col] = ""

    for i, idx in enumerate(sampleIdxs):
        r = linkageResults[i]
        df.at[idx, "Aglycone_Linkage_Type"] = r["linkage_type"]
        df.at[idx, "Root_Bridge_Element"] = r["bridge_element"]
        df.at[idx, "Root_CIP_Code"] = r["cip_code"]

    # 统计 (Stats)
    linkTypes = Counter(r["linkage_type"] for r in linkageResults
                        if r["linkage_type"])
    bridgeTypes = Counter(r["bridge_element"] for r in linkageResults
                          if r["bridge_element"])

    print(f"\n  Linkage Type Distribution (100 samples):")
    for lt, cnt in linkTypes.most_common():
        print(f"    {lt:15s}: {cnt}")
    print(f"\n  Bridge Element Distribution:")
    for be, cnt in bridgeTypes.most_common():
        label = {"O": "O-glycoside", "N": "N-glycoside",
                 "C": "C-glycoside", "S": "S-glycoside"}.get(be, be)
        print(f"    {label:15s}: {cnt}")

    # 审查输出 A: 5 个连接探测样本
    print(f"\n  === Audit A: 5 Linkage Probe Samples ===")
    shown = 0
    for i, idx in enumerate(sampleIdxs):
        r = linkageResults[i]
        if r["linkage_type"] and r["linkage_type"] != "?-O" and shown < 5:
            scaffold = str(df.at[idx, "Murcko_Scaffold"])[:40]
            seq = str(df.at[idx, seqCol])[:50]
            npClass = str(df.at[idx, "Detailed_NP_Class"])[:30] if \
                "Detailed_NP_Class" in df.columns else ""
            print(f"\n    Sample {shown+1} (Row {idx}):")
            print(f"      NP Class:         {npClass}")
            print(f"      Linkage Type:     {r['linkage_type']}")
            print(f"      CIP Code:         {r['cip_code']}")
            print(f"      Bridge:           {r['bridge_element']}")
            print(f"      Bonds Found:      {r['bonds_found']}")
            print(f"      Sugar Sequence:   {seq}")
            shown += 1

    # ============================
    # Task 2: Local LLM
    # ============================
    print("\n" + "=" * 70)
    print("  Task 2: Local LLM Rescue Engine")
    print("  本地 LLM 挖掘引擎")
    print("=" * 70)

    ollamaReady = checkOllamaAvailable()

    llmResults = {}
    if ollamaReady:
        print("  Starting Ollama inference...")
        for i, idx in enumerate(sampleIdxs):
            title = str(df.at[idx, "PMC_Title"]) if pd.notna(
                df.at[idx, "PMC_Title"]) else ""
            abstract = str(df.at[idx, "PMC_Abstract"]) if pd.notna(
                df.at[idx, "PMC_Abstract"]) else ""
            if title in ("nan", "None", ""):
                continue

            llmOut = callOllamaLocal(title, abstract)
            if llmOut:
                llmResults[idx] = llmOut

            if (i + 1) % 10 == 0:
                print(f"    [{i+1}/100] LLM calls, "
                      f"success={len(llmResults)}")
            time.sleep(0.3)
    else:
        print("  [WARN] Ollama not available at localhost:11434")
        print("  Falling back to NLP regex engine...")
        print("")
        print("  To enable local LLM, install Ollama:")
        print("    1. Download from https://ollama.com/download")
        print("    2. Run:  ollama pull llama3.2")
        print("    3. Re-run this script")
        print("")
        print("  Using NLP regex fallback for this trial...")

        # 回退到 NLP 正则引擎 (Fallback to regex)
        # 导入 NLP mining 函数
        sys.path.insert(0, os.path.join(BASE_DIR, "scripts"))
        from nlp_literature_mining import mineText

        for i, idx in enumerate(sampleIdxs):
            title = str(df.at[idx, "PMC_Title"]) if pd.notna(
                df.at[idx, "PMC_Title"]) else ""
            abstract = str(df.at[idx, "PMC_Abstract"]) if pd.notna(
                df.at[idx, "PMC_Abstract"]) else ""
            name = str(df.at[idx, "name"]) if pd.notna(
                df.at[idx, "name"]) else ""
            if title in ("nan", "None"):
                title = ""
            if abstract in ("nan", "None"):
                abstract = ""
            if name in ("nan", "None"):
                name = ""

            if not title and not name:
                continue

            result = mineText(title, abstract + " " + name)
            if result["sugars"]:
                llmResults[idx] = {
                    "compound_name": name[:50],
                    "llm_sugar_sequence": result["sugars"],
                    "bioactivities": result["bioactivities"],
                    "targets_and_diseases": result["targets"],
                }

    print(f"\n  LLM/NLP extraction results: {len(llmResults)}")

    # ============================
    # Task 3: v3 Firewall Merge
    # ============================
    print("\n" + "=" * 70)
    print("  Task 3: v3 Chemistry Firewall Merge")
    print("  v3 碳数防火墙合并")
    print("=" * 70)

    fwStats = Counter()
    mergedSamples = []
    blockedSamples = []

    for idx, llmData in llmResults.items():
        origSeq = str(df.at[idx, seqCol])
        llmSugars = llmData.get("llm_sugar_sequence", [])

        consensus, status = firewallMergeLlm(origSeq, llmSugars)
        fwStats[status] += 1

        if status == "MERGED":
            mergedSamples.append({
                "idx": int(idx),
                "original": origSeq[:60],
                "llm_sugars": llmSugars,
                "consensus": consensus[:60],
                "bioact": llmData.get("bioactivities", []),
            })
        elif "BLOCKED" in status:
            blockedSamples.append({
                "idx": int(idx),
                "original": origSeq[:60],
                "llm_sugars": llmSugars,
                "reason": status,
            })

    print(f"\n  Firewall v3 Stats:")
    for status, count in fwStats.most_common():
        print(f"    {status:30s}: {count}")

    # 审查输出 B: 3 成功 + 1 拦截
    print(f"\n  === Audit B: LLM + Firewall Results ===")
    print(f"\n  [APPROVED] Successfully merged samples:")
    for s in mergedSamples[:3]:
        print(f"\n    Row {s['idx']}:")
        print(f"      Original:    {s['original']}")
        print(f"      LLM Sugars:  {s['llm_sugars']}")
        print(f"      Consensus:   {s['consensus']}")
        if s['bioact']:
            print(f"      Bioactivity: {', '.join(s['bioact'][:3])}")

    if blockedSamples:
        print(f"\n  [BLOCKED] Firewall-rejected samples:")
        for s in blockedSamples[:2]:
            print(f"\n    Row {s['idx']}:")
            print(f"      Original:    {s['original']}")
            print(f"      LLM Sugars:  {s['llm_sugars']}")
            print(f"      Reason:      {s['reason']}")

    # ============================
    # Summary
    # ============================
    print("\n" + "=" * 70)
    print("  Joint Trial Summary")
    print("=" * 70)
    print(f"  Task 1 (Linkage Probe):    {sum(1 for r in linkageResults if r['linkage_type'])}/100 detected")
    print(f"  Task 2 (LLM/NLP Extract):  {len(llmResults)}/100 extracted")
    print(f"  Task 3 (Firewall Merge):   {fwStats.get('MERGED', 0)} merged, "
          f"{sum(v for k,v in fwStats.items() if 'BLOCKED' in k)} blocked")

    # Save
    outputPath = os.path.join(REPORT_DIR, "GlycoNP_Linkage_Enriched.csv")
    df.to_csv(outputPath, index=False, encoding="utf-8-sig")
    print(f"\n  Output: {outputPath}")
    print("=" * 70)


if __name__ == "__main__":
    main()
