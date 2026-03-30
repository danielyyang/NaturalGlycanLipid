"""
Local Ollama LLM Sugar Sequence Cleanup
========================================

本地 Ollama 大模型糖序列终极清扫 (工业级容灾):

Task 1: 本地 Ollama API 引擎 (llama3.2, temperature=0.0)
Task 2: 断点续传 (Checkpoint every 500 rows)
Task 3: v3 碳数防火墙 + 多重候选保护锁

Usage:
  python scripts/local_llm_cleanup.py
  python scripts/local_llm_cleanup.py --dry-run
  python scripts/local_llm_cleanup.py --limit 50
  python scripts/local_llm_cleanup.py --model llama3.2:latest
"""
import json
import os
import re
import sys
import time
from collections import Counter
from typing import Dict, List, Optional, Tuple

import pandas as pd
import requests

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
REPORT_DIR = os.path.join(BASE_DIR, "reports")

# =====================================================================
# Configuration (可配置项)
# =====================================================================

OLLAMA_URL = os.environ.get("OLLAMA_URL", "http://localhost:11434")
DEFAULT_MODEL = "llama3.2"

# 断点续传文件 (Checkpoint file)
CHECKPOINT_PATH = os.path.join(REPORT_DIR, "GlycoNP_LocalLLM_Checkpoint.csv")
CHECKPOINT_INTERVAL = 500

# API 调用间隔 (本地模型无需限流, 但留缓冲)
API_DELAY_SECONDS = 0.1
MAX_RETRIES = 3

GENERIC_PATTERN = re.compile(r'\b(Hex|Pen|dHex|HexA|Non|Oct|Hept)\b')

# =====================================================================
# 碳类映射 (Carbon-Class Map) — v3 防火墙
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
# Task 1: 本地 Ollama API 客户端
# =====================================================================

LLM_SYSTEM_PROMPT = """You are a glycochemistry expert. Extract sugar names and bioactivity from a scientific paper about glycosylated natural products.

Output ONLY a valid JSON object, no explanation, no markdown:
{"compound_name":"name","llm_sugar_sequence":["D-Glc","L-Rha"],"bioactivities":["cytotoxic"],"targets_and_diseases":["A549"]}

Standard sugar names to use:
Hexoses: D-Glc, D-Gal, D-Man, L-Glc, L-Gal, D-Tal, D-All, D-Fru
Pentoses: D-Xyl, L-Ara, D-Rib, D-Api, D-Lyx
Deoxysugars: L-Rha, L-Fuc, D-Fuc, L-Ole, D-Cym, D-The, D-Dig
Amino sugars: D-GlcN, D-GalN, D-GlcNAc, D-GalNAc
Uronic acids: D-GlcA, D-GalA, D-ManA
Special: Neu5Ac, KDO

Rules: Only include sugars EXPLICITLY mentioned. Return empty arrays if not found. Output ONLY JSON."""


def sanitizeJson(rawText: str) -> Optional[Dict]:
    """
    清理本地模型输出的 JSON (本地小模型比 Opus 更容易输出多余文字)。
    Sanitize JSON from local model output — small models often add
    explanatory text around the JSON.
    """
    if not rawText:
        return None

    text = rawText.strip()

    # 1. 剥离 markdown 代码块 (Strip ```json ... ```)
    text = re.sub(r'^```(?:json)?\s*\n?', '', text, flags=re.MULTILINE)
    text = re.sub(r'\n?```\s*$', '', text, flags=re.MULTILINE)
    text = text.strip()

    # 2. 如果有前后文字, 提取第一个 JSON 对象
    jsonMatch = re.search(r'\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}', text,
                          re.DOTALL)
    if jsonMatch:
        text = jsonMatch.group(0)

    # 3. 尝试解析
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass

    # 4. 替换单引号
    try:
        return json.loads(text.replace("'", '"'))
    except json.JSONDecodeError:
        pass

    # 5. 修复尾部逗号 (trailing comma before })
    try:
        fixed = re.sub(r',\s*}', '}', text)
        fixed = re.sub(r',\s*]', ']', fixed)
        return json.loads(fixed)
    except json.JSONDecodeError:
        pass

    return None


def callOllama(
    title: str,
    abstract: str,
    compoundName: str = "",
    model: str = DEFAULT_MODEL,
) -> Optional[Dict]:
    """
    调用本地 Ollama REST API。
    Call local Ollama /api/chat endpoint with retry.
    """
    userPrompt = (
        f"Paper Title: {title}\n"
        f"Abstract: {abstract[:1500] if abstract else 'N/A'}\n"
        f"Compound: {compoundName if compoundName else 'Unknown'}\n"
        f"Extract sugars and bioactivity as JSON."
    )

    payload = {
        "model": model,
        "messages": [
            {"role": "system", "content": LLM_SYSTEM_PROMPT},
            {"role": "user", "content": userPrompt},
        ],
        "stream": False,
        "format": "json",
        "options": {
            "temperature": 0.0,
            "num_predict": 300,
        },
    }

    url = f"{OLLAMA_URL}/api/chat"

    for attempt in range(MAX_RETRIES):
        try:
            resp = requests.post(url, json=payload, timeout=90)

            if resp.status_code == 200:
                data = resp.json()
                content = data.get("message", {}).get("content", "")
                return sanitizeJson(content)

            elif resp.status_code in (429, 500, 502, 503):
                wait = (2 ** attempt) + 1
                print(f"      [RETRY {attempt+1}/{MAX_RETRIES}] "
                      f"HTTP {resp.status_code}, wait {wait}s")
                time.sleep(wait)
                continue

            else:
                print(f"      [ERROR] HTTP {resp.status_code}: "
                      f"{resp.text[:150]}")
                return None

        except requests.exceptions.Timeout:
            wait = (2 ** attempt) + 2
            print(f"      [TIMEOUT] attempt {attempt+1}, wait {wait}s")
            time.sleep(wait)

        except requests.exceptions.ConnectionError:
            if attempt == 0:
                print(f"      [CONN ERROR] Ollama not reachable at {url}")
            return None

        except Exception as e:
            print(f"      [ERROR] {type(e).__name__}: {str(e)[:80]}")
            return None

    return None


def checkOllamaReady(model: str) -> bool:
    """检查 Ollama 是否运行且模型已加载。"""
    try:
        resp = requests.get(f"{OLLAMA_URL}/api/tags", timeout=5)
        if resp.status_code == 200:
            data = resp.json()
            models = [m.get("name", "") for m in data.get("models", [])]
            available = [m for m in models if model in m]
            if available:
                print(f"  Ollama OK! Model '{available[0]}' ready.")
                return True
            else:
                print(f"  [WARN] Model '{model}' not found.")
                print(f"  Available: {', '.join(models) if models else 'none'}")
                print(f"  Run: ollama pull {model}")
                return False
    except Exception:
        pass
    print(f"  [ERROR] Ollama not running at {OLLAMA_URL}")
    print(f"  Start Ollama first, then run: ollama pull {model}")
    return False


# =====================================================================
# Task 3: v3 碳数防火墙
# =====================================================================

def firewallMerge(origSeq: str, llmSugars: List[str]) -> Tuple[str, str]:
    """v3 防火墙: 碳类一致性 + 多重候选保护锁。"""
    if not llmSugars:
        return origSeq, "NO_LLM"
    if not GENERIC_PATTERN.search(str(origSeq)):
        return origSeq, "NO_GENERIC"

    genericTokens = GENERIC_PATTERN.findall(str(origSeq))
    tokenCounts = Counter(genericTokens)

    classToSugars = {}
    for sugar in llmSugars:
        sc = SUGAR_TO_CLASS.get(sugar)
        if sc:
            classToSugars.setdefault(sc, set()).add(sugar)

    blockedClasses = set()
    for label, count in tokenCounts.items():
        cands = classToSugars.get(label, set())
        if count >= 2 and len(cands) >= 2:
            blockedClasses.add(label)

    anyReplaced = False
    anyBlocked = False
    blockReason = ""

    def replaceToken(m):
        nonlocal anyReplaced, anyBlocked, blockReason
        token = m.group(0)
        if token not in {"Hex", "Pen", "dHex", "HexA", "Non", "Oct", "Hept"}:
            return token
        if token in blockedClasses:
            anyBlocked = True
            blockReason = "multi_candidate"
            return token
        cands = classToSugars.get(token, set())
        if len(cands) == 1:
            anyReplaced = True
            return f"{next(iter(cands))}(LLM)"
        for sugar in llmSugars:
            if sugar in SUGAR_CLASS_MAP.get(token, []):
                anyReplaced = True
                return f"{sugar}(LLM)"
        anyBlocked = True
        blockReason = "carbon_mismatch"
        return token

    newSeq = GENERIC_PATTERN.sub(replaceToken, str(origSeq))

    if anyReplaced and newSeq != origSeq:
        return newSeq, "MERGED"
    elif anyBlocked:
        return origSeq, f"BLOCKED:{blockReason}"
    return origSeq, "NO_COMPATIBLE"


# =====================================================================
# Task 2: 主循环 + 断点续传
# =====================================================================

def loadCheckpoint() -> Tuple[pd.DataFrame, set]:
    """加载断点文件。"""
    if os.path.exists(CHECKPOINT_PATH):
        ckpt = pd.read_csv(CHECKPOINT_PATH, dtype=str, encoding="utf-8-sig")
        processed = set(ckpt["row_idx"].astype(int)) if "row_idx" in ckpt.columns else set()
        print(f"  [RESUME] {len(processed):,} rows from checkpoint")
        return ckpt, processed
    return pd.DataFrame(), set()


def saveCheckpoint(results: List[Dict]) -> None:
    """保存断点。"""
    if results:
        pd.DataFrame(results).to_csv(
            CHECKPOINT_PATH, index=False, encoding="utf-8-sig")


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Local Ollama LLM Cleanup")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--model", type=str, default=DEFAULT_MODEL)
    args = parser.parse_args()

    model = args.model

    print("=" * 70)
    print("  Local Ollama LLM Sugar Sequence Cleanup")
    print("  本地 Ollama 糖序列终极清扫")
    print("=" * 70)
    print(f"  Ollama URL:   {OLLAMA_URL}")
    print(f"  Model:        {model}")
    print(f"  Temperature:  0.0 (deterministic)")
    print(f"  JSON format:  forced via 'format: json'")
    print(f"  Checkpoint:   {os.path.basename(CHECKPOINT_PATH)}")
    print(f"  Interval:     every {CHECKPOINT_INTERVAL} rows")

    if args.dry_run:
        print("\n  Checking Ollama...")
        checkOllamaReady(model)
        print("\n  [DRY RUN] Config verified. Exiting.")
        return

    # 检查 Ollama (Check Ollama)
    if not checkOllamaReady(model):
        print("\n  Aborting. Please start Ollama and pull model first.")
        return

    # 读取数据 (Load data)
    inputCsv = os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched_v2.csv")
    if not os.path.exists(inputCsv):
        inputCsv = os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched.csv")

    df = pd.read_csv(inputCsv, dtype=str, low_memory=False,
                     encoding="utf-8-sig")
    print(f"\n  Loaded: {len(df):,} rows")

    seqCol = "Consensus_Sugar_Sequence"
    if seqCol not in df.columns:
        seqCol = "Sugar_Sequence"

    # 筛选目标 (Select targets: fuzzy labels + has text)
    hasFuzzy = df[seqCol].str.contains(
        r'\b(?:Hex|Pen|dHex|HexA|Non|Oct|Hept)\b', na=False, regex=True)
    hasText = pd.Series(False, index=df.index)
    for col in ["PMC_Title", "name"]:
        if col in df.columns:
            hasText = hasText | (
                df[col].notna() &
                ~df[col].astype(str).isin(["", "nan", "None"]))

    targetIdxs = df.index[hasFuzzy & hasText].tolist()
    print(f"  Target records: {len(targetIdxs):,}")

    if args.limit > 0:
        targetIdxs = targetIdxs[:args.limit]
        print(f"  Limited to: {len(targetIdxs)}")

    # 断点续传 (Resume from checkpoint)
    ckptDf, processedRows = loadCheckpoint()
    allResults = ckptDf.to_dict("records") if not ckptDf.empty else []
    remainingIdxs = [i for i in targetIdxs if i not in processedRows]
    print(f"  Remaining: {len(remainingIdxs):,}")

    if not remainingIdxs:
        print("  All done! Proceeding to final merge...")
    else:
        # 预热测试 (Warmup test)
        print(f"\n  Warmup inference...")
        t0w = time.time()
        testResult = callOllama("Rutin glycoside isolation from buckwheat",
                                "A flavonoid di-glycoside with glucose and rhamnose.",
                                model=model)
        warmupTime = time.time() - t0w
        if testResult:
            print(f"  Warmup OK ({warmupTime:.1f}s): {testResult}")
        else:
            print(f"  [ERROR] Warmup failed! Check model.")
            return

        # 估算时间
        perCall = max(warmupTime, 1.0) + API_DELAY_SECONDS
        etaMin = len(remainingIdxs) * perCall / 60
        print(f"\n  Estimated time: ~{etaMin:.0f} min "
              f"({perCall:.1f}s/call x {len(remainingIdxs):,} records)")

    # 主循环 (Main loop)
    stats = Counter()
    t0 = time.time()
    mergedSamples = []
    blockedSamples = []

    for batchStart in range(0, len(remainingIdxs), CHECKPOINT_INTERVAL):
        batchEnd = min(batchStart + CHECKPOINT_INTERVAL, len(remainingIdxs))
        batchIdxs = remainingIdxs[batchStart:batchEnd]

        print(f"\n  --- Batch {batchStart//CHECKPOINT_INTERVAL + 1} "
              f"[{batchStart+1}-{batchEnd}/{len(remainingIdxs)}] ---")

        for i, idx in enumerate(batchIdxs):
            title = str(df.at[idx, "PMC_Title"]) if (
                "PMC_Title" in df.columns and
                pd.notna(df.at[idx, "PMC_Title"])) else ""
            abstract = str(df.at[idx, "PMC_Abstract"]) if (
                "PMC_Abstract" in df.columns and
                pd.notna(df.at[idx, "PMC_Abstract"])) else ""
            name = str(df.at[idx, "name"]) if pd.notna(
                df.at[idx, "name"]) else ""
            origSeq = str(df.at[idx, seqCol])

            for val in ("nan", "None"):
                if title == val: title = ""
                if abstract == val: abstract = ""
                if name == val: name = ""

            if not title and name:
                title = f"Compound: {name}"
            if not title:
                stats["skip_no_text"] += 1
                continue

            # Ollama 调用
            llmResult = callOllama(title, abstract, name, model=model)
            time.sleep(API_DELAY_SECONDS)

            row = {
                "row_idx": int(idx),
                "original_seq": origSeq[:80],
                "llm_sugars": "",
                "llm_bioact": "",
                "llm_targets": "",
                "consensus": origSeq,
                "fw_status": "",
            }

            if llmResult:
                sugars = llmResult.get("llm_sugar_sequence", [])
                if isinstance(sugars, str):
                    sugars = [s.strip() for s in sugars.split(",") if s.strip()]
                bioact = llmResult.get("bioactivities", [])
                if isinstance(bioact, str):
                    bioact = [bioact]
                targets = llmResult.get("targets_and_diseases", [])
                if isinstance(targets, str):
                    targets = [targets]

                row["llm_sugars"] = "; ".join(sugars) if sugars else ""
                row["llm_bioact"] = "; ".join(bioact) if bioact else ""
                row["llm_targets"] = "; ".join(targets) if targets else ""

                if sugars:
                    consensus, fwStatus = firewallMerge(origSeq, sugars)
                    row["consensus"] = consensus
                    row["fw_status"] = fwStatus
                    stats[fwStatus] += 1

                    if fwStatus == "MERGED" and len(mergedSamples) < 10:
                        mergedSamples.append(row.copy())
                    elif "BLOCKED" in fwStatus and len(blockedSamples) < 5:
                        blockedSamples.append(row.copy())
                else:
                    stats["llm_empty"] += 1
            else:
                stats["llm_fail"] += 1

            allResults.append(row)

            totalDone = batchStart + i + 1
            if totalDone % 25 == 0:
                elapsed = time.time() - t0
                rate = totalDone / elapsed if elapsed > 0 else 1
                eta = (len(remainingIdxs) - totalDone) / rate
                print(f"    [{totalDone:,}/{len(remainingIdxs):,}] "
                      f"merged={stats.get('MERGED',0)} "
                      f"blocked={sum(v for k,v in stats.items() if 'BLOCKED' in k)} "
                      f"fail={stats.get('llm_fail',0)} "
                      f"({elapsed:.0f}s, ETA {eta:.0f}s)")

        # 断点保存 (Checkpoint save)
        saveCheckpoint(allResults)
        print(f"    [SAVED] checkpoint: {len(allResults):,} rows")

    # ============================
    # 最终合并 (Final merge)
    # ============================
    print(f"\n" + "=" * 70)
    print(f"  Final Merge & Report")
    print(f"=" * 70)

    for col in ["LLM_Sugars", "LLM_Bioactivity", "LLM_Targets"]:
        if col not in df.columns:
            df[col] = ""

    mergedTotal = 0
    for r in allResults:
        idx = int(r["row_idx"])
        if idx not in df.index:
            continue
        if r.get("llm_sugars"):
            df.at[idx, "LLM_Sugars"] = r["llm_sugars"]
        if r.get("llm_bioact"):
            df.at[idx, "LLM_Bioactivity"] = r["llm_bioact"]
        if r.get("llm_targets"):
            df.at[idx, "LLM_Targets"] = r["llm_targets"]
        if r.get("fw_status") == "MERGED":
            df.at[idx, seqCol] = r["consensus"]
            mergedTotal += 1

    # 统计 (Stats)
    print(f"\n  Processing Stats:")
    for status, count in stats.most_common():
        print(f"    {status:30s}: {count:,}")

    print(f"\n  Consensus updated: {mergedTotal:,}")
    fuzzyAfter = df[seqCol].str.contains(
        r'\b(?:Hex|Pen|dHex|HexA|Non|Oct|Hept)\b',
        na=False, regex=True).sum()
    print(f"  Fuzzy remaining:   {fuzzyAfter:,}")

    # 样本展示
    if mergedSamples:
        print(f"\n  [APPROVED] Merged samples (first 3):")
        for s in mergedSamples[:3]:
            print(f"    Row {s['row_idx']}: {s['original_seq'][:50]}")
            print(f"      LLM: {s['llm_sugars']}")
            print(f"      ->   {s['consensus'][:50]}")

    if blockedSamples:
        print(f"\n  [BLOCKED] Firewall-rejected (first 2):")
        for s in blockedSamples[:2]:
            print(f"    Row {s['row_idx']}: {s['original_seq'][:50]}")
            print(f"      LLM: {s['llm_sugars']}")
            print(f"      Reason: {s['fw_status']}")

    # 保存 (Save)
    outPath = os.path.join(REPORT_DIR, "GlycoNP_LocalLLM_Final.csv")
    df.to_csv(outPath, index=False, encoding="utf-8-sig")
    print(f"\n  Output: {outPath}")
    print("=" * 70)


if __name__ == "__main__":
    main()
