"""
Opus LLM Sugar Sequence Cleanup — Industrial-Grade Pipeline
============================================================

终极 Opus LLM 清扫脚本 (工业级容灾机制):

Task 1: 强壮的 API 调用引擎 (Exponential Backoff + JSON Sanitization)
Task 2: 断点续传 (Checkpoint every 500 rows)
Task 3: v3 碳数防火墙合并

Usage:
  python scripts/opus_llm_cleanup.py
  python scripts/opus_llm_cleanup.py --dry-run     # 仅打印配置, 不执行
  python scripts/opus_llm_cleanup.py --limit 100    # 限制处理数量

Configuration:
  Set environment variables or edit constants below:
    OPUS_API_URL   = http://localhost:8080/v1/messages
    OPUS_API_KEY   = your-api-key-here
    OPUS_MODEL     = claude-sonnet-4-20250514
"""
import json
import os
import re
import sys
import time
import traceback
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

# API 端点 — 预留变量供用户配置
OPUS_API_URL = os.environ.get(
    "OPUS_API_URL", "http://localhost:8080/v1/messages")
OPUS_API_KEY = os.environ.get(
    "OPUS_API_KEY", "your-api-key-here")
OPUS_MODEL = os.environ.get(
    "OPUS_MODEL", "claude-sonnet-4-20250514")

# 断点续传文件 (Checkpoint file)
CHECKPOINT_PATH = os.path.join(REPORT_DIR, "GlycoNP_LLM_Checkpoint.csv")

# 批次大小 (Checkpoint interval)
CHECKPOINT_INTERVAL = 500

# API 调用间隔 (Rate limiting)
API_DELAY_SECONDS = 0.5

# 最大重试次数
MAX_RETRIES = 5

GENERIC_PATTERN = re.compile(r'\b(Hex|Pen|dHex|HexA|Non|Oct|Hept)\b')

# =====================================================================
# 碳类映射 (Carbon-Class Map)
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
# Task 1: 强壮的 Opus API 客户端
# =====================================================================

OPUS_SYSTEM_PROMPT = """You are a glycochemistry expert analyzing natural product literature.

Given the title and abstract of a scientific paper about a glycosylated natural product, extract:
1. The specific sugar monomers mentioned (use standard abbreviations)
2. Bioactivity information
3. Biological targets or diseases

Output ONLY valid JSON (no markdown, no explanation):
{
  "compound_name": "name of the compound",
  "llm_sugar_sequence": ["D-Glc", "L-Rha"],
  "bioactivities": ["cytotoxic", "antibacterial"],
  "targets_and_diseases": ["A549 cells", "S. aureus"]
}

Standard sugar abbreviations:
- Hexoses: D-Glc, D-Gal, D-Man, L-Glc, L-Gal, D-Tal, D-All, D-Fru
- Pentoses: D-Xyl, L-Ara, D-Rib, D-Api, D-Lyx
- Deoxysugars: L-Rha, L-Fuc, D-Fuc, L-Ole, D-Cym, D-The, D-Dig
- Amino sugars: D-GlcN, D-GalN, D-GlcNAc, D-GalNAc
- Uronic acids: D-GlcA, D-GalA, D-ManA
- Special: Neu5Ac, KDO

Rules:
- Only include sugars EXPLICITLY mentioned or clearly identifiable in the text
- If the paper discusses multiple compounds, focus on the PRIMARY glycoside
- Return empty arrays if information is not found
- Output ONLY the JSON object, nothing else"""


def sanitizeJsonResponse(rawText: str) -> Optional[Dict]:
    """
    清理并解析 Opus 返回的 JSON 响应。
    Sanitize and parse Opus API JSON response.

    处理常见格式问题:
    - 剥离 markdown 代码块 (```json ... ```)
    - 剥离前后缀文字
    - 处理单引号 vs 双引号
    """
    if not rawText:
        return None

    text = rawText.strip()

    # 去除 markdown 代码块 (Strip markdown code blocks)
    text = re.sub(r'^```(?:json)?\s*\n?', '', text, flags=re.MULTILINE)
    text = re.sub(r'\n?```\s*$', '', text, flags=re.MULTILINE)
    text = text.strip()

    # 尝试提取 JSON 对象 (Try to extract JSON object)
    jsonMatch = re.search(r'\{[^{}]*(?:\{[^{}]*\}[^{}]*)*\}', text,
                          re.DOTALL)
    if jsonMatch:
        text = jsonMatch.group(0)

    # 尝试解析 (Try parsing)
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        pass

    # 后备: 替换单引号 (Fallback: replace single quotes)
    try:
        fixed = text.replace("'", '"')
        return json.loads(fixed)
    except json.JSONDecodeError:
        pass

    return None


def callOpusApi(
    title: str,
    abstract: str,
    compoundName: str = "",
) -> Optional[Dict]:
    """
    带指数退避重试的 Opus API 调用。
    Call Opus API with exponential backoff retry.

    Returns:
        解析后的 JSON 字典, 或 None (失败时)
    """
    userPrompt = f"""Paper Title: {title}

Abstract: {abstract[:2000] if abstract else 'Not available'}

Compound Name (from database): {compoundName if compoundName else 'Unknown'}

Extract the sugar monomers and bioactivity from this paper."""

    payload = {
        "model": OPUS_MODEL,
        "max_tokens": 400,
        "temperature": 0.1,
        "system": OPUS_SYSTEM_PROMPT,
        "messages": [
            {"role": "user", "content": userPrompt},
        ],
    }

    headers = {
        "Content-Type": "application/json",
        "x-api-key": OPUS_API_KEY,
        "anthropic-version": "2023-06-01",
    }

    for attempt in range(MAX_RETRIES):
        try:
            resp = requests.post(
                OPUS_API_URL,
                json=payload,
                headers=headers,
                timeout=60,
            )

            # 成功 (Success)
            if resp.status_code == 200:
                respData = resp.json()
                # Anthropic Messages API 格式
                content = ""
                if "content" in respData:
                    for block in respData["content"]:
                        if block.get("type") == "text":
                            content += block.get("text", "")
                elif "completion" in respData:
                    content = respData["completion"]
                elif "choices" in respData:
                    # OpenAI 兼容格式
                    content = respData["choices"][0].get(
                        "message", {}).get("content", "")

                return sanitizeJsonResponse(content)

            # 可重试错误 (Retryable errors)
            elif resp.status_code in (429, 500, 502, 503, 529):
                waitTime = (2 ** attempt) + 1  # 1, 3, 5, 9, 17 秒
                retryAfter = resp.headers.get("retry-after")
                if retryAfter:
                    try:
                        waitTime = max(int(retryAfter), waitTime)
                    except ValueError:
                        pass
                print(f"      [RETRY {attempt+1}/{MAX_RETRIES}] "
                      f"HTTP {resp.status_code}, wait {waitTime}s")
                time.sleep(waitTime)
                continue

            # 不可重试错误 (Non-retryable errors)
            else:
                print(f"      [ERROR] HTTP {resp.status_code}: "
                      f"{resp.text[:200]}")
                return None

        except requests.exceptions.Timeout:
            waitTime = (2 ** attempt) + 1
            print(f"      [TIMEOUT] Attempt {attempt+1}/{MAX_RETRIES}, "
                  f"wait {waitTime}s")
            time.sleep(waitTime)

        except requests.exceptions.ConnectionError:
            waitTime = (2 ** attempt) + 2
            print(f"      [CONN ERROR] Attempt {attempt+1}/{MAX_RETRIES}, "
                  f"wait {waitTime}s")
            time.sleep(waitTime)

        except Exception as e:
            print(f"      [ERROR] {type(e).__name__}: {str(e)[:100]}")
            return None

    print(f"      [FAILED] Max retries exhausted")
    return None


# =====================================================================
# Task 3: v3 碳数防火墙合并
# =====================================================================

def firewallMerge(origSeq: str, llmSugars: List[str]) -> Tuple[str, str]:
    """
    v3 防火墙: 碳类一致性 + 多重候选保护锁。
    v3 Firewall: carbon-class consistency + multi-candidate backoff.
    """
    if not llmSugars:
        return origSeq, "NO_LLM"

    if not GENERIC_PATTERN.search(str(origSeq)):
        return origSeq, "NO_GENERIC"

    genericTokens = GENERIC_PATTERN.findall(str(origSeq))
    tokenCounts = Counter(genericTokens)

    classToSugars = {}
    for sugar in llmSugars:
        sugarClass = SUGAR_TO_CLASS.get(sugar)
        if sugarClass:
            classToSugars.setdefault(sugarClass, set()).add(sugar)

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
        if token not in {"Hex", "Pen", "dHex", "HexA", "Non", "Oct", "Hept"}:
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
    else:
        return origSeq, "NO_COMPATIBLE"


# =====================================================================
# Task 2: 主循环 — 断点续传
# =====================================================================

def loadCheckpoint() -> pd.DataFrame:
    """
    加载断点文件 (如存在)。
    Load checkpoint file if it exists.
    """
    if os.path.exists(CHECKPOINT_PATH):
        ckpt = pd.read_csv(CHECKPOINT_PATH, dtype=str, encoding="utf-8-sig")
        print(f"  [CHECKPOINT] Loaded {len(ckpt):,} processed rows "
              f"from {os.path.basename(CHECKPOINT_PATH)}")
        return ckpt
    return pd.DataFrame()


def saveCheckpoint(results: List[Dict]) -> None:
    """
    保存处理结果到断点文件。
    Save processing results to checkpoint file.
    """
    if not results:
        return
    ckptDf = pd.DataFrame(results)
    ckptDf.to_csv(CHECKPOINT_PATH, index=False, encoding="utf-8-sig")


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Opus LLM Sugar Sequence Cleanup")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print config and exit")
    parser.add_argument("--limit", type=int, default=0,
                        help="Limit number of records to process")
    args = parser.parse_args()

    print("=" * 70)
    print("  Opus LLM Sugar Sequence Cleanup")
    print("  Opus LLM 糖序列终极清扫")
    print("=" * 70)
    print(f"  API URL:   {OPUS_API_URL}")
    print(f"  Model:     {OPUS_MODEL}")
    print(f"  API Key:   {'***' + OPUS_API_KEY[-4:] if len(OPUS_API_KEY) > 8 else '[NOT SET]'}")
    print(f"  Checkpoint: {CHECKPOINT_PATH}")
    print(f"  Interval:  every {CHECKPOINT_INTERVAL} rows")

    if args.dry_run:
        print("\n  [DRY RUN] Config verified. Exiting.")
        return

    # 读取最新数据 (Load latest data)
    inputCsv = os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched_v2.csv")
    if not os.path.exists(inputCsv):
        inputCsv = os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched.csv")

    df = pd.read_csv(inputCsv, dtype=str, low_memory=False,
                     encoding="utf-8-sig")
    print(f"\n  Loaded: {len(df):,} rows from {os.path.basename(inputCsv)}")

    seqCol = "Consensus_Sugar_Sequence"
    if seqCol not in df.columns:
        seqCol = "Sugar_Sequence"

    # 筛选目标记录 (Select target records)
    hasFuzzy = df[seqCol].str.contains(
        r'\b(?:Hex|Pen|dHex|HexA|Non|Oct|Hept)\b', na=False, regex=True)
    hasDoi = (df["dois"].notna() &
              ~df["dois"].astype(str).isin(["", "nan", "None"]))
    hasText = False
    for col in ["PMC_Title", "name"]:
        if col in df.columns:
            hasText = hasText | (df[col].notna() &
                                 ~df[col].astype(str).isin(["", "nan", "None"]))

    targetMask = hasFuzzy & (hasDoi | hasText)
    targetIdxs = df.index[targetMask].tolist()
    print(f"  Target records (fuzzy + text): {len(targetIdxs):,}")

    if args.limit > 0:
        targetIdxs = targetIdxs[:args.limit]
        print(f"  Limited to: {len(targetIdxs):,}")

    # 加载断点 (Load checkpoint)
    ckptDf = loadCheckpoint()
    processedRows = set()
    allResults = []
    if not ckptDf.empty and "row_idx" in ckptDf.columns:
        processedRows = set(ckptDf["row_idx"].astype(int).tolist())
        allResults = ckptDf.to_dict("records")
        print(f"  Resuming: {len(processedRows):,} already processed")

    # 过滤已处理 (Filter already processed)
    remainingIdxs = [i for i in targetIdxs if i not in processedRows]
    print(f"  Remaining: {len(remainingIdxs):,}")

    if not remainingIdxs:
        print("  All records already processed!")
    else:
        # 连接测试 (Connection test)
        print(f"\n  Testing API connection...")
        testResp = callOpusApi("Test connection", "This is a test.")
        if testResp is None:
            print("  [ERROR] Cannot connect to Opus API!")
            print(f"  Please verify:")
            print(f"    1. API URL: {OPUS_API_URL}")
            print(f"    2. API Key is valid")
            print(f"    3. Service is running")
            print(f"\n  To configure, set environment variables:")
            print(f"    $env:OPUS_API_URL = 'http://your-api-endpoint'")
            print(f"    $env:OPUS_API_KEY = 'your-key'")
            return
        else:
            print(f"  API connection OK!")

    # 统计器 (Counters)
    stats = Counter()
    t0 = time.time()

    for batchStart in range(0, len(remainingIdxs), CHECKPOINT_INTERVAL):
        batchEnd = min(batchStart + CHECKPOINT_INTERVAL, len(remainingIdxs))
        batchIdxs = remainingIdxs[batchStart:batchEnd]

        print(f"\n  --- Batch {batchStart//CHECKPOINT_INTERVAL + 1} "
              f"[{batchStart+1}-{batchEnd}/{len(remainingIdxs)}] ---")

        for i, idx in enumerate(batchIdxs):
            # 收集文本 (Gather text)
            title = str(df.at[idx, "PMC_Title"]) if pd.notna(
                df.at[idx, "PMC_Title"]) else ""
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

            # 构造查询文本 (Build query text)
            if not title and name:
                title = f"Compound: {name}"
            if not title:
                stats["skip_no_text"] += 1
                continue

            # API调用 (Call API)
            llmResult = callOpusApi(title, abstract, name)
            time.sleep(API_DELAY_SECONDS)

            rowResult = {
                "row_idx": int(idx),
                "original_seq": origSeq[:80],
                "llm_raw": "",
                "llm_sugars": "",
                "llm_bioact": "",
                "llm_targets": "",
                "consensus": origSeq,
                "fw_status": "",
            }

            if llmResult:
                llmSugars = llmResult.get("llm_sugar_sequence", [])
                llmBioact = llmResult.get("bioactivities", [])
                llmTargets = llmResult.get("targets_and_diseases", [])

                rowResult["llm_raw"] = json.dumps(llmResult, ensure_ascii=False)
                rowResult["llm_sugars"] = "; ".join(llmSugars) if llmSugars else ""
                rowResult["llm_bioact"] = "; ".join(llmBioact) if llmBioact else ""
                rowResult["llm_targets"] = "; ".join(llmTargets) if llmTargets else ""

                # v3 防火墙合并 (Firewall merge)
                if llmSugars:
                    consensus, fwStatus = firewallMerge(origSeq, llmSugars)
                    rowResult["consensus"] = consensus
                    rowResult["fw_status"] = fwStatus
                    stats[fwStatus] += 1
                else:
                    stats["llm_empty"] += 1
            else:
                stats["api_fail"] += 1

            allResults.append(rowResult)

            # 每 50 条打印进度
            totalDone = batchStart + i + 1
            if totalDone % 50 == 0:
                elapsed = time.time() - t0
                rate = totalDone / elapsed if elapsed > 0 else 0
                eta = (len(remainingIdxs) - totalDone) / rate if rate > 0 else 0
                print(f"    [{totalDone:,}/{len(remainingIdxs):,}] "
                      f"merged={stats.get('MERGED',0)} "
                      f"blocked={sum(v for k,v in stats.items() if 'BLOCKED' in k)} "
                      f"({elapsed:.0f}s, ETA {eta:.0f}s)")

        # 保存断点 (Save checkpoint)
        saveCheckpoint(allResults)
        print(f"    [CHECKPOINT] Saved {len(allResults):,} rows")

    # ============================
    # 最终合并回主表 (Final merge)
    # ============================
    print(f"\n" + "=" * 70)
    print(f"  Final Merge & Statistics")
    print(f"=" * 70)

    # 初始化 LLM 列 (Initialize LLM columns)
    for col in ["LLM_Sugars", "LLM_Bioactivity", "LLM_Targets"]:
        if col not in df.columns:
            df[col] = ""

    mergedCount = 0
    for r in allResults:
        idx = int(r["row_idx"])
        if idx not in df.index:
            continue

        # 写入 LLM 提取结果
        if r.get("llm_sugars"):
            df.at[idx, "LLM_Sugars"] = r["llm_sugars"]
        if r.get("llm_bioact"):
            df.at[idx, "LLM_Bioactivity"] = r["llm_bioact"]
        if r.get("llm_targets"):
            df.at[idx, "LLM_Targets"] = r["llm_targets"]

        # 更新 Consensus
        if r.get("fw_status") == "MERGED":
            df.at[idx, seqCol] = r["consensus"]
            mergedCount += 1

    # 统计
    print(f"\n  Processing Stats:")
    for status, count in stats.most_common():
        print(f"    {status:30s}: {count:,}")

    print(f"\n  Consensus sequences updated: {mergedCount:,}")

    # 模糊标签前后对比
    fuzzyAfter = df[seqCol].str.contains(
        r'\b(?:Hex|Pen|dHex|HexA|Non|Oct|Hept)\b',
        na=False, regex=True).sum()
    print(f"  Fuzzy labels remaining:      {fuzzyAfter:,}")

    # 保存最终文件 (Save final output)
    outputPath = os.path.join(REPORT_DIR, "GlycoNP_Opus_Final.csv")
    df.to_csv(outputPath, index=False, encoding="utf-8-sig")
    print(f"\n  Output: {outputPath}")

    print("\n" + "=" * 70)
    print("  Opus LLM Cleanup Complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
