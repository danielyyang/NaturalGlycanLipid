"""
Literature-Powered Sugar Rescue Pipeline
=========================================

Task 1: Europe PMC bulk fetch (Title + Abstract)
Task 2: Claude API structured extraction (sugar sequence + bioactivity)
Task 3: Safe merge protocol (zero-pollution firewall)
Task 4: HTML debug report generation

Usage:
  set ANTHROPIC_API_KEY=sk-ant-...
  python scripts/literature_sugar_rescue.py --limit 500
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
from typing import Any, Dict, List, Optional

import pandas as pd
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
REPORT_DIR = os.path.join(BASE_DIR, "reports")
INPUT_CSV = os.path.join(REPORT_DIR, "GlycoNP_ZeroPollu.csv")
OUTPUT_CSV = os.path.join(REPORT_DIR, "GlycoNP_LLM_Rescued.csv")

GENERIC_PATTERN = re.compile(r'\b(Hex|Pen|dHex|HexA|Non|Oct|Hept)\b')


# =====================================================================
# Task 1: Europe PMC 批量抓取 (Europe PMC Batch Fetch)
# =====================================================================

def fetchEuropePmc(doi: str) -> Dict[str, str]:
    """
    通过 Europe PMC REST API 抓取 DOI 对应的 Title + Abstract。
    Fetch Title + Abstract via Europe PMC REST API.

    比 Crossref 更高的 Abstract 覆盖率, 特别是 PubMed-indexed 文献。
    """
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
            title = entry.get("title", "")
            abstract = entry.get("abstractText", "")
            # 清理 HTML/JATS tags (Clean HTML/JATS tags)
            title = re.sub(r'<[^>]+>', '', str(title)).strip()
            abstract = re.sub(r'<[^>]+>', '', str(abstract)).strip()
            return {"title": title, "abstract": abstract}
    except Exception:
        return {"title": "", "abstract": ""}


def runTask1_europePmcFetch(
    df: pd.DataFrame,
    targetIdxs: list,
) -> pd.DataFrame:
    """
    批量抓取 Europe PMC 数据。
    Batch fetch Title + Abstract from Europe PMC.
    """
    print("\n" + "=" * 70)
    print("  Task 1: Europe PMC Batch Fetch")
    print("  Europe PMC 文献批量抓取")
    print("=" * 70)
    t0 = time.time()

    # 初始化列 (Initialize columns)
    for col in ["PMC_Title", "PMC_Abstract"]:
        if col not in df.columns:
            df[col] = ""

    fetched = 0
    withAbstract = 0
    errors = 0

    for i, idx in enumerate(targetIdxs):
        doiStr = str(df.at[idx, "dois"]).strip()
        # 取第一个 DOI (Take first DOI)
        firstDoi = doiStr.split("|")[0].strip().split(";")[0].strip()
        if not firstDoi or firstDoi in ("nan", "None", ""):
            continue

        # 如果已经抓过, 跳过 (Skip if already fetched)
        existingTitle = str(df.at[idx, "PMC_Title"]) if pd.notna(
            df.at[idx, "PMC_Title"]) else ""
        if existingTitle and existingTitle not in ("", "nan", "None"):
            fetched += 1
            if str(df.at[idx, "PMC_Abstract"]) not in ("", "nan", "None"):
                withAbstract += 1
            continue

        time.sleep(0.2)  # 防封禁 (Rate limit)

        result = fetchEuropePmc(firstDoi)
        if result["title"]:
            df.at[idx, "PMC_Title"] = result["title"]
            fetched += 1
            if result["abstract"]:
                df.at[idx, "PMC_Abstract"] = result["abstract"][:3000]
                withAbstract += 1
        else:
            errors += 1

        if (i + 1) % 50 == 0:
            elapsed = time.time() - t0
            print(f"    [{i+1}/{len(targetIdxs)}] fetched={fetched}, "
                  f"abstracts={withAbstract}, errors={errors} "
                  f"({elapsed:.0f}s)")

    elapsed = time.time() - t0
    print(f"\n  Task 1 Results:")
    print(f"    Total queried: {len(targetIdxs):,}")
    print(f"    Titles fetched: {fetched:,}")
    print(f"    With abstract: {withAbstract:,} "
          f"({withAbstract/max(fetched,1)*100:.1f}%)")
    print(f"    Errors: {errors:,}")
    print(f"    Time: {elapsed:.0f}s")

    return df


# =====================================================================
# Task 2: Claude API 序列翻译 (Claude API Sequence Translation)
# =====================================================================

SYSTEM_PROMPT = (
    "You are an expert natural product chemist. Your core task is to "
    "resolve ambiguous sugar sequences (like 'hexose' or 'pentose') "
    "into specific monosaccharides based on the provided Title and "
    "Abstract, and extract bioactivity data. Output strictly in valid "
    "JSON format."
)


def buildUserPrompt(title: str, abstract: str, fuzzySeq: str) -> str:
    """构建 Claude User Prompt。 / Build Claude User Prompt."""
    textBlock = f"Title: {title}"
    if abstract:
        textBlock += f"\nAbstract: {abstract}"

    return f"""Analyze the following literature metadata for a natural product:
<paper>
{textBlock}
</paper>

The current cheminformatics system identified the sugar chain vaguely as: {fuzzySeq}

Extract the precise biochemical information and output this exact JSON schema:
{{
  "compound_name": "Specific name of the compound if clearly matched, else null",
  "llm_sugar_sequence": "The specific sugar sequence derived from text (e.g., 'D-Glc', 'L-Rha'). Do NOT use generic terms. Return null if unclear.",
  "bioactivities": ["Crucial: Extract ANY mentioned biological activities, assays, IC50/MIC values, or pharmacological effects. Return [] if none."],
  "targets_and_diseases": ["Crucial: Extract ANY mentioned cell lines, proteins, pathogens, or diseases. Return [] if none."],
  "evidence_quote": "Exact short quote from the text proving the sugar sequence or bioactivity."
}}"""


def callClaudeApi(
    systemPrompt: str,
    userPrompt: str,
    apiKey: str,
    model: str = "claude-sonnet-4-20250514",
    maxTokens: int = 500,
) -> Optional[Dict]:
    """
    调用 Claude Messages API。
    Call Claude Messages API.

    返回解析后的 JSON dict, 或 None (失败时)。
    Returns parsed JSON dict, or None on failure.
    """
    url = "https://api.anthropic.com/v1/messages"
    headers = {
        "Content-Type": "application/json",
        "x-api-key": apiKey,
        "anthropic-version": "2023-06-01",
    }
    payload = json.dumps({
        "model": model,
        "max_tokens": maxTokens,
        "system": systemPrompt,
        "messages": [{"role": "user", "content": userPrompt}],
    }).encode("utf-8")

    try:
        req = urllib.request.Request(url, data=payload, headers=headers,
                                     method="POST")
        with urllib.request.urlopen(req, timeout=30) as resp:
            respData = json.loads(resp.read().decode("utf-8"))
            # 提取 content text (Extract content text)
            contentBlocks = respData.get("content", [])
            text = ""
            for block in contentBlocks:
                if block.get("type") == "text":
                    text += block.get("text", "")

            if not text:
                return None

            # 尝试解析 JSON (Try parsing JSON)
            # Claude 可能输出带 ```json 的代码块
            cleaned = text.strip()
            if cleaned.startswith("```"):
                cleaned = re.sub(r'^```\w*\n?', '', cleaned)
                cleaned = re.sub(r'\n?```$', '', cleaned)
                cleaned = cleaned.strip()

            return json.loads(cleaned)
    except urllib.error.HTTPError as e:
        if e.code == 429:
            time.sleep(5)
            return None
        return None
    except (json.JSONDecodeError, Exception):
        return None


def runTask2_claudeParsing(
    df: pd.DataFrame,
    targetIdxs: list,
    apiKey: str,
) -> pd.DataFrame:
    """
    Claude API 批量调用: 序列翻译 + 活性提取。
    Claude API batch: sequence translation + bioactivity extraction.
    """
    print("\n" + "=" * 70)
    print("  Task 2: Claude API Sequence Translation")
    print("  Claude API 序列精准翻译引擎")
    print("=" * 70)
    t0 = time.time()

    # 初始化列 (Initialize columns)
    for col in ["LLM_Inferred_Sugars", "LLM_Bioactivity_Profile",
                 "LLM_Compound_Name", "LLM_Evidence_Quote",
                 "LLM_Raw_JSON"]:
        if col not in df.columns:
            df[col] = ""

    parsed = 0
    sugarResolved = 0
    bioactivityFound = 0
    errors = 0
    sampleResults = []  # 存储样本用于展示

    # 只处理有 Title 的记录 (Only process records with Title)
    validIdxs = [
        idx for idx in targetIdxs
        if pd.notna(df.at[idx, "PMC_Title"]) and
        str(df.at[idx, "PMC_Title"]) not in ("", "nan", "None")
    ]
    print(f"  Records with PMC_Title: {len(validIdxs):,}")

    for i, idx in enumerate(validIdxs):
        title = str(df.at[idx, "PMC_Title"])
        abstract = str(df.at[idx, "PMC_Abstract"]) if pd.notna(
            df.at[idx, "PMC_Abstract"]) else ""
        if abstract in ("nan", "None"):
            abstract = ""
        fuzzySeq = str(df.at[idx, "Sugar_Sequence"])

        userPrompt = buildUserPrompt(title, abstract, fuzzySeq)

        time.sleep(0.3)  # Claude API 限流
        result = callClaudeApi(SYSTEM_PROMPT, userPrompt, apiKey)

        if result:
            parsed += 1
            rawJson = json.dumps(result, ensure_ascii=False)
            df.at[idx, "LLM_Raw_JSON"] = rawJson[:2000]

            # 提取糖序列 (Extract sugar sequence)
            llmSugar = result.get("llm_sugar_sequence")
            if llmSugar and llmSugar != "null" and str(llmSugar) != "None":
                df.at[idx, "LLM_Inferred_Sugars"] = str(llmSugar)
                sugarResolved += 1

            # 提取活性 (Extract bioactivity)
            bioacts = result.get("bioactivities", [])
            targets = result.get("targets_and_diseases", [])
            combined = []
            if bioacts:
                combined.extend(bioacts)
            if targets:
                combined.extend([f"[Target] {t}" for t in targets])
            if combined:
                df.at[idx, "LLM_Bioactivity_Profile"] = "; ".join(combined)
                bioactivityFound += 1

            # 提取化合物名 (Extract compound name)
            compName = result.get("compound_name")
            if compName and str(compName) != "None" and compName != "null":
                df.at[idx, "LLM_Compound_Name"] = str(compName)

            # 提取证据 (Extract evidence)
            evidence = result.get("evidence_quote")
            if evidence and str(evidence) != "None":
                df.at[idx, "LLM_Evidence_Quote"] = str(evidence)[:500]

            # 保存样本 (Save sample)
            sampleResults.append({
                "idx": int(idx),
                "title": title[:100],
                "abstract_len": len(abstract),
                "fuzzy_seq": fuzzySeq,
                "llm_result": result,
            })
        else:
            errors += 1

        if (i + 1) % 20 == 0:
            elapsed = time.time() - t0
            print(f"    [{i+1}/{len(validIdxs)}] parsed={parsed}, "
                  f"sugars={sugarResolved}, bioact={bioactivityFound}, "
                  f"errors={errors} ({elapsed:.0f}s)")

    elapsed = time.time() - t0
    print(f"\n  Task 2 Results:")
    print(f"    Total sent to Claude: {len(validIdxs):,}")
    print(f"    Successfully parsed: {parsed:,}")
    print(f"    Sugar sequences resolved: {sugarResolved:,}")
    print(f"    Bioactivity found: {bioactivityFound:,}")
    print(f"    Errors: {errors:,}")
    print(f"    Time: {elapsed:.0f}s")

    # 保存完整样本 (Save full samples to JSON)
    samplesPath = os.path.join(REPORT_DIR, "claude_pilot_results.json")
    with open(samplesPath, "w", encoding="utf-8") as f:
        json.dump(sampleResults, f, ensure_ascii=False, indent=2)
    print(f"  Samples saved: {samplesPath}")

    # 打印 TOP-3 最丰富的 bioactivity 结果
    bioSamples = [s for s in sampleResults
                  if s["llm_result"].get("bioactivities")]
    bioSamples.sort(
        key=lambda x: len(x["llm_result"].get("bioactivities", [])),
        reverse=True)

    print(f"\n  {'='*50}")
    print(f"  TOP-3 Richest Bioactivity JSON Results:")
    print(f"  {'='*50}")
    for j, s in enumerate(bioSamples[:3]):
        print(f"\n  --- Sample {j+1} ---")
        print(f"  Title: {s['title']}")
        print(f"  Fuzzy: {s['fuzzy_seq']}")
        print(json.dumps(s["llm_result"], ensure_ascii=False, indent=4))

    return df


# =====================================================================
# Task 3: 零污染安全合并 (Zero-Pollution Safe Merge)
# =====================================================================

def runTask3_safeMerge(df: pd.DataFrame) -> pd.DataFrame:
    """
    Pandas 防火墙: 生成 Consensus_Sugar_Sequence。

    条件 A (防污染): 原序列已精确 → 保持原样
    条件 B (精准反哺): 原序列含模糊标签 + LLM 有结果 → 替换并标记
    """
    print("\n" + "=" * 70)
    print("  Task 3: Zero-Pollution Safe Merge")
    print("  零污染安全合并 (Pandas 防火墙)")
    print("=" * 70)

    seqCol = "Sugar_Sequence"
    llmCol = "LLM_Inferred_Sugars"
    consensusCol = "Consensus_Sugar_Sequence"

    # 初始化 (Initialize)
    df[consensusCol] = df[seqCol].copy()

    total = len(df)
    protectedCount = 0   # 条件 A: 防污染守住的
    mergedCount = 0      # 条件 B: 成功合并的
    noLlmCount = 0       # LLM 没有结果的

    for idx in df.index:
        origSeq = str(df.at[idx, seqCol]) if pd.notna(
            df.at[idx, seqCol]) else ""
        llmSugar = str(df.at[idx, llmCol]) if pd.notna(
            df.at[idx, llmCol]) else ""

        if llmSugar in ("", "nan", "None"):
            noLlmCount += 1
            continue

        # 条件 A: 原序列是否已经精确? (Is original already precise?)
        hasGeneric = bool(GENERIC_PATTERN.search(origSeq))

        if not hasGeneric:
            # 防火墙: 原序列已精确, 不污染
            protectedCount += 1
            continue

        # 条件 B: 原序列含模糊标签 + LLM 有值 → 替换
        df.at[idx, consensusCol] = f"{llmSugar} (LLM-Rescued)"
        mergedCount += 1

    print(f"  Firewall Stats:")
    print(f"    Condition A (protected, already precise): {protectedCount:,}")
    print(f"    Condition B (merged from LLM): {mergedCount:,}")
    print(f"    No LLM data: {noLlmCount:,}")

    # 对比统计 (Comparison stats)
    hasGenericOrig = df[seqCol].str.contains(
        r'\bHex\b|\bPen\b|\bdHex\b|\bNon\b|\bOct\b|\bHept\b',
        na=False, regex=True).sum()
    hasGenericConsensus = df[consensusCol].str.contains(
        r'\bHex\b|\bPen\b|\bdHex\b|\bNon\b|\bOct\b|\bHept\b',
        na=False, regex=True).sum()

    print(f"\n  Generic labels in Sugar_Sequence: {hasGenericOrig:,}")
    print(f"  Generic labels in Consensus: {hasGenericConsensus:,}")
    print(f"  Reduction: -{hasGenericOrig - hasGenericConsensus:,}")

    return df


# =====================================================================
# Task 4: HTML Debug Report
# =====================================================================

def runTask4_htmlReport(
    df: pd.DataFrame,
    targetIdxs: list,
    outputPath: str,
) -> None:
    """生成 HTML 调试报告。 / Generate HTML debug report."""
    print("\n" + "=" * 70)
    print("  Task 4: HTML Debug Report")
    print("  HTML 调试报告生成")
    print("=" * 70)

    htmlParts = ["""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>GlycoNP LLM Rescue Debug Report</title>
<style>
body { font-family: 'Segoe UI', Arial, sans-serif; margin: 20px;
       background: #f5f5f5; }
h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
table { border-collapse: collapse; width: 100%; margin: 10px 0;
        background: white; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
th { background: #2c3e50; color: white; padding: 10px; text-align: left;
     font-size: 12px; }
td { padding: 8px 10px; border-bottom: 1px solid #eee; font-size: 11px;
     max-width: 300px; word-wrap: break-word; }
tr:hover { background: #f0f8ff; }
.fuzzy { background: #fff3cd; color: #856404; padding: 2px 4px;
         border-radius: 3px; font-weight: bold; }
.precise { background: #d4edda; color: #155724; padding: 2px 4px;
           border-radius: 3px; }
.llm { background: #d1ecf1; color: #0c5460; padding: 2px 4px;
       border-radius: 3px; }
.bioact { background: #f8d7da; color: #721c24; padding: 2px 4px;
          border-radius: 3px; font-size: 10px; }
.consensus { background: #e2e3f1; color: #383d6e; padding: 2px 4px;
             border-radius: 3px; font-weight: bold; }
.stat-box { display: inline-block; background: white; padding: 15px 25px;
            margin: 5px; border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
.stat-num { font-size: 24px; font-weight: bold; color: #2c3e50; }
.stat-label { font-size: 12px; color: #7f8c8d; }
</style>
</head>
<body>
<h1>GlycoNP LLM Rescue Debug Report</h1>
"""]

    # 统计摘要 (Stats summary)
    subset = df.loc[targetIdxs]
    totalRows = len(subset)
    hasLlm = (subset["LLM_Inferred_Sugars"].notna() &
              ~subset["LLM_Inferred_Sugars"].astype(str).isin(
                  ["", "nan", "None"])).sum()
    hasBioact = (subset["LLM_Bioactivity_Profile"].notna() &
                 ~subset["LLM_Bioactivity_Profile"].astype(str).isin(
                     ["", "nan", "None"])).sum()
    hasConsensus = (subset["Consensus_Sugar_Sequence"].astype(str
                    ).str.contains("LLM-Rescued", na=False)).sum()

    htmlParts.append(f"""
<div>
  <div class="stat-box"><div class="stat-num">{totalRows}</div>
    <div class="stat-label">Total Records</div></div>
  <div class="stat-box"><div class="stat-num">{hasLlm}</div>
    <div class="stat-label">LLM Sugar Resolved</div></div>
  <div class="stat-box"><div class="stat-num">{hasBioact}</div>
    <div class="stat-label">Bioactivity Extracted</div></div>
  <div class="stat-box"><div class="stat-num">{hasConsensus}</div>
    <div class="stat-label">Consensus Merged</div></div>
</div>
<br>
<table>
<tr>
  <th>#</th>
  <th>Name</th>
  <th>Sugar_Sequence (Original)</th>
  <th>LLM_Inferred_Sugars</th>
  <th>Consensus_Sugar_Sequence</th>
  <th>LLM_Bioactivity</th>
  <th>Evidence</th>
  <th>DOI</th>
</tr>
""")

    for i, idx in enumerate(targetIdxs):
        row = df.loc[idx]
        name = html.escape(str(row.get("name", ""))[:60])
        origSeq = str(row.get("Sugar_Sequence", ""))
        llmSugar = str(row.get("LLM_Inferred_Sugars", ""))
        consensus = str(row.get("Consensus_Sugar_Sequence", ""))
        bioact = str(row.get("LLM_Bioactivity_Profile", ""))
        evidence = str(row.get("LLM_Evidence_Quote", ""))
        doi = str(row.get("dois", ""))[:40]

        # 高亮模糊标签 (Highlight generic labels)
        origHtml = GENERIC_PATTERN.sub(
            r'<span class="fuzzy">\1</span>',
            html.escape(origSeq))

        llmHtml = html.escape(llmSugar) if llmSugar not in (
            "", "nan", "None") else "-"
        if llmHtml != "-":
            llmHtml = f'<span class="llm">{llmHtml}</span>'

        consHtml = html.escape(consensus)
        if "LLM-Rescued" in consensus:
            consHtml = f'<span class="consensus">{consHtml}</span>'
        elif not GENERIC_PATTERN.search(consensus):
            consHtml = f'<span class="precise">{consHtml}</span>'

        bioHtml = ""
        if bioact and bioact not in ("", "nan", "None"):
            bioHtml = f'<span class="bioact">{html.escape(bioact[:120])}</span>'

        evidHtml = html.escape(evidence[:100]) if evidence not in (
            "", "nan", "None") else ""

        htmlParts.append(f"""<tr>
  <td>{i+1}</td>
  <td>{name}</td>
  <td>{origHtml}</td>
  <td>{llmHtml}</td>
  <td>{consHtml}</td>
  <td>{bioHtml}</td>
  <td style="font-size:10px;color:#666">{evidHtml}</td>
  <td style="font-size:10px">{html.escape(doi)}</td>
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
        description="Literature Sugar Rescue Pipeline")
    parser.add_argument("--limit", type=int, default=500,
                        help="Number of records to process")
    parser.add_argument("--model", type=str,
                        default="claude-sonnet-4-20250514",
                        help="Claude model name")
    args = parser.parse_args()

    print("=" * 70)
    print("  Literature-Powered Sugar Rescue Pipeline")
    print("  文献驱动的糖链序列精准抢救管线")
    print("=" * 70)

    # API key 检查 (API key check)
    apiKey = os.environ.get("ANTHROPIC_API_KEY", "")
    if not apiKey:
        print("\n  [CRITICAL] ANTHROPIC_API_KEY not set!")
        print("  Please set it: set ANTHROPIC_API_KEY=sk-ant-...")
        print("  Continuing with Task 1 (PMC fetch) only...")

    if not os.path.exists(INPUT_CSV):
        print(f"\n  [ERROR] Input not found: {INPUT_CSV}")
        return

    df = pd.read_csv(INPUT_CSV, dtype=str, low_memory=False,
                     encoding="utf-8-sig")
    total = len(df)
    print(f"  Loaded: {total:,} rows")

    # 选择目标集: 有 DOI 的记录 (Select targets: records with DOI)
    # 混合采样: 含模糊糖 + 精确糖 (以测试防火墙)
    hasDoi = (
        df["dois"].notna() &
        (~df["dois"].astype(str).isin(["", "nan", "None"]))
    )
    hasGeneric = df["Sugar_Sequence"].str.contains(
        r'\b(Hex|Pen|dHex|HexA|Non|Oct|Hept)\b', na=False, regex=True)

    # 优先选模糊糖记录, 再补精确糖记录
    fuzzyWithDoi = df.index[hasDoi & hasGeneric].tolist()
    preciseWithDoi = df.index[hasDoi & ~hasGeneric].tolist()

    # 80% 模糊 + 20% 精确 (混合测试防火墙)
    nFuzzy = min(int(args.limit * 0.8), len(fuzzyWithDoi))
    nPrecise = min(args.limit - nFuzzy, len(preciseWithDoi))

    import random
    random.seed(42)
    targetIdxs = (random.sample(fuzzyWithDoi, nFuzzy) +
                  random.sample(preciseWithDoi, nPrecise))
    random.shuffle(targetIdxs)

    print(f"  Target records: {len(targetIdxs):,} "
          f"(fuzzy={nFuzzy}, precise={nPrecise})")

    # ============ Task 1: Europe PMC Fetch ============
    df = runTask1_europePmcFetch(df, targetIdxs)

    # ============ Task 2: Claude API ============
    if apiKey:
        df = runTask2_claudeParsing(df, targetIdxs, apiKey)
    else:
        print("\n  [SKIP] Task 2: No API key, skipping Claude parsing")
        # 确保 LLM 列存在 (Ensure LLM columns exist for downstream tasks)
        for col in ["LLM_Inferred_Sugars", "LLM_Bioactivity_Profile",
                     "LLM_Compound_Name", "LLM_Evidence_Quote",
                     "LLM_Raw_JSON"]:
            if col not in df.columns:
                df[col] = ""

    # ============ Task 3: Safe Merge ============
    df = runTask3_safeMerge(df)

    # ============ Task 4: HTML Report ============
    htmlPath = os.path.join(REPORT_DIR,
                            "debug_sample_500_llm_rescued.html")
    runTask4_htmlReport(df, targetIdxs, htmlPath)

    # ============ Save ============
    df.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")
    print(f"\n  Output: {OUTPUT_CSV}")
    print("=" * 70)


if __name__ == "__main__":
    main()
