"""
sample_1k_stress_test.py — 从 GlycoNP v12 管线输出中分层抽样 1,000 个含糖分子
(Stratified 1K Sampling from GlycoNP v12 Pipeline Output)
================================================================

数据源: D:\Glycan_Database\reports\GlycoNP_Deep_Enriched_v12_backup_bond.csv
注意: 旧管线的 Glycan_SMILES/Sugar_Sequence 字段可能有误,
      本脚本仅取 canonical_smiles 重新跑识别

分层策略 (Stratification Strategy):
  - Tier S (Simple glycosides, 1 sugar): 400
  - Tier M (Multi-sugar, 2-3 sugars): 300
  - Tier C (Complex, 4+ sugars): 200
  - Tier N (No sugar detected by old pipeline): 100

Usage:
  python scripts/sample_1k_stress_test.py
"""
import os
import sys
import json
import time
import random
from typing import Dict, List

import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from rdkit import Chem
from lib.glycan_topology import find_mapped_sugar_units
from lib.monosaccharide_identifier import identify_monosaccharide_v10


# =====================================================================
# 配置 (Configuration)
# =====================================================================
CSV_PATH = os.path.join(
    os.path.dirname(__file__), "..",
    "reports", "GlycoNP_Deep_Enriched_v12_backup_bond.csv"
)
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "..", "reports")
SAMPLE_SIZES = {
    "S": 400,   # 简单糖苷: 1 个糖
    "M": 300,   # 多糖: 2-3 个糖
    "C": 200,   # 复杂: 4+ 个糖
    "N": 100,   # 旧管线未检出糖 (验证假阴性)
}
TOTAL = sum(SAMPLE_SIZES.values())
RANDOM_SEED = 42


def classifyBySugarCount(row: pd.Series) -> str:
    """根据 Total_Sugar_Count 列分层。
    Classify by Total_Sugar_Count column from old pipeline.

    设计意图: Total_Sugar_Count 是旧管线的数值型统计,
    比解析 Sugar_Sequence 文本更可靠。
    """
    try:
        count = int(row.get("Total_Sugar_Count", 0))
    except (ValueError, TypeError):
        count = 0

    if count == 0:
        return "N"   # 旧管线未检出
    elif count == 1:
        return "S"   # 简单糖苷
    elif count <= 3:
        return "M"   # 多糖
    else:
        return "C"   # 复杂


def runStressTest() -> None:
    print("=" * 70)
    print("  1K Stress Test — Stratified Sampling from GlycoNP v12")
    print("=" * 70)

    # === 1. 加载 CSV ===
    print("\n  Loading CSV...")
    df = pd.read_csv(CSV_PATH, low_memory=False)
    print("  Total rows: %d" % len(df))
    print("  Columns: %d" % len(df.columns))

    # 确保有 canonical_smiles
    if "canonical_smiles" not in df.columns:
        print("  ERROR: 'canonical_smiles' column not found!")
        return

    # 过滤无效 SMILES
    df = df[df["canonical_smiles"].notna()].copy()
    df = df[df["canonical_smiles"].str.len() > 0].copy()
    print("  Valid SMILES rows: %d" % len(df))

    # === 2. 分层 ===
    print("\n  Stratifying by old pipeline Sugar_Sequence...")
    df["tier"] = df.apply(classifyBySugarCount, axis=1)
    tierCounts = df["tier"].value_counts().to_dict()
    print("  Tier distribution:")
    for t in ["S", "M", "C", "N"]:
        print("    %s: %d" % (t, tierCounts.get(t, 0)))

    # === 3. 分层抽样 ===
    random.seed(RANDOM_SEED)
    sampled = []
    for tier, targetN in SAMPLE_SIZES.items():
        pool = df[df["tier"] == tier]
        actualN = min(targetN, len(pool))
        if actualN < targetN:
            print("  WARNING: Tier %s only has %d (need %d)" % (tier, actualN, targetN))
        sample = pool.sample(n=actualN, random_state=RANDOM_SEED)
        sampled.append(sample)
        print("  Sampled %d from Tier %s" % (actualN, tier))

    sampleDf = pd.concat(sampled, ignore_index=True)
    print("\n  Total sampled: %d" % len(sampleDf))

    # === 4. 跑新管线 ===
    print("\n  Running new pipeline on %d molecules..." % len(sampleDf))

    results = []
    passCount = 0
    failCount = 0
    errorCount = 0
    totalTime = 0.0
    furanoseCount = 0
    pseudosugarCount = 0

    for rowIdx, row in sampleDf.iterrows():
        smi = row["canonical_smiles"]
        identifier = row.get("identifier", "?")
        oldSeq = str(row.get("Sugar_Sequence", ""))
        tier = row["tier"]

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            results.append({
                "identifier": identifier, "tier": tier,
                "status": "ERROR", "detail": "SMILES parse failed",
                "old_seq": oldSeq, "new_sugars": [], "latency_ms": 0,
            })
            errorCount += 1
            continue

        t0 = time.perf_counter()

        try:
            units = find_mapped_sugar_units(mol)
            newSugars = []

            for u in units:
                ra = u.get("ring_atoms", [])
                if not ra:
                    continue

                # Pseudosugar check
                ringO = sum(1 for i in ra if mol.GetAtomWithIdx(i).GetAtomicNum() == 8)
                if ringO == 0:
                    newSugars.append("[PSEUDOSUGAR]")
                    pseudosugarCount += 1
                    continue

                if len(ra) == 5:
                    furanoseCount += 1

                result = identify_monosaccharide_v10(mol, ra)
                sugarName = result[0] if isinstance(result, tuple) else result
                newSugars.append(sugarName)

            latencyMs = (time.perf_counter() - t0) * 1000
            totalTime += latencyMs

            # 健康检查: 引擎是否崩溃、超时、或输出异常
            hasCrash = False
            for s in newSugars:
                if "ERROR" in str(s) or s is None:
                    hasCrash = True

            status = "OK" if not hasCrash else "WARN"
            if hasCrash:
                failCount += 1
            else:
                passCount += 1

            results.append({
                "identifier": identifier, "tier": tier,
                "status": status, "old_seq": oldSeq,
                "new_sugars": newSugars, "num_sugars": len(newSugars),
                "latency_ms": round(latencyMs, 1),
            })

        except Exception as e:
            latencyMs = (time.perf_counter() - t0) * 1000
            totalTime += latencyMs
            results.append({
                "identifier": identifier, "tier": tier,
                "status": "CRASH", "detail": str(e)[:100],
                "old_seq": oldSeq, "new_sugars": [],
                "latency_ms": round(latencyMs, 1),
            })
            errorCount += 1

        # 每 100 个打印进度
        processed = passCount + failCount + errorCount
        if processed % 100 == 0:
            print("    ... %d / %d (%.1f%%)" % (processed, len(sampleDf),
                                                 processed / len(sampleDf) * 100))

    # === 5. 汇总报告 ===
    total = passCount + failCount + errorCount
    avgLatency = totalTime / total if total else 0

    report = []
    report.append("=" * 70)
    report.append("  1K Stress Test Report")
    report.append("=" * 70)
    report.append("")
    report.append("  Total:      %d" % total)
    report.append("  OK:         %d  (%.1f%%)" % (passCount, passCount / total * 100))
    report.append("  WARN:       %d  (%.1f%%)" % (failCount, failCount / total * 100))
    report.append("  ERROR/CRASH:%d  (%.1f%%)" % (errorCount, errorCount / total * 100))
    report.append("")
    report.append("  Avg Latency: %.1f ms/mol" % avgLatency)
    report.append("  Total Time:  %.1f s" % (totalTime / 1000))
    report.append("  Max Latency: %.1f ms" % max(r["latency_ms"] for r in results))
    report.append("")
    report.append("  Furanose (5-ring): %d" % furanoseCount)
    report.append("  Pseudosugar:       %d" % pseudosugarCount)
    report.append("")

    # 按 Tier 统计
    report.append("-" * 70)
    report.append("  BY TIER")
    report.append("-" * 70)
    for tier in ["S", "M", "C", "N"]:
        tierResults = [r for r in results if r["tier"] == tier]
        tierOK = sum(1 for r in tierResults if r["status"] == "OK")
        tierErr = sum(1 for r in tierResults if r["status"] in ("ERROR", "CRASH"))
        tierTotal = len(tierResults)
        if tierTotal:
            report.append("  Tier %s: %d total, %d OK (%.1f%%), %d errors" % (
                tier, tierTotal, tierOK, tierOK / tierTotal * 100, tierErr))

    # 错误/崩溃列表
    crashes = [r for r in results if r["status"] in ("ERROR", "CRASH")]
    if crashes:
        report.append("")
        report.append("-" * 70)
        report.append("  ERRORS & CRASHES (first 20)")
        report.append("-" * 70)
        for r in crashes[:20]:
            report.append("  [%s] %s: %s" % (r["tier"], r["identifier"],
                                               r.get("detail", "unknown")))

    report.append("")
    report.append("=" * 70)

    reportText = "\n".join(report)
    print("\n" + reportText)

    # 保存报告
    reportPath = os.path.join(OUTPUT_DIR, "stress_test_1k_report.txt")
    with open(reportPath, "w", encoding="utf-8") as f:
        f.write(reportText)

    # 保存详细 JSON
    jsonPath = os.path.join(OUTPUT_DIR, "stress_test_1k_results.json")
    with open(jsonPath, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)

    print("\n  Report: %s" % reportPath)
    print("  JSON:   %s" % jsonPath)


if __name__ == "__main__":
    runStressTest()
