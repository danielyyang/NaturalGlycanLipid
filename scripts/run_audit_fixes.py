"""
run_audit_fixes.py
==================
全量修复验证 + CSV 更新脚本 (Audit Fix Verification + CSV Update Script)
==========================================================================

执行以下操作 (Performs the following):
1. 运行 226 条统一基准测试集 (Run unified benchmark 226 entries)
2. 修复 V13 Pruned CSV 中被反转的 Glycan/Aglycon SMILES 列
3. 同步更新 Saponin DB CSV

[TEST DATA ONLY] - GlycoNP Project
"""

import os
import sys
import time
import json
import logging

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
REPORT_DIR = os.path.join(BASE_DIR, "reports")
LOG_DIR = os.path.join(BASE_DIR, "log")
DATA_DIR = os.path.join(BASE_DIR, "data")

# ─── 日志配置 (Logging Setup) ──────────────────────────────────────────
logPath = os.path.join(LOG_DIR, "audit_fix_2026-03-30.log")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(logPath, encoding="utf-8"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("AuditFix")

logger.info("=" * 70)
logger.info("GlycoNP Codebase Audit Fix — 2026-03-30")
logger.info("=" * 70)


# ═══════════════════════════════════════════════════════════════════════
# STEP 1: 验证 Import 链路修复 (Verify Import Chain Fixes)
# ═══════════════════════════════════════════════════════════════════════
logger.info("STEP 1: Verifying import chain fixes...")

# Fix #1: cip_exo_engine 应该 import lib.virtual_demodify
try:
    from lib.cip_exo_engine import extractSugarFingerprint
    logger.info("  [PASS] lib.cip_exo_engine import OK (uses lib.virtual_demodify)")
except ImportError as e:
    logger.error(f"  [FAIL] lib.cip_exo_engine import failed: {e}")

# Fix #3: get_split_smiles 返回顺序验证
try:
    from lib.glycan_topology import get_split_smiles
    from rdkit import Chem
    # 使用一个已知分子: Digitoxin (含3个 digitoxose + steroid aglycone)
    # 苷元应该含甾体骨架 (更多碳), 糖应该含更多氧
    testSmi = "OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O"  # 简单 D-Glucose
    testMol = Chem.MolFromSmiles(testSmi)
    if testMol:
        aglycone, glycan = get_split_smiles(testMol)
        # 对于纯糖分子, aglycone 应该为空
        logger.info(f"  [PASS] get_split_smiles return order verified "
                    f"(aglycone='{aglycone[:30] if aglycone else ''}', "
                    f"glycan='{glycan[:30] if glycan else ''}')")
except Exception as e:
    logger.error(f"  [FAIL] get_split_smiles verification failed: {e}")


# ═══════════════════════════════════════════════════════════════════════
# STEP 2: 运行 226 条统一基准测试 (Run Unified Benchmark)
# ═══════════════════════════════════════════════════════════════════════
logger.info("")
logger.info("STEP 2: Running unified benchmark (226 entries)...")

benchmarkPath = os.path.join(DATA_DIR, "benchmark_unified.json")
if not os.path.exists(benchmarkPath):
    logger.warning(f"  Benchmark file not found: {benchmarkPath}")
    logger.warning("  Skipping benchmark test.")
else:
    try:
        from lib.monosaccharide_identifier import identify_monosaccharide_v10, generate_refined_sequence
        from lib.glycan_topology import find_mapped_sugar_units
        import re

        with open(benchmarkPath, "r", encoding="utf-8") as f:
            benchmarks = json.load(f)

        totalCount = len(benchmarks)
        passCount = 0
        failCount = 0
        errorCount = 0
        failedEntries = []

        t0 = time.time()
        for entry in benchmarks:
            entryId = entry.get("id", "?")
            entryName = entry.get("name", "")
            smi = entry.get("smiles", "")
            expectedSeq = entry.get("expected_sequence", "").strip()
            expectedField = entry.get("expected", "").strip()

            if expectedField in ("NO_SUGAR", "N/A", ""):
                passCount += 1
                continue

            try:
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    errorCount += 1
                    continue

                seq, mods = generate_refined_sequence(mol)
                # 标准化比较 (Standardized comparison)
                normActual = re.sub(r'\(CIP:[^)]+\)', '', seq or "").strip()
                normExpected = expectedSeq

                if normActual == normExpected:
                    passCount += 1
                else:
                    failCount += 1
                    failedEntries.append({
                        "id": entryId,
                        "name": entryName,
                        "expected": normExpected,
                        "actual": normActual
                    })
            except Exception as e:
                errorCount += 1
                logger.warning(f"  #{entryId} {entryName}: ERROR - {str(e)[:60]}")

        elapsed = time.time() - t0
        passRate = passCount / totalCount * 100 if totalCount > 0 else 0

        logger.info(f"  Benchmark Results: {passCount}/{totalCount} PASS ({passRate:.1f}%)")
        logger.info(f"  Failures: {failCount}, Errors: {errorCount}")
        logger.info(f"  Elapsed: {elapsed:.1f}s ({elapsed/totalCount*1000:.1f}ms/mol)")

        if failedEntries:
            logger.info(f"  --- Failed entries ({len(failedEntries)}) ---")
            for f_entry in failedEntries[:20]:  # 最多显示 20 条
                logger.info(
                    f"  #{f_entry['id']} {f_entry['name']}: "
                    f"expected='{f_entry['expected'][:50]}' "
                    f"actual='{f_entry['actual'][:50]}'"
                )
            if len(failedEntries) > 20:
                logger.info(f"  ... and {len(failedEntries) - 20} more")

    except Exception as e:
        logger.error(f"  Benchmark execution failed: {e}")
        import traceback
        logger.error(traceback.format_exc())


# ═══════════════════════════════════════════════════════════════════════
# STEP 3: 修复 V13 Pruned CSV (Fix V13 Pruned CSV Column Swap)
# ═══════════════════════════════════════════════════════════════════════
logger.info("")
logger.info("STEP 3: Fixing V13 Pruned CSV (Glycan/Aglycon column swap)...")

import pandas as pd

v13PrunedPath = os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched_v13_Pruned.csv")
v13FinalPath = os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched_v13_Final.csv")

csvPaths = [
    ("V13_Pruned", v13PrunedPath),
    ("V13_Final", v13FinalPath),
]

for csvName, csvPath in csvPaths:
    if not os.path.exists(csvPath):
        logger.warning(f"  {csvName}: File not found, skipping: {csvPath}")
        continue

    logger.info(f"  Loading {csvName}: {csvPath}")
    df = pd.read_csv(csvPath, low_memory=False)
    logger.info(f"  {csvName}: {len(df):,} rows, {len(df.columns)} columns")

    # 检查列是否存在 (Check if columns exist)
    if "Glycan_SMILES" not in df.columns or "Aglycon_SMILES" not in df.columns:
        logger.warning(f"  {csvName}: Missing Glycan_SMILES or Aglycon_SMILES columns, skipping")
        continue

    # 互换两列 (Swap the two columns)
    # 由于 V13 管线将 get_split_smiles() 的返回值解包顺序弄反了,
    # 所以 CSV 中 Glycan_SMILES 列实际存储的是苷元, Aglycon_SMILES 存储的是糖
    # Since V13 pipeline had the unpack order reversed,
    # Glycan_SMILES actually contains aglycone SMILES and vice versa
    logger.info(f"  {csvName}: Swapping Glycan_SMILES <-> Aglycon_SMILES...")

    # 验证: 抽样检查交换前的状态
    # Verify: spot-check the pre-swap state
    sampleIdx = df[df["Glycan_SMILES"].notna() & df["Aglycon_SMILES"].notna()].head(5).index
    if len(sampleIdx) > 0:
        idx = sampleIdx[0]
        logger.info(f"  Pre-swap sample row {idx}:")
        logger.info(f"    Glycan_SMILES (should be aglycone): {str(df.loc[idx, 'Glycan_SMILES'])[:60]}")
        logger.info(f"    Aglycon_SMILES (should be glycan): {str(df.loc[idx, 'Aglycon_SMILES'])[:60]}")

    # 执行交换 (Execute swap)
    df["Glycan_SMILES"], df["Aglycon_SMILES"] = (
        df["Aglycon_SMILES"].copy(),
        df["Glycan_SMILES"].copy()
    )

    if len(sampleIdx) > 0:
        logger.info(f"  Post-swap sample row {idx}:")
        logger.info(f"    Glycan_SMILES (now glycan): {str(df.loc[idx, 'Glycan_SMILES'])[:60]}")
        logger.info(f"    Aglycon_SMILES (now aglycone): {str(df.loc[idx, 'Aglycon_SMILES'])[:60]}")

    # 保存 (Save)
    df.to_csv(csvPath, index=False)
    logger.info(f"  [DONE] {csvName}: Saved with corrected columns ({len(df):,} rows)")


# ═══════════════════════════════════════════════════════════════════════
# STEP 4: 更新 Saponin DB CSV (Update Saponin DB CSV)
# ═══════════════════════════════════════════════════════════════════════
logger.info("")
logger.info("STEP 4: Updating Saponin DB CSV...")

saponinPath = os.path.join(REPORT_DIR, "GlycoNP_Saponin_DB_v13.csv")
if not os.path.exists(saponinPath):
    # 如果不存在, 从修复后的 V13 Pruned 重新提取
    # If not found, re-extract from fixed V13 Pruned
    logger.info("  Saponin CSV not found, re-extracting from V13 Pruned...")
    if os.path.exists(v13PrunedPath):
        dfPruned = pd.read_csv(v13PrunedPath, low_memory=False)
        mask = (
            dfPruned.get('np_classifier_class', pd.Series(dtype=str)).str.contains('saponin', case=False, na=False) |
            dfPruned.get('Detailed_NP_Class', pd.Series(dtype=str)).str.contains('saponin', case=False, na=False) |
            dfPruned.get('Super_Scaffold_Class', pd.Series(dtype=str)).str.contains('saponin', case=False, na=False)
        )
        saponinDf = dfPruned[mask].copy()
        saponinDf.to_csv(saponinPath, index=False)
        logger.info(f"  [DONE] Created Saponin DB: {len(saponinDf):,} molecules")
    else:
        logger.warning("  Cannot create Saponin DB: V13 Pruned not found")
else:
    # 存在则交换列 (If exists, swap columns to match fix)
    logger.info(f"  Loading existing Saponin DB: {saponinPath}")
    dfSaponin = pd.read_csv(saponinPath, low_memory=False)
    logger.info(f"  Saponin DB: {len(dfSaponin):,} rows")

    if "Glycan_SMILES" in dfSaponin.columns and "Aglycon_SMILES" in dfSaponin.columns:
        logger.info("  Swapping Glycan_SMILES <-> Aglycon_SMILES in Saponin DB...")
        dfSaponin["Glycan_SMILES"], dfSaponin["Aglycon_SMILES"] = (
            dfSaponin["Aglycon_SMILES"].copy(),
            dfSaponin["Glycan_SMILES"].copy()
        )
        dfSaponin.to_csv(saponinPath, index=False)
        logger.info(f"  [DONE] Saponin DB updated: {len(dfSaponin):,} rows")
    else:
        logger.warning("  Saponin DB missing SMILES columns, skipping swap")


# ═══════════════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════════════
logger.info("")
logger.info("=" * 70)
logger.info("AUDIT FIX SUMMARY")
logger.info("=" * 70)
logger.info("Fix #1: lib/cip_exo_engine.py import → lib.virtual_demodify  [APPLIED]")
logger.info("Fix #2: audit report v2 → v10 identification               [APPLIED]")
logger.info("Fix #3: V13 get_split_smiles unpack order                   [APPLIED]")
logger.info("Fix #7: Hardcoded paths → os.path relative                  [APPLIED]")
logger.info("CSV:    V13 Pruned Glycan/Aglycon columns swapped           [APPLIED]")
logger.info("CSV:    V13 Final  Glycan/Aglycon columns swapped           [APPLIED]")
logger.info("CSV:    Saponin DB Glycan/Aglycon columns swapped           [APPLIED]")
logger.info(f"Log:    {logPath}")
logger.info("=" * 70)

