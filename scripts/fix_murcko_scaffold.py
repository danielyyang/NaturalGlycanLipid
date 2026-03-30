"""
fix_murcko_scaffold.py
======================
修复 V13 CSV 中的 Murcko_Scaffold 列。
Fix Murcko_Scaffold column in V13 CSVs.

原因: V13 管线的 swap Bug 导致 Murcko_Scaffold 基于糖碎片(而非苷元)生成。
Reason: V13 pipeline's swap bug caused Murcko_Scaffold to be computed from
        glycan fragments instead of aglycone fragments.

修复: 用交换后的 Aglycon_SMILES (现在是真正的苷元) 重新计算拓扑骨架。
Fix: Recompute topology scaffold from the corrected Aglycon_SMILES column.

[TEST DATA ONLY]
"""
import os
import sys
import time
import logging

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import pandas as pd
from rdkit import Chem
from lib.feature_extractor import getTopologyScaffoldSmiles

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
REPORT_DIR = os.path.join(BASE_DIR, "reports")
LOG_DIR = os.path.join(BASE_DIR, "log")

logPath = os.path.join(LOG_DIR, "fix_murcko_2026-03-30.log")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(logPath, encoding="utf-8"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("FixMurcko")


def recomputeScaffold(smi: str) -> str:
    """从苷元 SMILES 重新计算拓扑骨架。
    Recompute topology scaffold from aglycone SMILES.
    """
    if pd.isna(smi) or not smi or str(smi) == "nan" or len(str(smi)) < 4:
        return ""
    try:
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            return ""
        return getTopologyScaffoldSmiles(mol) or ""
    except Exception:
        return ""


def fixFile(csvPath: str, name: str):
    """修复单个 CSV 文件的 Murcko_Scaffold 列。"""
    if not os.path.exists(csvPath):
        logger.warning(f"  {name}: File not found, skipping")
        return

    logger.info(f"Loading {name}: {csvPath}")
    df = pd.read_csv(csvPath, low_memory=False)
    logger.info(f"  {name}: {len(df):,} rows, {len(df.columns)} columns")

    if "Aglycon_SMILES" not in df.columns:
        logger.warning(f"  {name}: Missing Aglycon_SMILES column, skipping")
        return

    if "Murcko_Scaffold" not in df.columns:
        logger.warning(f"  {name}: Missing Murcko_Scaffold column, will create it")

    # 保存旧值用于对比 (Save old values for comparison)
    oldScaffold = df.get("Murcko_Scaffold", pd.Series(dtype=str)).copy()

    # 重新计算 (Recompute)
    logger.info(f"  Recomputing Murcko_Scaffold from Aglycon_SMILES...")
    t0 = time.time()
    df["Murcko_Scaffold"] = df["Aglycon_SMILES"].apply(recomputeScaffold)
    elapsed = time.time() - t0
    logger.info(f"  Recomputed in {elapsed:.1f}s")

    # 统计变化 (Statistics)
    nonEmpty = (df["Murcko_Scaffold"] != "").sum()
    logger.info(f"  Non-empty scaffolds: {nonEmpty:,}/{len(df):,} ({nonEmpty/len(df)*100:.1f}%)")

    if len(oldScaffold) == len(df):
        changedMask = (oldScaffold.fillna("") != df["Murcko_Scaffold"].fillna(""))
        changedCount = changedMask.sum()
        logger.info(f"  Changed scaffolds: {changedCount:,}/{len(df):,} ({changedCount/len(df)*100:.1f}%)")

        # 显示前 5 个变化样本 (Show first 5 changed samples)
        changedIdx = df[changedMask].head(5).index
        for idx in changedIdx:
            logger.info(f"    Row {idx}: OLD='{str(oldScaffold.iloc[idx])[:40]}' -> "
                        f"NEW='{str(df.loc[idx, 'Murcko_Scaffold'])[:40]}'")

    # 保存 (Save)
    df.to_csv(csvPath, index=False)
    logger.info(f"  [DONE] {name}: Saved ({len(df):,} rows)")


# ─── 主程序 ───
logger.info("=" * 70)
logger.info("Fix Murcko_Scaffold — Recompute from corrected Aglycon_SMILES")
logger.info("=" * 70)

fixFile(
    os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched_v13_Pruned.csv"),
    "V13_Pruned"
)
fixFile(
    os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched_v13_Final.csv"),
    "V13_Final"
)
fixFile(
    os.path.join(REPORT_DIR, "GlycoNP_Saponin_DB_v13.csv"),
    "Saponin_DB"
)

logger.info("")
logger.info("=" * 70)
logger.info("ALL DONE — Murcko_Scaffold recomputed from true aglycone SMILES")
logger.info(f"Log: {logPath}")
logger.info("=" * 70)
