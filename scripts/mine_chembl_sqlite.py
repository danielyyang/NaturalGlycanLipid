"""
GlycoNP × ChEMBL SQLite 直接查询引擎
GlycoNP × ChEMBL SQLite Direct Query Engine

跳过 CSV 中间层, 直接对本地 30GB ChEMBL SQLite 数据库执行联表查询。
Skip CSV intermediary — query local 30GB ChEMBL SQLite database directly.

核心 SQL: compound_structures → molecule_dictionary → activities → assays → target_dictionary
5-table JOIN: InChIKey → molregno → activities → assays → targets

输出: 为每个 GlycoNP 化合物聚合 bioactivity_summary 字段。
Output: bioactivity_summary field aggregated per GlycoNP compound.

使用方法 / Usage:
  python scripts/mine_chembl_sqlite.py
  python scripts/mine_chembl_sqlite.py --pilot 1000  # 测试前 1000 个化合物
"""
import argparse
import os
import sqlite3
import sys
import time
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


# =====================================================================
# 1. InChIKey 碰撞: 确定 GlycoNP ∩ ChEMBL 的化合物集合
# InChIKey collision: determine GlycoNP ∩ ChEMBL compound set
# =====================================================================

def collideInchiKeys(
    glycoKeys: Set[str],
    conn: sqlite3.Connection,
) -> Dict[str, int]:
    """
    从 ChEMBL SQLite 的 compound_structures 表中查找与 GlycoNP InChIKey
    碰撞的分子, 返回 InChIKey → molregno 映射。

    Find GlycoNP InChIKey matches in ChEMBL compound_structures table.
    Returns InChIKey → molregno mapping.

    设计意图: 先做碰撞筛选, 再对命中分子执行昂贵的活性查询。
    这样避免了在 24M 行的 activities 表上全量扫描。
    Design: collision first, then expensive activity queries only for hits.
    """
    print("  [1/3] InChIKey collision with ChEMBL...")
    t0 = time.time()

    # 分批查询以避免 SQLite "too many SQL variables" 限制 (max ~999)
    # Batch query to avoid SQLite variable limit
    BATCH_SIZE = 500
    keyList = sorted(glycoKeys)
    keyToMolregno: Dict[str, int] = {}

    for i in range(0, len(keyList), BATCH_SIZE):
        batch = keyList[i:i + BATCH_SIZE]
        placeholders = ",".join(["?" for _ in batch])
        query = f"""
        SELECT standard_inchi_key, molregno
        FROM compound_structures
        WHERE standard_inchi_key IN ({placeholders})
        """
        cur = conn.execute(query, batch)
        for row in cur.fetchall():
            keyToMolregno[row[0]] = row[1]

    elapsed = time.time() - t0
    hitRate = len(keyToMolregno) / len(glycoKeys) * 100 if glycoKeys else 0
    print(f"    GlycoNP keys: {len(glycoKeys):,}")
    print(f"    ChEMBL hits:  {len(keyToMolregno):,} ({hitRate:.1f}%)")
    print(f"    Time: {elapsed:.1f}s")

    return keyToMolregno


# =====================================================================
# 2. 活性数据批量提取: 5 表联查
# Batch activity extraction: 5-table JOIN
# =====================================================================

def extractActivitiesForHits(
    hitMolregnos: Set[int],
    conn: sqlite3.Connection,
) -> pd.DataFrame:
    """
    对命中的 molregno 集合, 执行 5 表联查提取活性数据。
    For hit molregnos, execute 5-table JOIN to extract activity data.

    SQL 联查逻辑:
      activities (molregno, assay_id, standard_type, pchembl_value)
      → assays (assay_id → tid)
      → target_dictionary (tid → pref_name, target_type)

    过滤条件 / Filters:
      - standard_type IN ('IC50', 'Ki', 'EC50', 'Kd', 'MIC') — 可量化指标
      - pchembl_value IS NOT NULL — 保证可比较性

    Args:
        hitMolregnos: 命中分子的 molregno 集合
        conn: ChEMBL SQLite 连接

    Returns:
        DataFrame with columns: molregno, target, standard_type, pchembl_value
    """
    print("  [2/3] Extracting activity data (5-table JOIN)...")
    t0 = time.time()

    # 分批查询 (Batch query)
    BATCH_SIZE = 500
    molregnoList = sorted(hitMolregnos)
    allRows: List[dict] = []

    for i in range(0, len(molregnoList), BATCH_SIZE):
        batch = molregnoList[i:i + BATCH_SIZE]
        placeholders = ",".join(["?" for _ in batch])

        # 核心 SQL: 5 表联查
        # Core SQL: 5-table JOIN
        query = f"""
        SELECT
            act.molregno,
            td.pref_name    AS target,
            td.target_type,
            act.standard_type,
            act.standard_value,
            act.standard_units,
            act.pchembl_value,
            ass.assay_type
        FROM activities act
        JOIN assays ass          ON act.assay_id = ass.assay_id
        JOIN target_dictionary td ON ass.tid = td.tid
        WHERE act.molregno IN ({placeholders})
          AND act.standard_type IN ('IC50', 'Ki', 'EC50', 'Kd', 'MIC')
          AND act.pchembl_value IS NOT NULL
          AND td.pref_name IS NOT NULL
        """

        cur = conn.execute(query, batch)
        columns = [desc[0] for desc in cur.description]
        for row in cur.fetchall():
            allRows.append(dict(zip(columns, row)))

        if (i // BATCH_SIZE + 1) % 10 == 0:
            print(f"    ... queried {i + len(batch):,}/{len(molregnoList):,} "
                  f"molregnos, {len(allRows):,} activity rows")

    actDf = pd.DataFrame(allRows)
    elapsed = time.time() - t0

    if not actDf.empty:
        print(f"    Activity rows: {len(actDf):,}")
        print(f"    Unique targets: {actDf['target'].nunique():,}")
        print(f"    Activity types: {actDf['standard_type'].value_counts().to_dict()}")
    else:
        print(f"    No activity data found for hit molecules")

    print(f"    Time: {elapsed:.1f}s")
    return actDf


# =====================================================================
# 3. 聚合为 bioactivity_summary 字符串
# Aggregate into bioactivity_summary string
# =====================================================================

def aggregateBioactivity(
    df: pd.DataFrame,
    keyToMolregno: Dict[str, int],
    actDf: pd.DataFrame,
    maxTargetsPerCompound: int = 5,
) -> pd.DataFrame:
    """
    将活性数据聚合为每个化合物的 bioactivity_summary 字段。
    Aggregate activity data into a bioactivity_summary field per compound.

    格式 / Format:
      "COX-2:IC50(pC=5.2) | Alpha-glucosidase:Ki(pC=6.8) | ..."

    同时生成详细列:
      - ChEMBL_Targets: "COX-2 | Alpha-glucosidase | ..."
      - ChEMBL_Activity_Types: "IC50, Ki"
      - ChEMBL_Best_pChEMBL: 6.8 (最大 pChEMBL 值)
      - ChEMBL_N_Activities: 活性数据条数
    """
    print("  [3/3] Aggregating bioactivity summary...")
    t0 = time.time()

    # 反转映射: molregno → inchikey
    # Reverse mapping: molregno → inchikey
    molregnoToKey: Dict[int, str] = {v: k for k, v in keyToMolregno.items()}

    # 初始化列 (Initialize columns)
    df["bioactivity_summary"] = ""
    df["ChEMBL_Targets"] = ""
    df["ChEMBL_Activity_Types"] = ""
    df["ChEMBL_Best_pChEMBL"] = np.nan
    df["ChEMBL_N_Activities"] = 0

    if actDf.empty:
        print("    No activity data to aggregate")
        return df

    # 按 molregno 分组 (Group by molregno)
    grouped = actDf.groupby("molregno")

    filled = 0
    for idx in df.index:
        ik = str(df.at[idx, "standard_inchi_key"])
        if ik not in keyToMolregno:
            continue

        molregno = keyToMolregno[ik]
        if molregno not in grouped.groups:
            continue

        compAct = grouped.get_group(molregno)

        # 靶点去重聚合 (Unique targets with best pChEMBL)
        targetSummaries = []
        targetPchembl = {}
        for _, actRow in compAct.iterrows():
            target = str(actRow["target"])
            stdType = str(actRow["standard_type"])
            pchembl = float(actRow["pchembl_value"]) if actRow["pchembl_value"] else 0

            if target not in targetPchembl or pchembl > targetPchembl[target][1]:
                targetPchembl[target] = (stdType, pchembl)

        # 按 pChEMBL 降序排序 (Sort by pChEMBL descending)
        sortedTargets = sorted(targetPchembl.items(),
                               key=lambda x: x[1][1], reverse=True)

        for target, (stdType, pchembl) in sortedTargets[:maxTargetsPerCompound]:
            targetSummaries.append(f"{target}:{stdType}(pC={pchembl:.1f})")

        # 填充列 (Fill columns)
        df.at[idx, "bioactivity_summary"] = " | ".join(targetSummaries)
        df.at[idx, "ChEMBL_Targets"] = " | ".join(
            [t for t, _ in sortedTargets[:maxTargetsPerCompound]])
        df.at[idx, "ChEMBL_Activity_Types"] = ", ".join(
            sorted(set(str(r["standard_type"]) for _, r in compAct.iterrows())))
        df.at[idx, "ChEMBL_Best_pChEMBL"] = compAct["pchembl_value"].astype(
            float).max()
        df.at[idx, "ChEMBL_N_Activities"] = len(compAct)

        filled += 1

    elapsed = time.time() - t0
    print(f"    Filled {filled:,} compounds with bioactivity data")
    print(f"    Time: {elapsed:.1f}s")

    return df


# =====================================================================
# 4. Top-20 靶点报告 (Top-20 Targets Report)
# =====================================================================

def printTargetReport(df: pd.DataFrame):
    """打印 Top-20 靶点汇总 / Print top-20 target summary."""
    hasBioact = df["ChEMBL_Targets"].notna() & (df["ChEMBL_Targets"] != "")
    nWithBioact = hasBioact.sum()
    print(f"\n  === Bioactivity Summary ===")
    print(f"  Compounds with activity: {nWithBioact:,} / {len(df):,} "
          f"({nWithBioact/len(df)*100:.1f}%)")

    if nWithBioact == 0:
        return

    # 统计靶点 (Count targets)
    from collections import Counter
    targetCounter: Counter = Counter()
    for targets in df.loc[hasBioact, "ChEMBL_Targets"]:
        for t in str(targets).split(" | "):
            t = t.strip()
            if t and t != "nan":
                targetCounter[t] += 1

    print(f"\n  Top 20 Targets:")
    print(f"  {'Rank':>4}  {'Target':<50} {'Count':>6}")
    print(f"  {'----':>4}  {'-----':<50} {'-----':>6}")
    for i, (target, count) in enumerate(targetCounter.most_common(20), 1):
        print(f"  {i:4d}  {target:<50} {count:6d}")

    # pChEMBL 分布 (pChEMBL distribution)
    pvals = pd.to_numeric(df["ChEMBL_Best_pChEMBL"], errors="coerce").dropna()
    if not pvals.empty:
        print(f"\n  pChEMBL Distribution:")
        print(f"    Mean: {pvals.mean():.2f}")
        print(f"    Median: {pvals.median():.2f}")
        print(f"    Max: {pvals.max():.2f}")
        print(f"    >7 (IC50 < 100nM): {(pvals >= 7).sum():,}")
        print(f"    >6 (IC50 < 1uM):   {(pvals >= 6).sum():,}")


# =====================================================================
# Main
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        description="GlycoNP x ChEMBL SQLite Direct Mining")
    parser.add_argument("--input", type=str, default=None,
                        help="GlycoNP pipeline CSV path")
    parser.add_argument("--db", type=str, default=None,
                        help="ChEMBL SQLite database path")
    parser.add_argument("--output", type=str, default=None,
                        help="Output enriched CSV path")
    parser.add_argument("--pilot", type=int, default=0,
                        help="Only process first N compounds (0 = all)")
    args = parser.parse_args()

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")

    # 解析路径 (Resolve paths)
    inputPath = args.input or os.path.join(
        reportDir, "GlycoNP_Pipeline_Full_Cleaned.csv")
    dbPath = args.db or os.path.join(
        baseDir, "data", "chembl_36", "chembl_36_sqlite", "chembl_36.db")
    outputPath = args.output or os.path.join(
        reportDir, "GlycoNP_ChEMBL_Enriched.csv")

    print("=" * 70)
    print("  GlycoNP x ChEMBL SQLite Direct Mining Engine")
    print("  糖缀合物 × ChEMBL 本地 SQLite 直接查询引擎")
    print("=" * 70)
    print(f"  Input:    {inputPath}")
    print(f"  ChEMBL:   {dbPath}")
    print(f"  Output:   {outputPath}")

    # 验证文件存在 (Verify files exist)
    if not os.path.exists(inputPath):
        print(f"\n  [ERROR] Input not found: {inputPath}")
        return
    if not os.path.exists(dbPath):
        print(f"\n  [ERROR] ChEMBL SQLite not found: {dbPath}")
        print(f"  Download: https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/")
        return

    t0 = time.time()

    # Step 0: 加载 GlycoNP 数据 (Load GlycoNP data)
    df = pd.read_csv(inputPath, low_memory=False, dtype=str,
                     encoding="utf-8-sig")
    print(f"\n  Loaded: {len(df):,} rows")

    if args.pilot > 0:
        df = df.head(args.pilot).copy()
        print(f"  Pilot mode: {len(df):,} rows")

    # 提取 InChIKey 集合 (Extract InChIKey set)
    glycoKeys = set(df["standard_inchi_key"].dropna().astype(str))
    glycoKeys.discard("nan")
    glycoKeys.discard("")
    print(f"  Valid InChIKeys: {len(glycoKeys):,}")

    # Step 1-3: 连接 SQLite 并执行查询 (Connect and query)
    conn = sqlite3.connect(dbPath)
    try:
        keyToMolregno = collideInchiKeys(glycoKeys, conn)

        if not keyToMolregno:
            print("\n  [WARN] No InChIKey hits! Nothing to enrich.")
            conn.close()
            return

        hitMolregnos = set(keyToMolregno.values())
        actDf = extractActivitiesForHits(hitMolregnos, conn)
    finally:
        conn.close()

    df = aggregateBioactivity(df, keyToMolregno, actDf)

    # 报告 (Report)
    printTargetReport(df)

    # 保存 (Save)
    df.to_csv(outputPath, index=False, encoding="utf-8-sig")
    elapsed = time.time() - t0
    print(f"\n  Output: {outputPath}")
    print(f"  Total time: {elapsed:.0f}s")
    print("=" * 70)


if __name__ == "__main__":
    main()
