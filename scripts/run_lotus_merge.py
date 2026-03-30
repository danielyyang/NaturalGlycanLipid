"""
Task 2 实质性合并: LOTUS 离线数据 → GlycoNP 回填
Task 2 Merge: LOTUS offline data → GlycoNP backfill

使用 non_isomeric_inchikey_block1 作为二维万能钥匙进行跨库匹配。
Uses non_isomeric_inchikey_block1 as 2D universal key for cross-DB matching.

注意: 不覆盖原始 standard_inchi_key (3D 空间信息有价值)。
NOTE: Original standard_inchi_key is NEVER overwritten (3D info is valuable).
"""
import os
import sys
import time

import pandas as pd
import numpy as np

# ---- 配置 (Configuration) ----
BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
REPORT_DIR = os.path.join(BASE_DIR, "reports")

IMPUTED_CSV = os.path.join(REPORT_DIR, "GlycoNP_Imputed.csv")
LOTUS_GZ = os.path.join(BASE_DIR, "data", "230106_frozen_metadata.csv.gz")
OUTPUT_CSV = os.path.join(REPORT_DIR, "GlycoNP_Imputed_Merged.csv")

# LOTUS 列映射 (LOTUS column mapping)
LOTUS_KEY_COL = "structure_inchikey"
LOTUS_FIELDS = {
    "reference_doi": "dois",
    "organism_name": "organisms",
    "structure_nameTraditional": "name",
    "organism_taxonomy_02kingdom": "LOTUS_kingdom",
    "organism_taxonomy_03phylum": "LOTUS_phylum",
    "organism_taxonomy_06family": "LOTUS_family",
}


def main():
    print("=" * 70)
    print("  Task 2: LOTUS Merge & Backfill")
    print("  使用二维钥匙进行离线知识库回填")
    print("=" * 70)

    # ---- Step 0: 验证文件 (Verify files) ----
    if not os.path.exists(IMPUTED_CSV):
        print(f"  [ERROR] Imputed CSV not found: {IMPUTED_CSV}")
        return
    if not os.path.exists(LOTUS_GZ):
        print(f"  [ERROR] LOTUS file not found: {LOTUS_GZ}")
        return

    t0 = time.time()

    # ---- Step 1: 加载 GlycoNP (Load GlycoNP) ----
    print(f"\n  [1/5] Loading GlycoNP imputed data...")
    df = pd.read_csv(IMPUTED_CSV, dtype=str, low_memory=False,
                     encoding="utf-8-sig")
    total = len(df)
    print(f"    Rows: {total:,}")

    # 记录补全前统计 (Record pre-merge stats)
    def countValid(series):
        return series.notna().sum() - series.astype(str).isin(
            ["", "nan", "None"]).sum()

    beforeStats = {
        "dois": countValid(df["dois"]),
        "organisms": countValid(df["organisms"]),
        "name": countValid(df["name"]),
    }

    print(f"    Before merge:")
    for field, count in beforeStats.items():
        print(f"      {field}: {count:,} / {total:,} "
              f"({count/total*100:.1f}% valid)")

    # ---- Step 2: 加载 LOTUS (Load LOTUS) ----
    print(f"\n  [2/5] Loading LOTUS: {LOTUS_GZ}")
    useCols = [LOTUS_KEY_COL] + list(LOTUS_FIELDS.keys())
    lotusDf = pd.read_csv(LOTUS_GZ, dtype=str, low_memory=False,
                          usecols=useCols)
    print(f"    LOTUS rows: {len(lotusDf):,}")
    print(f"    Unique InChIKeys: {lotusDf[LOTUS_KEY_COL].nunique():,}")

    # ---- Step 3: 生成 LOTUS Block-1 (Generate LOTUS Block-1) ----
    print(f"\n  [3/5] Generating LOTUS 2D Block-1 keys...")
    lotusDf["_lotus_block1"] = (
        lotusDf[LOTUS_KEY_COL].astype(str).str[:14]
    )
    lotusDf = lotusDf.dropna(subset=["_lotus_block1"])
    lotusDf = lotusDf[lotusDf["_lotus_block1"] != "nan"]

    # 聚合: 同一个 Block-1 可能有多条 LOTUS 记录
    # Aggregate: same Block-1 may have multiple LOTUS records
    # 策略: DOI 拼接, organism/name 取第一个非空值
    # Strategy: concatenate DOIs, take first non-null for organism/name
    print(f"    Aggregating by Block-1...")

    def aggFunc(group):
        result = {}
        # DOI: 收集所有唯一 DOI, 用 " | " 拼接
        dois = group["reference_doi"].dropna().unique()
        dois = [d for d in dois if d and d != "nan"]
        result["_lotus_doi"] = " | ".join(dois[:5]) if dois else ""

        # Organism: 取第一个非空
        for col in ["organism_name"]:
            vals = group[col].dropna()
            vals = vals[vals.astype(str) != "nan"]
            result[f"_lotus_{col}"] = vals.iloc[0] if len(vals) > 0 else ""

        # Name: 取第一个非空
        for col in ["structure_nameTraditional"]:
            vals = group[col].dropna()
            vals = vals[vals.astype(str) != "nan"]
            result[f"_lotus_{col}"] = vals.iloc[0] if len(vals) > 0 else ""

        # Taxonomy: 取第一个非空
        for col in ["organism_taxonomy_02kingdom",
                     "organism_taxonomy_03phylum",
                     "organism_taxonomy_06family"]:
            vals = group[col].dropna()
            vals = vals[vals.astype(str) != "nan"]
            result[f"_lotus_{col}"] = vals.iloc[0] if len(vals) > 0 else ""

        return pd.Series(result)

    lotusAgg = lotusDf.groupby("_lotus_block1").apply(aggFunc).reset_index()
    print(f"    Aggregated LOTUS entries: {len(lotusAgg):,}")

    # ---- Step 4: Merge (Left join) ----
    print(f"\n  [4/5] Merging on non_isomeric_inchikey_block1...")
    merged = df.merge(
        lotusAgg,
        left_on="non_isomeric_inchikey_block1",
        right_on="_lotus_block1",
        how="left",
    )
    hitCount = merged["_lotus_block1"].notna().sum()
    print(f"    Block-1 hits: {hitCount:,} / {total:,} ({hitCount/total*100:.1f}%)")

    # ---- Step 5: Backfill (只填空值, 不覆盖) ----
    print(f"\n  [5/5] Backfilling empty fields...")

    def backfill(targetCol, sourceCol, label):
        isBlank = (
            merged[targetCol].isna() |
            merged[targetCol].astype(str).isin(["", "nan", "None"])
        )
        hasFill = (
            merged[sourceCol].notna() &
            (~merged[sourceCol].astype(str).isin(["", "nan", "None"]))
        )
        fillMask = isBlank & hasFill
        nFilled = fillMask.sum()
        merged.loc[fillMask, targetCol] = merged.loc[fillMask, sourceCol]
        print(f"    {label}: +{nFilled:,} rows backfilled")
        return nFilled

    doiFilled = backfill("dois", "_lotus_doi", "DOI")
    orgFilled = backfill("organisms", "_lotus_organism_name", "Organisms")
    nameFilled = backfill("name", "_lotus_structure_nameTraditional", "Name")

    # 新增 taxonomy 列 (Add taxonomy columns — these are always new)
    for col, newName in [
        ("_lotus_organism_taxonomy_02kingdom", "LOTUS_kingdom"),
        ("_lotus_organism_taxonomy_03phylum", "LOTUS_phylum"),
        ("_lotus_organism_taxonomy_06family", "LOTUS_family"),
    ]:
        if col in merged.columns:
            merged[newName] = merged[col]

    # 清理临时列 (Clean temp columns)
    dropCols = [c for c in merged.columns if c.startswith("_lotus_")]
    merged.drop(columns=dropCols, inplace=True, errors="ignore")

    # ---- 最终统计 (Final statistics) ----
    afterStats = {
        "dois": countValid(merged["dois"]),
        "organisms": countValid(merged["organisms"]),
        "name": countValid(merged["name"]),
    }

    print("\n" + "=" * 70)
    print("  数据缺口分析对比报告 / Gap Analysis Before vs After")
    print("=" * 70)
    print(f"\n  {'Field':<12} {'Before':>10} {'After':>10} {'Gained':>10} "
          f"{'Before%':>9} {'After%':>9}")
    print(f"  {'-'*12} {'-'*10} {'-'*10} {'-'*10} {'-'*9} {'-'*9}")

    for field in ["dois", "organisms", "name"]:
        before = beforeStats[field]
        after = afterStats[field]
        gained = after - before
        beforePct = before / total * 100
        afterPct = after / total * 100
        print(f"  {field:<12} {before:>10,} {after:>10,} {'+' + str(gained) if gained >= 0 else str(gained):>10} "
              f"{beforePct:>8.1f}% {afterPct:>8.1f}%")

    print(f"\n  总行数 (Total rows): {total:,}")
    print(f"  Block-1 命中: {hitCount:,} ({hitCount/total*100:.1f}%)")
    print(f"  DOI 抢救: +{doiFilled:,} 条文献来源")
    print(f"  Organisms 抢救: +{orgFilled:,} 条物种分类")
    print(f"  Name 抢救: +{nameFilled:,} 条化合物名称")

    # 检查新增的 LOTUS taxonomy 覆盖率
    for col in ["LOTUS_kingdom", "LOTUS_phylum", "LOTUS_family"]:
        if col in merged.columns:
            validCount = countValid(merged[col])
            print(f"  {col}: {validCount:,} ({validCount/total*100:.1f}%)")

    # ---- 保存 (Save) ----
    merged.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")
    elapsed = time.time() - t0
    print(f"\n  Output: {OUTPUT_CSV}")
    print(f"  Time: {elapsed:.0f}s")
    print("=" * 70)


if __name__ == "__main__":
    main()
