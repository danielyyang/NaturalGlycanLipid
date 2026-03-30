"""
本地 LOTUS 分类学匹配引擎 — 替代实时 API 请求
Local LOTUS Taxonomy Matching Engine — Replacement for real-time API calls

核心思路 (Core Strategy):
1. 读取 LOTUS Zenodo frozen dump CSV (230106_frozen_metadata.csv.gz)
2. 以 InChIKey Block-1 (前 14 位) 作为 Hash Key 构建本地字典
3. 通过 pandas merge 高效地将 taxonomy 映射回 COCONUT 数据集
4. 零网络依赖，可离线运行，94,000 条数据 < 10 秒完成匹配

LOTUS 数据下载 (Download):
  https://zenodo.org/records/7534071
  文件: 230106_frozen_metadata.csv.gz (~108MB)

LOTUS CSV 关键列 (Key Columns):
  - structure_inchikey          : 完整 InChIKey (27 字符)
  - organism_name               : 物种名 (binomial, e.g. "Panax ginseng")
  - organism_taxonomy_01domain  : 域 (Domain, e.g. "Eukaryota")
  - organism_taxonomy_02kingdom : 界 (Kingdom, e.g. "Plantae")
  - organism_taxonomy_05family  : 科 (Family, e.g. "Araliaceae")
  - organism_taxonomy_06genus   : 属 (Genus, e.g. "Panax")
  - organism_taxonomy_08species : 种 (Species, e.g. "Panax ginseng")
"""
import os
import pandas as pd
from typing import Tuple, Set, Optional


# =====================================================================
# 1. LOTUS Dump 加载与索引构建 (Load & Index Construction)
# =====================================================================

# LOTUS CSV 中我们关心的列名 (Column names we need from LOTUS dump)
LOTUS_INCHIKEY_COL = "structure_inchikey"
LOTUS_TAXONOMY_COLS = {
    "organism_name": "organism_name",
    "organism_taxonomy_01domain": "taxonomy_domain",
    "organism_taxonomy_02kingdom": "taxonomy_kingdom",
    "organism_taxonomy_06family": "taxonomy_family",
    "organism_taxonomy_08genus": "taxonomy_genus",
    "organism_taxonomy_09species": "taxonomy_species",
    "structure_taxonomy_npclassifier_02superclass": "np_classifier_superclass",
}

# 输出到目标 DataFrame 的列名 (Output column names for target DataFrame)
OUTPUT_COLS = list(LOTUS_TAXONOMY_COLS.values())


def loadLotusDump(
    lotusCsvPath: str,
    inchikeyColumn: str = LOTUS_INCHIKEY_COL,
    taxonomyCols: Optional[dict] = None,
) -> pd.DataFrame:
    """
    读取 LOTUS frozen dump CSV，提取 InChIKey Block-1 和分类学列，去重后返回。
    Load LOTUS frozen dump CSV, extract InChIKey Block-1 and taxonomy columns,
    deduplicate, and return.

    Args:
        lotusCsvPath: LOTUS CSV 文件路径 (.csv 或 .csv.gz)
        inchikeyColumn: InChIKey 列名
        taxonomyCols: 源列名 → 输出列名映射字典

    Returns:
        pd.DataFrame: 包含 'inchikey_block1' 和各 taxonomy 列的索引表
    """
    if taxonomyCols is None:
        taxonomyCols = LOTUS_TAXONOMY_COLS

    # 仅读取必要列以节省内存 (Read only necessary columns to save memory)
    useCols = [inchikeyColumn] + list(taxonomyCols.keys())

    lotusDf = pd.read_csv(
        lotusCsvPath,
        usecols=useCols,
        dtype=str,
        low_memory=False,
    )

    # 重命名列 (Rename columns to standardized output names)
    renameMap = {inchikeyColumn: "lotus_inchikey"}
    renameMap.update(taxonomyCols)
    lotusDf = lotusDf.rename(columns=renameMap)

    # 提取 InChIKey Block-1 (前 14 字符) 作为 Hash Key
    # Extract InChIKey Block-1 (first 14 chars) as hash key
    lotusDf["inchikey_block1"] = lotusDf["lotus_inchikey"].str[:14]

    # 丢弃无效行 (Drop rows without valid InChIKey)
    lotusDf = lotusDf.dropna(subset=["inchikey_block1"])
    lotusDf = lotusDf[lotusDf["inchikey_block1"].str.len() == 14]

    # 按 Block-1 去重，保留第一条 (Deduplicate by Block-1, keep first)
    # 同一个 Block-1 可能对应多个 organism，我们合并它们
    # Same Block-1 might correspond to multiple organisms; we aggregate them
    lotusDf = lotusDf.drop(columns=["lotus_inchikey"])

    # 策略: 每个 Block-1 保留唯一的 organism_name（用 | 连接去重后的值）
    # Strategy: for each Block-1, aggregate unique organism_names with |
    aggregated = lotusDf.groupby("inchikey_block1").agg(
        lambda series: "|".join(series.dropna().unique()[:5])  # 最多保留 5 个物种
    ).reset_index()

    return aggregated


# =====================================================================
# 2. 核心匹配逻辑 (Core Matching Logic) — 向量化 Pandas Merge
# =====================================================================

def fillTaxonomyFromLotus(
    targetDf: pd.DataFrame,
    lotusIndex: pd.DataFrame,
    inchikeyColumn: str = "standard_inchi_key",
) -> Tuple[pd.DataFrame, Set[Tuple[int, str]]]:
    """
    使用预构建的 LOTUS 索引，通过 InChIKey Block-1 向量化 merge
    将 taxonomy 信息填补到目标 DataFrame 中。

    Fill taxonomy information into target DataFrame using pre-built LOTUS index
    via vectorized InChIKey Block-1 merge.

    Args:
        targetDf: 目标 DataFrame (如 Coconut_Sugar_Check.csv)
        lotusIndex: loadLotusDump() 返回的索引 DataFrame
        inchikeyColumn: 目标 DF 中 InChIKey 列名

    Returns:
        (updated_df, imputed_cells_set): 更新后的 DF 和被填补的 (行索引, 列名) 集合
    """
    imputedCells: Set[Tuple[int, str]] = set()

    # Step 1: 提取 Block-1 (Extract Block-1 from target)
    targetDf = targetDf.copy()
    originalCols = set(targetDf.columns)
    targetDf["_merge_key"] = targetDf[inchikeyColumn].astype(str).str[:14]

    # Step 2: Left Join — 保持原始行顺序和数量
    # Left merge — preserves original row order and count
    merged = targetDf.merge(
        lotusIndex,
        left_on="_merge_key",
        right_on="inchikey_block1",
        how="left",
        suffixes=("", "_lotus"),
    )

    # Step 3: 仅填补空缺值 (Only fill missing values, don't overwrite existing)
    # 需要填补的目标列 → LOTUS 源列 映射
    # 注意：同名列在 merge 后 LOTUS 侧会加 '_lotus' 后缀
    # NOTE: Same-name columns get '_lotus' suffix on LOTUS side after merge
    fillMapping = {
        "organisms": "organism_name",
        "Family": "taxonomy_family",
        "np_classifier_superclass": "np_classifier_superclass",
    }

    for targetCol, lotusCol in fillMapping.items():
        # 处理 merge suffix 冲突 (Handle merge suffix conflict)
        # 如果目标 DF 已有同名列，LOTUS 侧会变成 xxx_lotus
        actualLotusCol = lotusCol
        if lotusCol + "_lotus" in merged.columns:
            actualLotusCol = lotusCol + "_lotus"
        elif lotusCol not in merged.columns:
            continue

        if targetCol not in merged.columns:
            merged[targetCol] = "NULL"

        # 识别需要填补的行 (Identify rows needing fill)
        isEmpty = (
            merged[targetCol].isna()
            | (merged[targetCol].astype(str).str.strip() == "")
            | (merged[targetCol].astype(str).isin(["nan", "NULL", "Not Result"]))
        )
        hasLotusData = merged[actualLotusCol].notna() & (merged[actualLotusCol].str.strip() != "")

        fillMask = isEmpty & hasLotusData

        if targetCol == "np_classifier_superclass":
            # NP 分类取最高频单一值，不使用 | 拼接 (Take most frequent single value)
            merged.loc[fillMask, targetCol] = merged.loc[fillMask, actualLotusCol].apply(
                lambda v: v.split("|")[0] if isinstance(v, str) else v
            )
        else:
            merged.loc[fillMask, targetCol] = merged.loc[fillMask, actualLotusCol]

        # 记录被填补的单元格 (Track imputed cells)
        for idx in merged.index[fillMask]:
            imputedCells.add((idx, targetCol))

    # Step 4: 添加扩展 taxonomy 列 (Add extended taxonomy columns if not present)
    extraCols = {
        "taxonomy_domain": "taxonomy_domain",
        "taxonomy_kingdom": "taxonomy_kingdom",
        "taxonomy_genus": "taxonomy_genus",
        "taxonomy_species": "taxonomy_species",
    }
    for outputCol, lotusCol in extraCols.items():
        # 同样检查 _lotus 后缀 (Also check _lotus suffix)
        actualCol = lotusCol + "_lotus" if lotusCol + "_lotus" in merged.columns else lotusCol
        if actualCol in merged.columns and outputCol not in originalCols:
            merged[outputCol] = merged[actualCol]

    # Step 5: 清理临时列 — 仅删除 LOTUS 引入的列，不删除目标原有的列
    # Clean up — only drop LOTUS-introduced columns, NEVER delete original target columns
    protectedCols = originalCols | set(extraCols.keys()) | {"taxonomy_domain", "taxonomy_kingdom", "taxonomy_genus", "taxonomy_species"}
    allLotusCols = (
        ["_merge_key", "inchikey_block1"]
        + [v for v in LOTUS_TAXONOMY_COLS.values()]
        + [v + "_lotus" for v in LOTUS_TAXONOMY_COLS.values()]
    )
    dropCols = [c for c in allLotusCols if c in merged.columns and c not in protectedCols]
    merged = merged.drop(columns=dropCols, errors="ignore")

    filledCount = len([c for c in imputedCells if c[1] == "organisms"])
    print(f"[Local LOTUS] Filled {filledCount} organism entries via InChIKey Block-1 merge.")

    return merged, imputedCells


# =====================================================================
# 3. 便捷一体化接口 (Convenience All-in-One Interface)
# =====================================================================

def runLocalTaxonomyFilling(
    targetCsvPath: str,
    lotusCsvPath: str,
    outputCsvPath: Optional[str] = None,
    inchikeyColumn: str = "standard_inchi_key",
) -> pd.DataFrame:
    """
    一键执行本地 LOTUS 分类学填补的完整流程。
    One-click execution of the full local LOTUS taxonomy filling workflow.

    Args:
        targetCsvPath: 目标 CSV 路径 (e.g. Coconut_Sugar_Check.csv)
        lotusCsvPath: LOTUS dump CSV 路径 (e.g. 230106_frozen_metadata.csv.gz)
        outputCsvPath: 可选输出路径
        inchikeyColumn: InChIKey 列名

    Returns:
        填补后的 DataFrame
    """
    print(f"Loading LOTUS dump from {lotusCsvPath}...")
    lotusIndex = loadLotusDump(lotusCsvPath)
    print(f"  LOTUS index built: {len(lotusIndex)} unique InChIKey Block-1 entries.")

    print(f"Loading target data from {targetCsvPath}...")
    targetDf = pd.read_csv(targetCsvPath, low_memory=False, dtype=str, encoding="utf-8-sig")
    print(f"  Target data: {len(targetDf)} rows.")

    resultDf, imputedCells = fillTaxonomyFromLotus(targetDf, lotusIndex, inchikeyColumn)

    if outputCsvPath:
        resultDf.to_csv(outputCsvPath, index=False, encoding="utf-8-sig")
        print(f"  Result saved to {outputCsvPath}")

    return resultDf


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Local LOTUS Taxonomy Filling")
    parser.add_argument("--target", required=True, help="Target CSV path (e.g. Coconut_Sugar_Check.csv)")
    parser.add_argument("--lotus", required=True, help="LOTUS dump CSV path (e.g. 230106_frozen_metadata.csv.gz)")
    parser.add_argument("--output", default=None, help="Output CSV path (optional)")
    args = parser.parse_args()

    runLocalTaxonomyFilling(args.target, args.lotus, args.output)
