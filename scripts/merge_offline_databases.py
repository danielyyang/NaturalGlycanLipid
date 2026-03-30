"""
离线数据库推荐与通用合并工具
Offline Database Recommender & Universal Merge Utility

为 GlycoNP 项目推荐 5 个可离线下载的开源数据库, 并提供通用合并函数。
Recommends 5 downloadable open-source databases and provides a universal merge function.

推荐数据库清单 / Recommended Databases:
  1. LOTUS  — 天然产物-物种映射 (NP-Organism mapping), InChIKey merge
  2. COCONUT — 天然产物综合库 (Comprehensive NP DB), InChIKey merge
  3. NPAtlas — 微生物源天然产物 (Microbial NP), InChIKey merge
  4. GlyConnect — 糖链结构数据库 (Glycan structure DB), name/composition merge
  5. GBIF Backbone Taxonomy — 物种分类学层级 (Taxonomy hierarchy), organism name merge

数据下载指南 / Download Instructions:
  详见本文件底部的 DOWNLOAD_GUIDE 字典。

使用方法 / Usage:
  python scripts/merge_offline_databases.py --input reports/GlycoNP_ChEMBL_Enriched.csv \
    --lotus data/external/lotus_frozen.csv.gz \
    --npatlas data/external/npatlas_download.tsv

  或单独调用 mergeLocalTsv() 函数:
  from merge_offline_databases import mergeLocalTsv
  mergedDf = mergeLocalTsv(mainDf, "path/to.tsv", mergeKey="standard_inchi_key")
"""
import argparse
import os
import sys
import time
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


# =====================================================================
# 离线数据库下载指南 (Download Guide)
# =====================================================================

DOWNLOAD_GUIDE: Dict[str, dict] = {
    "LOTUS": {
        "description_zh": "目前最大的开放天然产物-物种映射数据库 (约 80 万条)",
        "description_en": "Largest open NP-organism mapping DB (~800k entries)",
        "url": "https://lotus.naturalproducts.net/download",
        "file": "lotus_frozen.csv.gz (~300MB)",
        "merge_key": "standard_inchi_key (InChIKey 列名: inchikey 或 structure_inchikey)",
        "fills_columns": ["organisms (物种)", "organism_taxonomy (分类学)", "LOTUS_class"],
        "download_steps": [
            "1. 访问 https://lotus.naturalproducts.net/download",
            "2. 下载 'Full frozen LOTUS dataset' (CSV 格式, ~300MB)",
            "3. 放入 data/external/ 目录",
        ],
    },
    "COCONUT": {
        "description_zh": "最大的开放天然产物综合数据库 (约 40 万化合物)",
        "description_en": "Largest open NP comprehensive DB (~400k compounds)",
        "url": "https://coconut.naturalproducts.net/download",
        "file": "COCONUT_DB.csv 或 COCONUT_DB.sdf",
        "merge_key": "standard_inchi_key",
        "fills_columns": ["NP metadata", "chemical class", "superclass"],
        "download_steps": [
            "1. 访问 https://coconut.naturalproducts.net/download",
            "2. 选择 'Complete COCONUT' → 下载 CSV 格式",
            "3. 放入 data/external/ 目录",
        ],
    },
    "NPAtlas": {
        "description_zh": "微生物源天然产物专库 (细菌+真菌, ~33k 化合物)",
        "description_en": "Microbial NP DB (bacteria + fungi, ~33k compounds)",
        "url": "https://www.npatlas.org/download",
        "file": "npatlas_download.tsv (~30MB)",
        "merge_key": "standard_inchi_key (InChIKey 列名: compound_inchikey)",
        "fills_columns": ["source_organism (微生物种属)", "genus", "origin_type"],
        "download_steps": [
            "1. 访问 https://www.npatlas.org/download",
            "2. 下载 'Full Database' TSV 文件",
            "3. 放入 data/external/ 目录",
        ],
    },
    "GlyConnect": {
        "description_zh": "糖链结构-蛋白关联数据库 (Swiss 维护)",
        "description_en": "Glycan structure-protein association DB (maintained by SIB/ExPASy)",
        "url": "https://glyconnect.expasy.org/downloads",
        "file": "glyconnect_compositions.csv",
        "merge_key": "glycan_composition 或 sugar_sequence (按组成匹配)",
        "fills_columns": ["glycan_function", "protein_association"],
        "download_steps": [
            "1. 访问 https://glyconnect.expasy.org/downloads",
            "2. 下载 'Compositions' 数据集",
            "3. 放入 data/external/ 目录",
        ],
    },
    "GBIF_Backbone": {
        "description_zh": "全球生物多样性分类学骨架 (界/门/纲/目/科/属/种)",
        "description_en": "Global Biodiversity Taxonomy Backbone (Kingdom→Species hierarchy)",
        "url": "https://hosted-datasets.gbif.org/datasets/backbone/current/",
        "file": "backbone/Taxon.tsv (~1.5GB, gzipped ~200MB)",
        "merge_key": "organism name (模糊匹配 canonical_name 列)",
        "fills_columns": ["kingdom", "phylum", "class", "order", "family", "genus"],
        "download_steps": [
            "1. 访问 https://hosted-datasets.gbif.org/datasets/backbone/current/",
            "2. 下载 'simple.txt.gz' 或完整的 backbone.zip",
            "3. 解压后取 Taxon.tsv, 放入 data/external/ 目录",
        ],
    },
}


# =====================================================================
# 通用合并函数 (Universal Merge Function)
# =====================================================================

def mergeLocalTsv(
    mainDf: pd.DataFrame,
    downloadedTsvPath: str,
    mergeKey: str = "standard_inchi_key",
    externalKeyCol: Optional[str] = None,
    keepCols: Optional[List[str]] = None,
    prefix: str = "",
    how: str = "left",
    chunkSize: int = 500_000,
) -> Tuple[pd.DataFrame, int]:
    """
    通用离线数据合并函数 — 将下载好的 TSV/CSV 文件与主数据集合并。
    Universal offline data merge function — merge downloaded TSV/CSV with main dataset.

    设计意图: 不依赖任何外部 API, 纯本地文件操作。
    支持大文件分块读取, 避免内存溢出。
    Design: zero network dependency, pure local file I/O.
    Supports chunked reading for large files to prevent OOM.

    Args:
        mainDf: 主数据集 (GlycoNP Pipeline output)
        downloadedTsvPath: 下载的 TSV/CSV 文件路径
        mergeKey: 主数据集中的合并键列名
        externalKeyCol: 外部文件中的合并键列名 (为 None 时与 mergeKey 相同)
        keepCols: 只保留外部文件中这些列 (为 None 时保留全部)
        prefix: 合并后列名前缀 (避免列名冲突)
        how: merge 方式 ("left" 保留全部主数据, "inner" 只保留匹配行)
        chunkSize: 大文件分块读取大小

    Returns:
        (mergedDf, hitCount): 合并后的 DataFrame 和命中行数
    """
    if not os.path.exists(downloadedTsvPath):
        print(f"  [SKIP] File not found: {downloadedTsvPath}")
        return mainDf, 0

    extKeyCol = externalKeyCol or mergeKey
    sep = "\t" if downloadedTsvPath.endswith(".tsv") else ","

    print(f"  Loading: {downloadedTsvPath}")
    t0 = time.time()

    # 探测文件大小, 决定是否分块 (Detect file size for chunking decision)
    fileSize = os.path.getsize(downloadedTsvPath)
    fileSizeMb = fileSize / (1024 * 1024)
    print(f"    File size: {fileSizeMb:.1f} MB")

    # 读取外部数据 (Load external data)
    try:
        if fileSizeMb > 500:
            # 大文件: 分块读取, 仅保留匹配键 (Large: chunked, keep only matches)
            mainKeys = set(mainDf[mergeKey].dropna().astype(str))
            chunks = []
            for chunk in pd.read_csv(
                downloadedTsvPath, sep=sep, dtype=str, low_memory=False,
                chunksize=chunkSize,
            ):
                if extKeyCol in chunk.columns:
                    filtered = chunk[chunk[extKeyCol].isin(mainKeys)]
                    if not filtered.empty:
                        chunks.append(filtered)
            extDf = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
        else:
            # 小文件: 直接读取 (Small: read all at once)
            extDf = pd.read_csv(
                downloadedTsvPath, sep=sep, dtype=str, low_memory=False,
            )
    except Exception as e:
        print(f"    [ERROR] Failed to read: {e}")
        return mainDf, 0

    if extDf.empty or extKeyCol not in extDf.columns:
        print(f"    [WARN] No data or missing key column '{extKeyCol}'")
        return mainDf, 0

    print(f"    External rows: {len(extDf):,} | Key column: '{extKeyCol}'")

    # 筛选列 (Filter columns)
    if keepCols:
        keepCols = [c for c in keepCols if c in extDf.columns]
        extDf = extDf[[extKeyCol] + keepCols]

    # 列名冲突处理 (Handle column name conflicts)
    if extKeyCol != mergeKey:
        extDf = extDf.rename(columns={extKeyCol: mergeKey})

    # 去重 (Deduplicate external data on merge key)
    extDf = extDf.drop_duplicates(subset=[mergeKey])

    # 添加前缀 (Add prefix to avoid column name conflicts)
    if prefix:
        renameCols = {c: f"{prefix}_{c}" for c in extDf.columns if c != mergeKey}
        extDf = extDf.rename(columns=renameCols)

    # 合并 (Merge)
    originalLen = len(mainDf)
    mergedDf = mainDf.merge(extDf, on=mergeKey, how=how)

    # 计算命中数 (Count hits)
    newCols = [c for c in mergedDf.columns if c not in mainDf.columns]
    hitCount = 0
    if newCols:
        hitCount = mergedDf[newCols[0]].notna().sum()

    elapsed = time.time() - t0
    print(f"    Hits: {hitCount:,} / {originalLen:,} ({hitCount/originalLen*100:.1f}%)")
    print(f"    New columns: {newCols[:5]}")
    print(f"    Time: {elapsed:.1f}s")

    return mergedDf, hitCount


# =====================================================================
# 批量合并 (Batch Merge)
# =====================================================================

def batchMerge(
    mainDf: pd.DataFrame,
    externalDir: str,
) -> pd.DataFrame:
    """
    按优先级依次合并所有可用的离线数据库。
    Merge all available offline databases in priority order.
    """
    print(f"\n  Checking for offline databases in: {externalDir}")

    # LOTUS — 最优先, 补全 organisms
    # LOTUS — highest priority, fills organisms
    lotusPath = os.path.join(externalDir, "lotus_frozen.csv.gz")
    if not os.path.exists(lotusPath):
        lotusPath = os.path.join(externalDir, "lotus_frozen.csv")
    if os.path.exists(lotusPath):
        mainDf, hits = mergeLocalTsv(
            mainDf, lotusPath, mergeKey="standard_inchi_key",
            externalKeyCol="structure_inchikey",
            keepCols=["organism_name", "organism_taxonomy_01kingdom",
                       "organism_taxonomy_02phylum",
                       "organism_taxonomy_06family"],
            prefix="LOTUS",
        )
        print(f"    [LOTUS] {hits:,} organisms enriched")
    else:
        print(f"    [LOTUS] Not found: {lotusPath}")

    # NPAtlas — 微生物数据
    # NPAtlas — microbial NP data
    npatlasPath = os.path.join(externalDir, "npatlas_download.tsv")
    if os.path.exists(npatlasPath):
        mainDf, hits = mergeLocalTsv(
            mainDf, npatlasPath, mergeKey="standard_inchi_key",
            externalKeyCol="compound_inchikey",
            keepCols=["genus", "species", "origin_type",
                       "compound_names"],
            prefix="NPAtlas",
        )
        print(f"    [NPAtlas] {hits:,} microbial NPs enriched")
    else:
        print(f"    [NPAtlas] Not found: {npatlasPath}")

    return mainDf


# =====================================================================
# 下载指南打印 (Download Guide Printer)
# =====================================================================

def printDownloadGuide():
    """打印完整的离线数据库下载指南。 / Print complete download guide."""
    print("\n" + "=" * 70)
    print("  离线数据库下载指南 / Offline Database Download Guide")
    print("=" * 70)

    for name, info in DOWNLOAD_GUIDE.items():
        print(f"\n  [{name}]")
        print(f"    {info['description_zh']}")
        print(f"    {info['description_en']}")
        print(f"    URL: {info['url']}")
        print(f"    File: {info['file']}")
        print(f"    Merge Key: {info['merge_key']}")
        print(f"    Fills: {', '.join(info['fills_columns'])}")
        print(f"    Steps:")
        for step in info['download_steps']:
            print(f"      {step}")

    print(f"\n  放置位置: data/external/ 目录")
    print(f"  All files should be placed in: data/external/")
    print("=" * 70)


# =====================================================================
# Main
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Offline Database Merge Utility")
    parser.add_argument("--input", type=str, default=None)
    parser.add_argument("--output", type=str, default=None)
    parser.add_argument("--guide", action="store_true",
                        help="Print download guide and exit")
    parser.add_argument("--lotus", type=str, default=None,
                        help="Path to LOTUS CSV")
    parser.add_argument("--npatlas", type=str, default=None,
                        help="Path to NPAtlas TSV")
    args = parser.parse_args()

    if args.guide:
        printDownloadGuide()
        return

    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    reportDir = os.path.join(baseDir, "reports")
    externalDir = os.path.join(baseDir, "data", "external")

    inputPath = args.input or os.path.join(
        reportDir, "GlycoNP_ChEMBL_Enriched.csv")
    outputPath = args.output or os.path.join(
        reportDir, "GlycoNP_Fully_Enriched.csv")

    print("=" * 70)
    print("  GlycoNP Offline Database Merge Utility")
    print("  糖缀合物离线数据库合并工具")
    print("=" * 70)

    if not os.path.exists(inputPath):
        print(f"  [ERROR] Input not found: {inputPath}")
        print(f"  Run mine_chembl_sqlite.py first.")
        return

    t0 = time.time()
    df = pd.read_csv(inputPath, low_memory=False, dtype=str,
                     encoding="utf-8-sig")
    print(f"  Loaded: {len(df):,} rows")

    # 尝试批量合并 (Attempt batch merge)
    os.makedirs(externalDir, exist_ok=True)
    df = batchMerge(df, externalDir)

    # 如果指定了单独路径 (If specific paths specified)
    if args.lotus:
        df, hits = mergeLocalTsv(
            df, args.lotus, mergeKey="standard_inchi_key",
            externalKeyCol="structure_inchikey",
            prefix="LOTUS",
        )
    if args.npatlas:
        df, hits = mergeLocalTsv(
            df, args.npatlas, mergeKey="standard_inchi_key",
            externalKeyCol="compound_inchikey",
            prefix="NPAtlas",
        )

    # 保存 (Save)
    df.to_csv(outputPath, index=False, encoding="utf-8-sig")
    elapsed = time.time() - t0
    print(f"\n  Output: {outputPath}")
    print(f"  Total time: {elapsed:.0f}s")

    printDownloadGuide()


if __name__ == "__main__":
    main()
