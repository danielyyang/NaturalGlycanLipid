"""
Mock LOTUS 数据生成器 — 用真实 InChIKey 构造模拟分类学数据
Mock LOTUS Data Generator — Build simulated taxonomy from real InChIKeys

功能 (Features):
1. 读取 Coconut_Sugar_Check.csv 中真实的 InChIKey
2. 随机抽取 1,000 个，提取 Block-1 (前 14 位)
3. 为每个分配虚构但格式正确的 LOTUS 分类学字段
4. 部分分子分配 2~3 个物种（测试 | 拼接聚合逻辑）
5. 保存到 data/mock_lotus_metadata.csv
6. 自动调用 runLocalTaxonomyFilling 进行匹配测试
"""
import os
import sys
import time
import random
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# =====================================================================
# [TEST DATA ONLY] — 虚构分类学词库 (Synthetic taxonomy vocabulary)
# =====================================================================

MOCK_DOMAINS = ["Eukaryota", "Bacteria", "Archaea"]

MOCK_KINGDOMS = ["Plantae", "Fungi", "Animalia", "Bacteria", "Chromista"]

# 科名-属名-种名的配套组合 (Family-Genus-Species triplets)
MOCK_TAXONOMY_POOL = [
    ("Araliaceae",       "Panax",          "Panax ginseng"),
    ("Araliaceae",       "Panax",          "Panax notoginseng"),
    ("Araliaceae",       "Panax",          "Panax quinquefolius"),
    ("Fabaceae",         "Glycyrrhiza",    "Glycyrrhiza glabra"),
    ("Fabaceae",         "Astragalus",     "Astragalus membranaceus"),
    ("Fabaceae",         "Trifolium",      "Trifolium pratense"),
    ("Rutaceae",         "Citrus",         "Citrus sinensis"),
    ("Rutaceae",         "Citrus",         "Citrus limon"),
    ("Poaceae",          "Oryza",          "Oryza sativa"),
    ("Poaceae",          "Zea",            "Zea mays"),
    ("Brassicaceae",     "Arabidopsis",    "Arabidopsis thaliana"),
    ("Brassicaceae",     "Brassica",       "Brassica napus"),
    ("Rosaceae",         "Rosa",           "Rosa chinensis"),
    ("Rosaceae",         "Malus",          "Malus domestica"),
    ("Asteraceae",       "Artemisia",      "Artemisia annua"),
    ("Asteraceae",       "Helianthus",     "Helianthus annuus"),
    ("Lamiaceae",        "Salvia",         "Salvia miltiorrhiza"),
    ("Lamiaceae",        "Mentha",         "Mentha piperita"),
    ("Solanaceae",       "Solanum",        "Solanum tuberosum"),
    ("Solanaceae",       "Capsicum",       "Capsicum annuum"),
    ("Ginkgoaceae",      "Ginkgo",         "Ginkgo biloba"),
    ("Apiaceae",         "Angelica",       "Angelica sinensis"),
    ("Berberidaceae",    "Berberis",       "Berberis vulgaris"),
    ("Zingiberaceae",    "Curcuma",        "Curcuma longa"),
    ("Theaceae",         "Camellia",       "Camellia sinensis"),
    ("Streptomycetaceae","Streptomyces",   "Streptomyces griseus"),
    ("Pseudomonadaceae", "Pseudomonas",    "Pseudomonas aeruginosa"),
    ("Saccharomycetaceae","Saccharomyces", "Saccharomyces cerevisiae"),
    ("Aspergillaceae",   "Aspergillus",    "Aspergillus niger"),
    ("Hominidae",        "Homo",           "Homo sapiens"),
]


def generateMockLotus(coconutCsvPath: str, outputPath: str, sampleSize: int = 1000) -> str:
    """
    从真实 COCONUT 数据中抽取 InChIKey 并生成 Mock LOTUS CSV。

    Args:
        coconutCsvPath: Coconut_Sugar_Check.csv 路径
        outputPath: 输出文件路径
        sampleSize: 抽样数量

    Returns:
        输出文件路径
    """
    print(f"Reading real InChIKeys from {coconutCsvPath}...")
    df = pd.read_csv(coconutCsvPath, usecols=["standard_inchi_key"], dtype=str, low_memory=False)
    df = df.dropna(subset=["standard_inchi_key"])

    # 去重并取有效的 InChIKey (27 字符格式)
    # Deduplicate and filter valid InChIKeys
    validKeys = df["standard_inchi_key"].str.strip()
    validKeys = validKeys[validKeys.str.len() == 27].unique()
    print(f"  Found {len(validKeys)} unique valid InChIKeys.")

    # 随机抽样 (Random sample)
    actualSize = min(sampleSize, len(validKeys))
    sampledKeys = random.sample(list(validKeys), actualSize)
    print(f"  Sampled {actualSize} InChIKeys for mock LOTUS generation.")

    # 生成 Mock 记录 (Generate mock records)
    # 约 30% 的分子会有 2~3 条记录（多物种），模拟 LOTUS 中一个化合物来源于多个物种的情况
    # ~30% of molecules will have 2-3 records (multi-species), simulating LOTUS multi-source
    random.seed(42)
    records = []

    for inchikey in sampledKeys:
        # 决定该分子有几个物种来源 (Decide how many species for this molecule)
        numSpecies = random.choices([1, 2, 3], weights=[0.7, 0.2, 0.1], k=1)[0]

        # 随机选择 taxonomy triplets
        chosenTaxa = random.sample(MOCK_TAXONOMY_POOL, min(numSpecies, len(MOCK_TAXONOMY_POOL)))

        for family, genus, species in chosenTaxa:
            domain = random.choice(MOCK_DOMAINS)
            # 根据 family 粗略推断 kingdom
            # Rough kingdom inference from family
            if family in ("Streptomycetaceae", "Pseudomonadaceae"):
                kingdom = "Bacteria"
                domain = "Bacteria"
            elif family in ("Saccharomycetaceae", "Aspergillaceae"):
                kingdom = "Fungi"
            elif family == "Hominidae":
                kingdom = "Animalia"
            else:
                kingdom = "Plantae"
                domain = "Eukaryota"

            records.append({
                "structure_inchikey": inchikey,
                "organism_name": species,
                "organism_taxonomy_01domain": domain,
                "organism_taxonomy_02kingdom": kingdom,
                "organism_taxonomy_05family": family,
                "organism_taxonomy_06genus": genus,
                "organism_taxonomy_08species": species,
            })

    mockDf = pd.DataFrame(records)
    os.makedirs(os.path.dirname(outputPath), exist_ok=True)
    mockDf.to_csv(outputPath, index=False, encoding="utf-8-sig")
    print(f"  Mock LOTUS CSV saved: {outputPath} ({len(mockDf)} records for {actualSize} molecules)")
    return outputPath


def main():
    baseDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    coconutCsv = os.path.join(baseDir, "reports", "Coconut_Sugar_Check.csv")
    mockOutput = os.path.join(baseDir, "data", "mock_lotus_metadata.csv")
    resultOutput = os.path.join(baseDir, "reports", "Coconut_Sugar_MockTax.csv")

    if not os.path.exists(coconutCsv):
        print(f"ERROR: {coconutCsv} not found!")
        return

    # Step 1: 生成 Mock LOTUS CSV
    print("=" * 70)
    print("[Step 1] Generating Mock LOTUS CSV...")
    print("=" * 70)
    generateMockLotus(coconutCsv, mockOutput, sampleSize=1000)

    # Step 2: 使用 taxonomy_lotus_matcher.py 进行匹配测试
    print()
    print("=" * 70)
    print("[Step 2] Running Local Taxonomy Filling with Mock Data...")
    print("=" * 70)

    from lib.taxonomy_lotus_matcher import loadLotusDump, fillTaxonomyFromLotus

    startTime = time.time()

    # 加载 mock LOTUS 索引 (Load mock LOTUS index)
    lotusIndex = loadLotusDump(mockOutput)
    loadTime = time.time() - startTime
    print(f"  LOTUS index load time: {loadTime:.3f}s")
    print(f"  LOTUS index entries: {len(lotusIndex)}")

    # 加载目标数据 (Load target data)
    targetDf = pd.read_csv(coconutCsv, low_memory=False, dtype=str, encoding="utf-8-sig")

    # 记录匹配前的空缺数 (Count missing before)
    beforeMissing = (
        targetDf["organisms"].isna()
        | targetDf["organisms"].astype(str).str.strip().isin(["", "nan", "NULL", "Not Result"])
    ).sum()

    mergeStart = time.time()
    resultDf, imputedCells = fillTaxonomyFromLotus(targetDf, lotusIndex)
    mergeTime = time.time() - mergeStart

    # 记录匹配后的结果 (Count results after)
    afterMissing = (
        resultDf["organisms"].isna()
        | resultDf["organisms"].astype(str).str.strip().isin(["", "nan", "NULL", "Not Result"])
    ).sum()

    orgImputedCount = len([c for c in imputedCells if c[1] == "organisms"])
    famImputedCount = len([c for c in imputedCells if c[1] == "Family"])

    # 保存结果 (Save results)
    resultDf.to_csv(resultOutput, index=False, encoding="utf-8-sig")

    # 输出统计 (Print statistics)
    totalTime = time.time() - startTime
    print()
    print("=" * 70)
    print("[Results] Mock LOTUS Taxonomy Filling Statistics")
    print("=" * 70)
    print(f"  Total rows:              {len(resultDf)}")
    print(f"  Missing before fill:     {beforeMissing}")
    print(f"  Organisms filled:        {orgImputedCount}")
    print(f"  Families filled:         {famImputedCount}")
    print(f"  Missing after fill:      {afterMissing}")
    print(f"  Fill rate (organism):    {orgImputedCount / max(beforeMissing, 1) * 100:.1f}%")
    print(f"  ---")
    print(f"  LOTUS index load time:   {loadTime:.3f}s")
    print(f"  Merge time (94K rows):   {mergeTime:.3f}s")
    print(f"  Total time:              {totalTime:.3f}s")
    print(f"  Result saved to:         {resultOutput}")
    print()

    # 展示几条填补示例 (Show a few fill examples)
    filledRows = resultDf.loc[
        resultDf.index.isin([c[0] for c in imputedCells if c[1] == "organisms"])
    ]
    if not filledRows.empty:
        print("[Sample Filled Rows]")
        sample = filledRows.head(5)[["standard_inchi_key", "organisms", "Family"]].copy()
        sample["inchikey_block1"] = sample["standard_inchi_key"].str[:14]
        print(sample.to_string(index=False))


if __name__ == "__main__":
    main()
