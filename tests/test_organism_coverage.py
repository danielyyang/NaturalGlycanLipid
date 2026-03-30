"""
诊断脚本：测试 COCONUT 数据库中 organism 搜索覆盖率
Diagnostic script: test organism search coverage in COCONUT database

采样前 50 条缺失 organism 的记录，依次尝试 4 个搜索层级：
1. COCONUT API (InChIKey)
2. LOTUS/Wikidata SPARQL (InChIKey)
3. PubChem Taxonomy (iupac_name / name)
4. Wikipedia NLP (name)

输出每层命中数和总搜索成功率。
"""
import os
import sys
import time
import pandas as pd

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from lib.taxonomy_online_resolver import (
    get_organism_from_coconut,
    get_organism_from_wikidata,
    get_organism_from_pubchem,
    get_organism_from_wikipedia,
)

def main():
    csvPath = os.path.join(os.path.dirname(__file__), "..", "reports", "Coconut_Sugar_Check.csv")
    if not os.path.exists(csvPath):
        print(f"❌ 文件未找到: {csvPath}")
        return

    df = pd.read_csv(csvPath, low_memory=False, dtype=str, encoding='utf-8-sig')
    print(f"📊 总行数: {len(df)}")

    # 统计已有 organism 的行数 (Count rows with existing organisms)
    hasOrganism = df['organisms'].notna() & (df['organisms'].str.strip() != '') & (df['organisms'].str.strip() != 'nan')
    missingOrganism = ~hasOrganism

    print(f"✅ 已有 organism: {hasOrganism.sum()} ({hasOrganism.sum()/len(df)*100:.1f}%)")
    print(f"❌ 缺失 organism: {missingOrganism.sum()} ({missingOrganism.sum()/len(df)*100:.1f}%)")
    print()

    # 采样缺失行 (Sample missing rows for testing)
    SAMPLE_SIZE = 30
    missingDf = df[missingOrganism].head(SAMPLE_SIZE).copy()
    print(f"🔬 采样 {len(missingDf)} 条缺失记录进行搜索测试...\n")

    # 按层级统计 (Track hits per tier)
    tierHits = {"COCONUT": 0, "LOTUS/Wikidata": 0, "PubChem": 0, "Wikipedia": 0}
    totalFound = 0
    results = []

    for idx, row in missingDf.iterrows():
        name = str(row.get('name', ''))
        iupac = str(row.get('iupac_name', ''))
        inchikey = str(row.get('standard_inchi_key', ''))
        queryName = name if name and name != 'nan' else iupac

        result = {"name": name[:40], "inchikey": inchikey[:20], "tier": "MISS", "organism": "—"}

        # Tier 1: COCONUT
        if inchikey and inchikey != 'nan':
            org = get_organism_from_coconut(inchikey)
            if org:
                tierHits["COCONUT"] += 1
                totalFound += 1
                result["tier"] = "COCONUT"
                result["organism"] = org[:50]
                results.append(result)
                continue

        # Tier 2: LOTUS/Wikidata
        if inchikey and inchikey != 'nan':
            org = get_organism_from_wikidata(inchikey)
            if org:
                tierHits["LOTUS/Wikidata"] += 1
                totalFound += 1
                result["tier"] = "LOTUS"
                result["organism"] = org[:50]
                results.append(result)
                continue

        # Tier 3: PubChem
        if queryName and queryName != 'nan':
            org = get_organism_from_pubchem(queryName)
            if org:
                tierHits["PubChem"] += 1
                totalFound += 1
                result["tier"] = "PubChem"
                result["organism"] = org[:50]
                results.append(result)
                continue

        # Tier 4: Wikipedia
        if queryName and queryName != 'nan':
            org = get_organism_from_wikipedia(queryName)
            if org:
                tierHits["Wikipedia"] += 1
                totalFound += 1
                result["tier"] = "Wikipedia"
                result["organism"] = org[:50]
                results.append(result)
                continue

        results.append(result)

    # 输出结果 (Print results)
    print("=" * 90)
    print(f"{'#':>3} | {'化合物名':40} | {'搜索层':12} | {'Organism'}")
    print("-" * 90)
    for i, r in enumerate(results, 1):
        tierDisplay = r['tier']
        nameDisplay = r['name'] if r['name'] != 'nan' else '(无名称)'
        print(f"{i:3} | {nameDisplay:40} | {tierDisplay:12} | {r['organism']}")

    print("=" * 90)
    print(f"\n📊 搜索结果摘要 (Search Results Summary):")
    print(f"  采样数: {len(missingDf)}")
    print(f"  搜索命中总数: {totalFound} ({totalFound/len(missingDf)*100:.1f}%)")
    print()
    print(f"  按层级命中数 (Hits per Tier):")
    for tier, hits in tierHits.items():
        print(f"    {tier:20}: {hits:3} ({hits/len(missingDf)*100:.1f}%)")
    print()
    print(f"  搜索失败数: {len(missingDf) - totalFound} ({(len(missingDf) - totalFound)/len(missingDf)*100:.1f}%)")


if __name__ == "__main__":
    main()
