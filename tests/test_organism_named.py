"""
Targeted test: sample NAMED compounds missing organism to assess API hit rate.
"""
import os, sys, pandas as pd
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from lib.taxonomy_online_resolver import (
    get_organism_from_coconut, get_organism_from_wikidata,
    get_organism_from_pubchem, get_organism_from_wikipedia
)

def main():
    csvPath = os.path.join(os.path.dirname(__file__), "..", "reports", "Coconut_Sugar_Check.csv")
    df = pd.read_csv(csvPath, low_memory=False, dtype=str, encoding='utf-8-sig')

    hasOrg = df['organisms'].notna() & (df['organisms'].str.strip() != '') & (df['organisms'].str.strip() != 'nan')
    missing = df[~hasOrg].copy()
    hasName = missing['name'].notna() & (missing['name'].str.strip() != '') & (missing['name'].str.strip() != 'nan')
    namedMissing = missing[hasName]

    print(f"Total rows: {len(df)}")
    print(f"Has organism: {hasOrg.sum()} ({hasOrg.sum()/len(df)*100:.1f}%)")
    print(f"Missing organism: {len(missing)}")
    print(f"Missing but with name: {len(namedMissing)}")
    print(f"Missing without any name: {len(missing) - len(namedMissing)}")
    print()

    SAMPLE_SIZE = 20
    sample = namedMissing.head(SAMPLE_SIZE)
    tierHits = {"COCONUT": 0, "LOTUS": 0, "PubChem": 0, "Wikipedia": 0}
    total = 0

    for i, (idx, row) in enumerate(sample.iterrows(), 1):
        name = str(row['name'])
        iupac = str(row.get('iupac_name', ''))
        ik = str(row.get('standard_inchi_key', ''))
        qn = name if name != 'nan' else iupac
        foundOrg = None
        tier = 'MISS'

        if ik and ik != 'nan':
            foundOrg = get_organism_from_coconut(ik)
            if foundOrg:
                tier = 'COCONUT'

        if not foundOrg and ik and ik != 'nan':
            foundOrg = get_organism_from_wikidata(ik)
            if foundOrg:
                tier = 'LOTUS'

        if not foundOrg and qn and qn != 'nan':
            foundOrg = get_organism_from_pubchem(qn)
            if foundOrg:
                tier = 'PubChem'

        if not foundOrg and qn and qn != 'nan':
            foundOrg = get_organism_from_wikipedia(qn)
            if foundOrg:
                tier = 'Wikipedia'

        if foundOrg:
            tierHits[tier] += 1
            total += 1

        orgDisplay = str(foundOrg)[:50] if foundOrg else "-"
        nameDisplay = name[:48]
        print(f"{i:2}. [{tier:10}] {nameDisplay:50} => {orgDisplay}")

    print()
    print(f"=== RESULT: {total}/{len(sample)} found ({total/len(sample)*100:.0f}%) ===")
    for t, h in tierHits.items():
        print(f"  {t}: {h}")

if __name__ == "__main__":
    main()
