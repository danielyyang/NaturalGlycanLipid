"""
PubChem InChIKey 批量验证与数据补全 (PubChem InChIKey Batch Enrichment)
=====================================================================
用 InChIKey 前 14 位 (molecular layer) 查询 PubChem PUG REST:
  1. CID + IsomericSMILES + IUPACName
  2. 物种 (Taxonomy) 信息 — 如果 PubChem 有 BioAssay 关联
  3. 同义名 (Synonyms) — 用于数据补全

用法 (Usage):
    python scripts/pubchem_enrichment.py [--limit 500] [--delay 0.3]

[TEST DATA ONLY]
"""
import argparse
import json
import os
import sys
import time
from typing import Optional, Dict, Tuple

import pandas as pd
import requests
from tqdm import tqdm

sys.path.insert(0, r"d:\Glycan_Database")

BASE_DIR = r"d:\Glycan_Database"
REPORT_DIR = os.path.join(BASE_DIR, "reports")
INPUT_CSV = os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched_v12.csv")
OUTPUT_CSV = os.path.join(REPORT_DIR, "GlycoNP_Deep_Enriched_v12.csv")

# PubChem PUG REST endpoints
PUG_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def queryPubchemByInchikey(
    inchikey: str,
    delay: float = 0.3,
) -> Optional[Dict]:
    """用 InChIKey 查 PubChem 获取基础记录 (两步法: InChIKey→CID→属性).
    Query PubChem by InChIKey using 2-step approach: InChIKey→CID→Properties.

    两步法原因: PUG REST /property/ endpoint 对 InChIKey 直查支持有限,
    部分属性 (如 IUPACName) 会返回 400。先拿 CID 再查属性更稳定。
    Reason: PUG REST /property/ endpoint has limited support for InChIKey
    direct queries. Getting CID first then querying by CID is more reliable.

    Args:
        inchikey: 完整 InChIKey (27 chars) 或前 14 位 (connectivity layer)
        delay: 请求间隔秒数 (避免 429 限流)

    Returns:
        字典: {"cid": int, "iupac_name": str, "isomeric_smiles": str,
               "synonyms": list, "molecular_formula": str}
        或 None (查询失败)
    """
    if not inchikey or str(inchikey) in ("nan", "None", ""):
        return None

    # 清洗 InChIKey (clean up)
    inchikey = str(inchikey).strip()
    if len(inchikey) < 14:
        return None

    try:
        # Step 1: InChIKey → CID
        cidUrl = f"{PUG_BASE}/compound/inchikey/{inchikey}/cids/JSON"
        time.sleep(delay)
        cidResp = requests.get(cidUrl, timeout=15)

        if cidResp.status_code == 404:
            return None
        if cidResp.status_code == 429:
            time.sleep(5)
            cidResp = requests.get(cidUrl, timeout=15)
        if cidResp.status_code != 200:
            return None

        cidData = cidResp.json()
        cidList = cidData.get("IdentifierList", {}).get("CID", [])
        if not cidList:
            return None
        cid = cidList[0]

        # Step 2: CID → 属性 (Properties)
        propUrl = f"{PUG_BASE}/compound/cid/{cid}/property/IsomericSMILES,MolecularFormula,IUPACName/JSON"
        time.sleep(delay)
        propResp = requests.get(propUrl, timeout=15)

        result = {
            "cid": cid,
            "iupac_name": "",
            "isomeric_smiles": "",
            "molecular_formula": "",
            "synonyms": [],
        }

        if propResp.status_code == 200:
            propData = propResp.json()
            props = propData.get("PropertyTable", {}).get("Properties", [])
            if props:
                prop = props[0]
                result["iupac_name"] = prop.get("IUPACName", "")
                result["isomeric_smiles"] = prop.get("IsomericSMILES", "")
                result["molecular_formula"] = prop.get("MolecularFormula", "")

        # Step 3: 获取同义名 (Get synonyms)
        try:
            time.sleep(delay)
            synUrl = f"{PUG_BASE}/compound/cid/{cid}/synonyms/JSON"
            synResp = requests.get(synUrl, timeout=10)
            if synResp.status_code == 200:
                synData = synResp.json()
                synList = synData.get("InformationList", {}).get("Information", [])
                if synList:
                    result["synonyms"] = synList[0].get("Synonym", [])[:10]
        except Exception:
            pass

        return result

    except Exception:
        return None


def enrichFromPubchem(
    df: pd.DataFrame,
    limit: Optional[int] = None,
    delay: float = 0.3,
) -> pd.DataFrame:
    """用 InChIKey 批量查 PubChem 并补全数据.
    Batch query PubChem by InChIKey and enrich the DataFrame.

    补全策略 (Enrichment strategy):
    1. 如果 iupac_name 为空 → 从 PubChem 补全
    2. 如果 name 为空 → 从 synonyms[0] 补全
    3. 记录 PubChem CID 用于后续引用
    4. 如果 SMILES 不一致 → 标记但不覆盖 (保留我们的)

    Args:
        df: 输入 DataFrame (必须含 InChIKey 列)
        limit: 最多查询行数 (默认全部)
        delay: 请求间隔秒数

    Returns:
        补全后的 DataFrame
    """
    print("\n" + "=" * 70)
    print("  PubChem InChIKey Batch Enrichment (2-Step CID Method)")
    print("=" * 70)

    # 自动检测 InChIKey 列名, 含 2D InChIKey 降级
    # Auto-detect InChIKey column, with 2D InChIKey fallback
    inchikeyCol = None
    for candidate in ["inchikey", "standard_inchi_key", "InChIKey", "inchi_key",
                       "non_isomeric_inchikey_block1"]:
        if candidate in df.columns:
            inchikeyCol = candidate
            break
    if inchikeyCol is None:
        print("  [ERROR] No InChIKey column found!")
        return df
    print(f"  InChIKey column: {inchikeyCol}")

    # 确保新列存在 (Ensure new columns exist)
    for col in ["PubChem_CID", "PubChem_IUPAC", "PubChem_Synonyms", "PubChem_SMILES_Match"]:
        if col not in df.columns:
            df[col] = ""

    # 选取有 InChIKey 的行 (Select rows with InChIKey)
    hasInchi = df[inchikeyCol].notna() & (df[inchikeyCol].astype(str) != "nan") & (df[inchikeyCol].astype(str) != "")
    targetIdx = df.index[hasInchi]
    if limit:
        targetIdx = targetIdx[:limit]

    print(f"  Total InChIKey rows: {hasInchi.sum():,}")
    print(f"  Query target: {len(targetIdx):,}")
    print(f"  Delay: {delay}s per request")
    print(f"  Estimated time: {len(targetIdx) * delay * 2 / 60:.0f} min")

    enriched = 0
    namesFilled = 0
    smilesMatch = 0
    smilesMismatch = 0
    errorCount = 0

    for idx in tqdm(targetIdx, desc="  PubChem", ncols=80):
        inchi = str(df.at[idx, inchikeyCol])
        try:
            result = queryPubchemByInchikey(inchi, delay=delay)
            if result is None:
                continue

            # CID
            if result["cid"]:
                df.at[idx, "PubChem_CID"] = str(result["cid"])

            # IUPAC Name 补全 (IUPAC Name fill)
            if result["iupac_name"]:
                if "iupac_name" not in df.columns:
                    df["iupac_name"] = ""
                currentIupac = str(df.at[idx, "iupac_name"]) if pd.notna(df.at[idx, "iupac_name"]) else ""
                if not currentIupac or currentIupac == "nan":
                    df.at[idx, "iupac_name"] = result["iupac_name"]
                df.at[idx, "PubChem_IUPAC"] = result["iupac_name"]

            # Name 补全 (Name fill from synonyms)
            if result["synonyms"]:
                df.at[idx, "PubChem_Synonyms"] = "; ".join(result["synonyms"][:5])
                currentName = str(df.at[idx, "name"]) if pd.notna(df.at[idx, "name"]) else ""
                if not currentName or currentName == "nan":
                    df.at[idx, "name"] = result["synonyms"][0]
                    namesFilled += 1

            # SMILES 比对 (SMILES comparison)
            pcSmi = result.get("isomeric_smiles", "")
            ourSmi = str(df.at[idx, "canonical_smiles"]) if pd.notna(df.at[idx, "canonical_smiles"]) else ""
            if pcSmi and ourSmi:
                from rdkit import Chem
                pcMol = Chem.MolFromSmiles(pcSmi)
                ourMol = Chem.MolFromSmiles(ourSmi)
                if pcMol and ourMol:
                    pcCanon = Chem.MolToSmiles(pcMol, isomericSmiles=True)
                    ourCanon = Chem.MolToSmiles(ourMol, isomericSmiles=True)
                    if pcCanon == ourCanon:
                        df.at[idx, "PubChem_SMILES_Match"] = "Match"
                        smilesMatch += 1
                    else:
                        df.at[idx, "PubChem_SMILES_Match"] = "Mismatch"
                        smilesMismatch += 1

            enriched += 1

        except Exception:
            errorCount += 1

    print(f"\n  Results:")
    print(f"    Enriched:       {enriched:,}")
    print(f"    Names filled:   {namesFilled:,}")
    print(f"    SMILES match:   {smilesMatch:,}")
    print(f"    SMILES mismatch:{smilesMismatch:,}")
    print(f"    Errors:         {errorCount:,}")

    return df


def main():
    parser = argparse.ArgumentParser(description="PubChem InChIKey Batch Enrichment")
    parser.add_argument("--limit", type=int, default=None, help="Max rows to query")
    parser.add_argument("--delay", type=float, default=0.3, help="Delay between requests (seconds)")
    args = parser.parse_args()

    print(f"Loading: {INPUT_CSV}")
    df = pd.read_csv(INPUT_CSV, low_memory=False)
    print(f"Rows: {len(df):,}")

    df = enrichFromPubchem(df, limit=args.limit, delay=args.delay)

    print(f"\nSaving: {OUTPUT_CSV}")
    df.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")
    print("Done!")


if __name__ == "__main__":
    main()
