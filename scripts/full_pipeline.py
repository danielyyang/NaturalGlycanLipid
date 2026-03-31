"""
GlycoNP Full Pipeline
=====================
Unified Multiprocessing Pipeline for Glycoside Analysis.
Handles:
  1. Glycan/Aglycon SMILES Cleavage
  2. Murcko Scaffold Generation
  3. Sugar Sequence & Topological Feature Extraction
  4. Glycosidic Bond Linkage Details (reducing-end & inter-sugar)
  5. NLP-based Sequence Rescue

[TEST DATA ONLY]
"""
import os
import sys
import time
import json
import re
import argparse
import pandas as pd
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
from collections import Counter
from rdkit import Chem

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from lib.glycan_topology import get_split_smiles, find_mapped_sugar_units
from lib.feature_extractor import getTopologyScaffoldSmiles
from lib.monosaccharide_identifier import generate_refined_sequence

# Import core robust functions directly from pipeline_utils instead of legacy scripts
from lib.pipeline_utils import (
    extractSugarFromName,
    crossValidateRareSugars,
    detectAllGlycosidicBonds,
    inferOrganismType,
    molToBase64Png,
    molToHighlightedBase64Png,
    GENERIC_PAT,
)

def run_molecule(row):
    """
    [EN] Worker function to process a single molecule for multiprocessing Pool.
         Executes the core pipeline constraints: Cleavage -> Topology -> NLP Rescue.
         Avoids global states to ensure thread safety during high-throughput execution.
         
    [CN] 多线程池处理单分子的工作函数。
         严格按照 剥离 -> 拓扑解析 -> NLP保底 的生命周期执行。
         隔离全局变量以确保在大规模高通量并发时的线程安全。
    """
    res = {
        'identifier': row.get('identifier', ''),
        'canonical_smiles': row.get('canonical_smiles', ''),
        'name': row.get('name', ''),
        'synonyms': row.get('synonyms', ''),
        'Consensus_Sugar_Sequence': '',
        'Aglycon_SMILES': '',
        'Glycan_SMILES': '',
        'Glycan_Modifications': '',
        'Murcko_Scaffold': '',
        'Total_Sugar_Count': 0,
        'Max_Chain_Length': 0,
        'Sub_Chain_Lengths': '',
        'Glycan-Aglycone_Bond_Detail': '[]',
        'Aglycone_Linkage_Type': '',
        'Organism_Type': '',
    }

    try:
        smi = str(row.get('canonical_smiles', ''))
        if not smi or smi in ("nan", "None", ""):
            return res
            
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            return res
            
        # 1. Glycosidic Cleavage (Glycan & Aglycone)
        aglycon_smi, glycan_smi = get_split_smiles(mol)
        
        # Pure sugar retention logic from V12/pipeline_utils
        if glycan_smi and glycan_smi == smi:
            from lib.pipeline_utils import _checkPureSugarMolecule
            if _checkPureSugarMolecule(mol):
                aglycon_smi = ""
            else:
                glycan_smi = ""
                aglycon_smi = ""
                
        res['Glycan_SMILES'] = glycan_smi
        res['Aglycon_SMILES'] = aglycon_smi
        
        # 2. Murcko Scaffold Generation
        if aglycon_smi and len(aglycon_smi) > 4:
            try:
                am = Chem.MolFromSmiles(aglycon_smi)
                if am:
                    res['Murcko_Scaffold'] = getTopologyScaffoldSmiles(aglycon_smi)
            except Exception:
                pass
                
        # 3. Sugar Sequence & Topology Parsing
        units = find_mapped_sugar_units(mol)
        res['Total_Sugar_Count'] = len(units)
        
        raw_seq, mods_str = generate_refined_sequence(mol)
        res['Glycan_Modifications'] = mods_str
        
        # Calculate sub lengths
        chains = [c.strip() for c in raw_seq.split(";")]
        sub_lens = []
        for c in chains:
            if not c or c in ("Non_Cyclic_Invalid", "Invalid", "Error"):
                continue
            # basic heuristic to count monomers in chain string
            sugars = re.findall(r"(?<![A-Za-z])([DL]-(?:d|me|6d)?[A-Z][a-z]+[A-Za-z]*|Hex|Pen|dHex|HexA|Oct|Hept|Non)(?![a-z])", re.sub(r"\[.*?\]-?", "", c))
            sub_lens.append(len(sugars))
            
        res['Sub_Chain_Lengths'] = str(sub_lens)
        res['Max_Chain_Length'] = max(sub_lens) if sub_lens else 0
        
        # 4. NLP Rescue
        name_extracted_sugar = extractSugarFromName(
            row.get('name'), 
            row.get('iupac_name'), 
            row.get('synonyms'), 
            row.get('PMC_Abstract')
        )
        final_seq = crossValidateRareSugars(raw_seq, name_extracted_sugar)
        
        # NLP Backfill for generic tags
        if name_extracted_sugar and bool(GENERIC_PAT.search(final_seq)):
            final_seq = re.sub(r'\b(Hex|Pen|dHex)\b', f"{name_extracted_sugar}(Rescue)", final_seq)
            
        # Anti-regression logic: 
        # [EN] Protect previous non-generic specific identifications from degrading into Hex.
        # [CN] 防回退逻辑：如果新的 2D 重算结果由于手性缺失降级为泛指糖 (Hex)，
        #      但原有管线数据保留了具体糖名（如 D-Glc），则强制沿用有效的高质量就数据。
        oldSeq = str(row.get('Sugar_Sequence', ''))
        if oldSeq and oldSeq not in ("nan", "None", ""):
            newHasGeneric = bool(GENERIC_PAT.search(final_seq))
            oldHasGeneric = bool(GENERIC_PAT.search(oldSeq))
            if newHasGeneric and not oldHasGeneric:
                final_seq = oldSeq
                
        res['Consensus_Sugar_Sequence'] = final_seq

        # 5. Glycosidic Bond Analysis
        rootBond, detailsJson = detectAllGlycosidicBonds(mol, units)
        res['Glycan-Aglycone_Bond_Detail'] = detailsJson
        res['Aglycone_Linkage_Type'] = rootBond
        
        # 6. Organism Type Inference
        res['Organism_Type'] = inferOrganismType(
            row.get('organisms'), 
            row.get('LOTUS_family') or row.get('LOTUS_Family'), 
            row.get('LOTUS_kingdom'), 
            row.get('LOTUS_phylum')
        )

        return res
    except Exception as e:
        return res


def print_stats(df):
    """
    [EN] Print sequence distribution summary. Extracts all sugar tags to calculate Generic vs Specific ratios.
    [CN] 序列分布统计台打印。通过正则抽提所有糖标签，计算泛指糖与具体糖的占比以供质量审核。
    """
    print("\n" + "=" * 70)
    print("  Pipeline Output Statistics")
    print("=" * 70)
    PAT = (r'Neu5Ac|Neu5Gc|KDO|'
           r'[DL]-[A-Z][a-z]+[A-Z]?[a-z]*(?:\([^)]*\))?|'
           r'Hex|dHex|Pen|HexA|HexN|HexNAc|Non|Oct|Hept')
    allTokens = []
    if "Consensus_Sugar_Sequence" in df.columns:
        seq_col = "Consensus_Sugar_Sequence"
    elif "Sugar_Sequence" in df.columns:
        seq_col = "Sugar_Sequence"
    else:
        return
        
    for seq in df[seq_col].dropna():
        allTokens.extend(re.findall(PAT, str(seq)))
        
    tc = Counter(allTokens)
    total = sum(tc.values())
    for rank, (sugar, count) in enumerate(tc.most_common(20), 1):
        pct = count / total * 100 if total > 0 else 0
        print(f"  {rank:2d}. {sugar:<35} {count:>8,} ({pct:5.1f}%)")
        
    genericCount = sum(v for k, v in tc.items() 
                       if k in ("Hex","Pen","dHex","HexA","Non","Oct","Hept"))
    print(f"\n  Generic:        {genericCount:>8,}")
    print(f"  Total:          {total:>8,}")


def main():
    parser = argparse.ArgumentParser(description="GlycoNP Unified Full Pipeline")
    parser.add_argument("--input", "-i", type=str, default="reports/GlycoNP_Deep_Enriched_Input.csv", help="Input CSV path")
    parser.add_argument("--recovered", "-r", type=str, default="reports/Recovered_Sugars_Input.csv", help="Recovered sugars CSV path (optional to merge)")
    parser.add_argument("--output", "-o", type=str, default="reports/GlycoNP_Deep_Enriched_Final.csv", help="Output CSV path")
    parser.add_argument("--cores", "-c", type=int, default=max(1, cpu_count() - 2), help="Number of CPU cores to use")
    args = parser.parse_args()
    
    BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    input_csv = os.path.join(BASE_DIR, args.input) if not os.path.isabs(args.input) else args.input
    recovered_csv = os.path.join(BASE_DIR, args.recovered) if not os.path.isabs(args.recovered) else args.recovered
    output_csv = os.path.join(BASE_DIR, args.output) if not os.path.isabs(args.output) else args.output

    print("=" * 70)
    print("  GlycoNP Unified Full Pipeline")
    print("=" * 70)
    
    if not os.path.exists(input_csv):
        print(f"[!] Error: Input file not found: {input_csv}")
        print("    If you are migrating from V12/V13, please rename your CSV file to match,")
        print("    or explicitly run with: python scripts/full_pipeline.py -i <your_input.csv> -o <your_output.csv>")
        sys.exit(1)

    print(f"Loading Base Data: {input_csv}")
    df_base = pd.read_csv(input_csv, low_memory=False)
    
    if os.path.exists(recovered_csv):
        print(f"Loading Recovered Sugars: {recovered_csv}")
        df_rec = pd.read_csv(recovered_csv, low_memory=False)
        print(f"Merging {len(df_rec)} additional rows...")
        df_combined = pd.concat([df_base, df_rec], ignore_index=True).drop_duplicates(subset=['canonical_smiles'])
    else:
        df_combined = df_base
        
    print(f"Total Unique Molecules: {len(df_combined):,}")

    records = df_combined.to_dict('records')
    
    print(f"\nRunning Parallel Execution on {len(records):,} molecules using {args.cores} cores...")
    t0 = time.time()
    
    results = []
    # Multiprocessing pool execution
    with Pool(args.cores) as p:
        for res in tqdm(p.imap(run_molecule, records, chunksize=1000), total=len(records), ncols=80):
            if res.get('canonical_smiles'):
                results.append(res)
                
    elapsed = time.time() - t0
    print(f"\nProcessing completed in {elapsed:.1f}s ({elapsed/60:.2f} min)")
    
    df_results = pd.DataFrame(results)
    
    # Merge pipeline results back to the original database metadata
    # We drop from df_combined any columns that df_results will overwrite (except identifiers)
    columns_to_drop = [c for c in df_results.columns 
                       if c in df_combined.columns and c not in ['identifier', 'canonical_smiles']]
    df_combined_clean = df_combined.drop(columns=columns_to_drop)
    
    # Use left join in case some were skipped instead of inner? Inner is fine since run_molecule always returns res.
    df_final = pd.merge(df_combined_clean, df_results, on=['identifier', 'canonical_smiles'], how='inner')
    
    # Ensure backward compatible column mapping
    if "Consensus_Sugar_Sequence" in df_final.columns:
        df_final["Sugar_Sequence"] = df_final["Consensus_Sugar_Sequence"]
    if "Glycan-Aglycone_Bond_Detail" in df_final.columns:
        df_final["Bond_Detail"] = df_final["Glycan-Aglycone_Bond_Detail"]
        
    # -----------------------------------------------------------------
    # Dataset-Wide Text-Mining NLP Rescue (Rescue Generic Hex)
    # -----------------------------------------------------------------
    try:
        from scripts.rescue_generic_sugars import rescueSequence
        print("\nRunning Dataset-Wide NLP Rescue (Rescue Generic Hex via Text-Mining)...")
        
        rescuedSeqs = []
        rsCounter = {"A": 0, "MISS": 0}
        
        def safe_get(r, col):
            val = r.get(col)
            return val if pd.notna(val) else None
            
        for idx, r in tqdm(df_final.iterrows(), total=len(df_final), desc="Rescuing Generics", ncols=80):
            seq = str(r.get("Sugar_Sequence", ""))
            if GENERIC_PAT.search(seq):
                newS = rescueSequence(
                    seq,
                    safe_get(r, "name"),
                    safe_get(r, "iupac_name"),
                    safe_get(r, "synonyms"),
                    rsCounter
                )
                rescuedSeqs.append(newS)
            else:
                rescuedSeqs.append(seq)
                
        df_final["Sugar_Sequence"] = rescuedSeqs
        df_final["Consensus_Sugar_Sequence"] = rescuedSeqs
        
        print(f"  -> Rescued by Text-Mining:  {rsCounter['A']:,}")
        print(f"  -> Unrescued (MISS):        {rsCounter['MISS']:,}")
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"\n[WARN] Statistical rescue skipped or failed: {e}")

    print_stats(df_final)
    
    # Save the output CSV
    out_dir = os.path.dirname(output_csv)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    df_final.to_csv(output_csv, index=False, encoding='utf-8-sig')
    print(f"\nSaved Unified DB: {output_csv}")
    print("=" * 70)


if __name__ == "__main__":
    main()
