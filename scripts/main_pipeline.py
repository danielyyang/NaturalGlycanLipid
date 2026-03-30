import pandas as pd
import os
import sys
import argparse
from rdkit import Chem
from tqdm import tqdm

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from lib import glycan_topology
# [REMOVED] Dead import — lib.visualizer no longer exists
# Visualization is now handled by lib.molecular_visualizer
from lib import taxonomy_online_resolver
from lib import monosaccharide_identifier
from rdkit.Chem.Scaffolds import MurckoScaffold
import shutil

def phase1_initialization_and_deduplication(raw_csv_path, sugar_csv_path, check_csv_path):
    """
    Phase 1: 读取全量原始数据，利用 contain_sugar 提取含糖子集，存为 Coconut_Sugar.csv。
    使用 InChIKey 进行精准去重，生成并后续操作于 Coconut_Sugar_Check.csv。
    """
    print("--- [Phase 1] Data Initialization & Deduplication ---")
    
    if not os.path.exists(raw_csv_path):
        raise FileNotFoundError(f"Raw data file not found: {raw_csv_path}")
        
    print(f"Loading raw data from {raw_csv_path}...")
    df_raw = pd.read_csv(raw_csv_path, low_memory=False, dtype=str, encoding='utf-8-sig')
    
    # 粗筛含糖数据
    df_sugar = df_raw[df_raw['contains_sugar'].astype(str).str.lower() == 'true'].copy()
    print(f"Filtered {len(df_sugar)} rows containing sugar. Saving to {sugar_csv_path}...")
    df_sugar.to_csv(sugar_csv_path, index=False, encoding='utf-8-sig')
    
    # 精准去重 (保留第一个)
    is_duplicate = df_sugar.duplicated(subset=['standard_inchi_key'], keep='first')
    df_dedup = df_sugar[~is_duplicate].copy()
    
    print(f"Deduplicated {is_duplicate.sum()} rows. Remaining {len(df_dedup)} unique compounds. Saving to {check_csv_path}...")
    df_dedup.to_csv(check_csv_path, index=False, encoding='utf-8-sig')
    
    return df_dedup

def main():
    parser = argparse.ArgumentParser(description="Top-Down Glycan Pipeline")
    parser.add_argument("--limit", type=int, default=0, help="Number of rows to process (0 for all)")
    args = parser.parse_args()
    
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    data_dir = os.path.join(base_dir, "data")
    reports_dir = os.path.join(base_dir, "reports")
    os.makedirs(reports_dir, exist_ok=True)
    
    raw_csv = os.path.join(data_dir, "Coconut.csv")
    sugar_csv = os.path.join(reports_dir, "Coconut_Sugar.csv")
    check_csv = os.path.join(reports_dir, "Coconut_Sugar_Check.csv")
    
    # Phase 1
    df = phase1_initialization_and_deduplication(raw_csv, sugar_csv, check_csv)
    
    if args.limit > 0:
        df = df.head(args.limit).copy()
        print(f"Limiting execution to {args.limit} rows for testing.")
        
    # Prepare to append new columns
    # New columns: Glycan_SMILES, Aglycan_SMILE_ALL, AminoAcid_SMILES, NUCLEOTIDES_SMILES, Family, Aglycan_Scaffold_IDs, murcko_framework, Sugar_Functional_Group
    new_columns = ['Glycan_SMILES', 'Aglycan_SMILE_ALL', 'AminoAcid_SMILES', 'NUCLEOTIDES_SMILES', 'Family', 'Aglycan_Scaffold_IDs', 'murcko_framework', 'Sugar_Functional_Group']
    for col in new_columns:
        if col not in df.columns:
            df[col] = None
            
    # Execute Phase 2: Core Matrix processing
    print("--- [Phase 2] Core Structure Cleavage ---")
    
    res_rows = []
    
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Processing Compounds"):
        row_dict = dict(row)
        smiles = str(row_dict.get('canonical_smiles', ''))
        
        # Reset new columns
        for col in new_columns:
            row_dict[col] = "NULL"
            
        if not smiles or smiles == 'nan':
            res_rows.append(row_dict)
            continue
            
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            res_rows.append(row_dict)
            continue
            
        # [Phase 2 Cleavage Logic]
        try:
            # Enforce Strict Topology checks & split
            aglycan_smi, glycan_smi, is_false_positive = sugar_utils.strict_cleavage(mol)
            
            if is_false_positive:
                row_dict['Glycan_SMILES'] = "FALSE_POSITIVE"
                row_dict['Aglycan_SMILE_ALL'] = "FALSE_POSITIVE"
                # We skip secondary fragments if it is a false positive
            else:
                row_dict['Aglycan_SMILE_ALL'] = aglycan_smi if aglycan_smi else "NULL"
                row_dict['Glycan_SMILES'] = glycan_smi if glycan_smi else "NULL"
                
                # --- [Phase 3] Secondary Fragments ---
            # Nucleotides: strict search in Glycan Part
            nucl_smi = "NULL"
            if glycan_smi:
                glycan_mol = Chem.MolFromSmiles(glycan_smi)
                if glycan_mol:
                    nucl_atoms = sugar_utils.find_nucleotides(glycan_mol)
                    if nucl_atoms:
                        nucl_smi = Chem.MolFragmentToSmiles(glycan_mol, atomsToUse=list(nucl_atoms), isomericSmiles=True)
            row_dict['NUCLEOTIDES_SMILES'] = nucl_smi
            
            # Amino Acids: strict search in Aglycan Part
            aa_smi = "NULL"
            if aglycan_smi:
                aglycan_mol = Chem.MolFromSmiles(aglycan_smi)
                if aglycan_mol:
                    peptides = sugar_utils.extract_amino_acids_and_peptides(aglycan_mol)
                    if peptides:
                        # Combine all extracted aminos
                        aa_smiles_list = [p['smiles'] for p in peptides if p['smiles']]
                        if aa_smiles_list:
                            aa_smi = "; ".join(aa_smiles_list)
            row_dict['AminoAcid_SMILES'] = aa_smi
            
        except Exception as e:
            print(f"Error on {smiles}: {e}")
            pass
            
        res_rows.append(row_dict)
        
    df_result = pd.DataFrame(res_rows)
    fp_count = len(df_result[df_result['Glycan_SMILES'] == 'FALSE_POSITIVE'])
    print(f"Phase 2 & 3 complete. {fp_count} false positives tagged.")
    
    csv_phase3 = os.path.join(reports_dir, "Coconut_Sugar_Phase3.csv")
    df_result.to_csv(csv_phase3, index=False, encoding='utf-8-sig')
    print(f"Intermediate result saved to {csv_phase3}")
    
    # --- [Phase 4] Taxonomy Imputation ---
    print("--- [Phase 4] Taxonomy Imputation ---")
    df_result, tax_imputed = taxonomy_online_resolver.fill_taxonomy_online(df_result)
    
    # --- [Phase 5] Scaffold & Sequence ---
    print("--- [Phase 5] Scaffold & Sequence ---")
    
    tqdm.pandas(desc="Calculating Scaffolds and Sequences")
    
    def process_phase5(row):
        # 1. Aglycan Scaffold Fingerprint & Framework
        aglycan_smiles = str(row.get('Aglycan_SMILE_ALL', ''))
        scaf_id = "NULL"
        murcko = "NULL"
        
        if aglycan_smiles and aglycan_smiles not in ("NULL", "FALSE_POSITIVE"):
            try:
                mol = Chem.MolFromSmiles(aglycan_smiles)
                if mol:
                    from rdkit.Chem import rdFingerprintGenerator
                    mg = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
                    fp = mg.GetFingerprint(mol)
                    scaf_id = fp.ToBase64()
                    
                    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
                    if scaffold.GetNumAtoms() > 0:
                        murcko = Chem.MolToSmiles(scaffold)
            except:
                pass
                
        # 2. Sequence Sugars
        glycan_smiles = str(row.get('Glycan_SMILES', ''))
        seq = "NULL"
        mods = "NULL"
        
        if glycan_smiles and glycan_smiles not in ("NULL", "FALSE_POSITIVE"):
            try:
                gmol = Chem.MolFromSmiles(glycan_smiles)
                if gmol:
                    seq_tmp, mods_tmp = sugar_sequence.generate_refined_sequence(gmol)
                    if seq_tmp: seq = seq_tmp
                    if mods_tmp: mods = mods_tmp
            except:
                pass
                
        return pd.Series([scaf_id, murcko, seq, mods])
        
    df_result[['Aglycan_Scaffold_IDs', 'murcko_framework', 'Sugar_Sequence', 'Sugar_Functional_Group']] = df_result.progress_apply(process_phase5, axis=1)
    
    csv_phase5 = os.path.join(reports_dir, "Coconut_Sugar_Phase5.csv")
    df_result.to_csv(csv_phase5, index=False, encoding='utf-8-sig')
    print(f"Intermediate result saved to {csv_phase5}")
    
    # --- [Phase 6] Chemical Classification Imputation & Export ---
    print("--- [Phase 6] Chemical Super Class Imputation ---")
    
    # Handle NaNs originally
    df_result['chemical_super_class'] = df_result['chemical_super_class'].fillna('No Classified')
    df_result.loc[df_result['chemical_super_class'].str.strip() == '', 'chemical_super_class'] = 'No Classified'
    
    df_result, class_imputed = taxonomy_online_resolver.fill_classification(df_result)
    
    # P6 二级分类填补: np_classifier_superclass
    # (P6 Secondary classification: fill np_classifier_superclass via same strategy)
    print("--- [Phase 6b] NP Classifier Superclass Imputation ---")
    df_result, np_imputed = taxonomy_online_resolver.fill_np_classifier(df_result)
    
    excel_classification = os.path.join(reports_dir, "Coconut_Classification.xlsx")
    print(f"Exporting classification sheets to {excel_classification}...")
    
    def save_to_excel(df_dict, filepath):
        if not df_dict:
            df_dict = {'Empty': pd.DataFrame()}
        with pd.ExcelWriter(filepath, engine='openpyxl') as writer:
            for sheet_name, sheet_df in df_dict.items():
                safe_name = str(sheet_name).replace('/', '_').replace('\\', '_').replace('[', '').replace(']', '')[:31]
                if not safe_name: safe_name = 'Sheet'
                sheet_df.to_excel(writer, sheet_name=safe_name, index=False)
                
    grouped_class = {name: group for name, group in df_result.groupby('chemical_super_class')}
    save_to_excel(grouped_class, excel_classification)
    
    csv_phase6 = os.path.join(reports_dir, "Coconut_Sugar_Phase6.csv")
    df_result.to_csv(csv_phase6, index=False, encoding='utf-8-sig')
    
    # --- [Phase 7] Visualization & Export ---
    print("--- [Phase 7] Visualization & Export ---")
    
    img_base_dir = os.path.join(base_dir, "images")
    if os.path.exists(img_base_dir):
        print("Cleaning up old image directory...")
        shutil.rmtree(img_base_dir)
    os.makedirs(img_base_dir, exist_ok=True)
    
    viz = StructureVisualizer()
    
    for idx, row in tqdm(df_result.iterrows(), total=len(df_result), desc="Drawing images"):
        smiles = str(row.get('canonical_smiles', ''))
        identifier = str(row.get('identifier', f'Row_{idx}'))
        super_class = str(row.get('chemical_super_class', 'Unknown'))
        
        # We will use 'Unknown' subclass if np_classifier_superclass isn't available
        sub_class = str(row.get('np_classifier_superclass', 'Unknown'))
        if not sub_class or sub_class == 'nan': sub_class = 'Unknown'
        
        if smiles and smiles != 'nan':
            # Check if this row is a False Positive
            is_fp = (str(row.get('Glycan_SMILES', '')) == "FALSE_POSITIVE")
            
            if is_fp:
                img_dir_path = os.path.join(img_base_dir, "TEST_False_Positives")
            else:
                safe_sup = "".join([c for c in super_class if c.isalnum() or c in (' ', '_', '-')]).strip()
                safe_sub = "".join([c for c in sub_class if c.isalnum() or c in (' ', '_', '-')]).strip()
                if not safe_sup: safe_sup = "No_Classified"
                if not safe_sub: safe_sub = "No_Classified"
                img_dir_path = os.path.join(img_base_dir, safe_sup, safe_sub)
                
            os.makedirs(img_dir_path, exist_ok=True)
            
            img_name = f"{identifier}.png"
            img_path = os.path.join(img_dir_path, img_name)
            
            try:
                # Rendering is handled internally by StructureVisualizer
                # Based on user direction, Aglycan=Blue, Sugar=Red, Sub=Yellow, Nuc/AA=LightGreen
                viz.analyze_glycolipid(smiles, img_path)
            except Exception as e:
                pass

    print("✅ All 7 phases of reconstruction pipeline complete! You can view outputs in the `reports/` and `images/` directory.")
    
    
if __name__ == "__main__":
    main()
