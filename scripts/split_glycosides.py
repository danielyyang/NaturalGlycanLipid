
import os
import sys
import pandas as pd
import logging
from rdkit import Chem
from tqdm import tqdm
import argparse
import shutil

# Dynamic Path Setup
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from lib import sugar_utils
from lib.visualizer import StructureVisualizer

# Config
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

INPUT_FILE = r"d:\Lipid Database\data\processed\Coconut_Sugars_Unique.xlsx"
OUTPUT_FILE = r"d:\Lipid Database\data\processed\Coconut_Sugar_Sort.csv"
IMAGE_BASE_DIR = r"d:\Lipid Database\images"

def cleanup_images():
    """Remove existing images to ensure fresh generation."""
    if os.path.exists(IMAGE_BASE_DIR):
        logger.info(f"Cleaning up image directory: {IMAGE_BASE_DIR}")
        # Iterating to avoid deleting the base dir if it's a mount or something, 
        # but here it's just a folder. 
        # Safety: only delete subfolders relevant to sheets? 
        # For now, user said "delete coverage", so we wipe the folders.
        for item in os.listdir(IMAGE_BASE_DIR):
            item_path = os.path.join(IMAGE_BASE_DIR, item)
            if os.path.isdir(item_path):
                try:
                    shutil.rmtree(item_path)
                    logger.info(f"Removed {item_path}")
                except Exception as e:
                    logger.warning(f"Could not remove {item_path}: {e}")

def process_file(limit: int = 0, mode: str = 'all'):
    """
    mode: 'all' (default), 'extract_only' (CSV), 'viz_only' (Images from CSV)
    """
    logger.info(f"Starting process in mode: {mode}")

    if mode in ['all', 'extract_only']:
        # --- PHASE 1: Data Extraction (Excel to CSV) ---
        logger.info(f"Loading data from {INPUT_FILE}...")
        if not os.path.exists(INPUT_FILE):
            logger.error("Input file not found!")
            return

        xl = pd.ExcelFile(INPUT_FILE)
        
        # We only clean images if we are ALSO generating them in this run, OR if explicitly requested?
        # User said "Organize CSV then generate images". 
        # If running extract_only, we don't touch images. 
        # If running all, we clean.
        if mode == 'all':
             cleanup_images()

        all_data = []

        for sheet_name in xl.sheet_names:
            logger.info(f"Processing sheet for Data: {sheet_name}")
            df = xl.parse(sheet_name)
            
            # Add Sheet Source
            df['Sheet_Source'] = sheet_name
            
            # Find SMILES Column
            smiles_col = None
            for col in df.columns:
                if str(col).lower() in ['smiles', 'canonical_smiles']:
                    smiles_col = col
                    break
            
            if not smiles_col:
                logger.warning(f"No SMILES column found in {sheet_name}, skipping.")
                continue
                
            rows_to_process = limit if limit > 0 else len(df)
            
            # Lists for new columns
            aglycone_list = []
            glycan_list = []
            is_sugar_list = []
            
            for idx, row in tqdm(df.head(rows_to_process).iterrows(), total=rows_to_process, desc=f"Extracting {sheet_name}"):
                smiles = row[smiles_col]
                agly_smi = ""
                glyc_smi = ""
                is_sugar = False
                
                if isinstance(smiles, str) and smiles:
                    try:
                        mol = Chem.MolFromSmiles(smiles)
                        if mol:
                            sugar_units, _ = sugar_utils.get_sugar_units(mol)
                            is_sugar = len(sugar_units) > 0
                            if is_sugar:
                                agly_smi, glyc_smi = sugar_utils.get_split_smiles(mol)
                    except:
                        pass # Logging inside loop slows things down too much for big data
                
                aglycone_list.append(agly_smi)
                glycan_list.append(glyc_smi)
                is_sugar_list.append(is_sugar)
                
            # Create processed DF partition
            df_processed = df.head(rows_to_process).copy()
            df_processed['Aglycone_SMILES'] = aglycone_list
            df_processed['Glycan_SMILES'] = glycan_list
            df_processed['Is_Glycoside'] = is_sugar_list
            
            all_data.append(df_processed)

        # Combine and Save CSV
        if all_data:
            logger.info("Combining all sheets...")
            final_df = pd.concat(all_data, ignore_index=True)
            
            # Reorder columns
            cols = list(final_df.columns)
            target_col = 'canonical_smiles'
            if target_col not in cols:
                target_col = [c for c in cols if 'smiles' in str(c).lower()][0] if [c for c in cols if 'smiles' in str(c).lower()] else cols[0]
            
            if target_col in cols:
                target_idx = cols.index(target_col)
                new_cols = ['Glycan_SMILES', 'Aglycone_SMILES', 'Sheet_Source', 'Is_Glycoside']
                for nc in new_cols:
                    if nc in cols: cols.remove(nc)
                final_cols = cols[:target_idx+1] + new_cols + cols[target_idx+1:]
                final_df = final_df[final_cols]
                
            logger.info(f"Saving {len(final_df)} rows to {OUTPUT_FILE}...")
            final_df.to_csv(OUTPUT_FILE, index=False)
            logger.info("CSV Generation Done.")
        else:
            logger.warning("No data processed.")

    if mode in ['all', 'viz_only']:
        # --- PHASE 2: Visualization (CSV to Images) ---
        logger.info(f"Starting Visualization Phase...")
        
        if not os.path.exists(OUTPUT_FILE):
             logger.error(f"CSV file {OUTPUT_FILE} not found! Cannot generate images.")
             return
             
        # If viz_only, we might want to clean up too? 
        # User said "Generate all images". Safe to clean up if we are doing a full run.
        if mode == 'viz_only':
            cleanup_images()
            
        df = pd.read_csv(OUTPUT_FILE)
        viz = StructureVisualizer()
        
        # We need to know which column has SMILES.
        smiles_col = None
        for col in df.columns:
            if str(col).lower() in ['smiles', 'canonical_smiles']:
                smiles_col = col
                break
        
        if not smiles_col:
             logger.error("No SMILES column in CSV.")
             return

        # Iterate only rows that are Glycosides
        df_glycosides = df[df['Is_Glycoside'] == True]
        total_glycosides = len(df_glycosides)
        logger.info(f"Found {total_glycosides} glycosides to visualize.")
        
        for idx, row in tqdm(df_glycosides.iterrows(), total=total_glycosides, desc="Generating Images"):
            smiles = row[smiles_col]
            sheet_name = row['Sheet_Source']
            
            # Ensure dir exists
            sheet_img_dir = os.path.join(IMAGE_BASE_DIR, sheet_name)
            os.makedirs(sheet_img_dir, exist_ok=True)
            
            # Filename
            safe_sheet_name = "".join([c for c in str(sheet_name) if c.isalnum() or c in (' ','_','-')]).strip()
            img_name = f"struct_{safe_sheet_name}_{idx}.png"
            img_path = os.path.join(sheet_img_dir, img_name)
            
            try:
                # We can pass the already split parts to viz if we wanted, but viz logic does it fresh from SMILES.
                # Visualization from SMILES is safer for layout.
                viz.analyze_glycolipid(smiles, img_path)
            except Exception as e:
                logger.error(f"Failed to visualize {idx}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=0, help="Limit rows per sheet for testing")
    parser.add_argument("--mode", type=str, default='all', choices=['all', 'extract_only', 'viz_only'], help="Operation mode")
    args = parser.parse_args()
    
    process_file(args.limit, args.mode)
