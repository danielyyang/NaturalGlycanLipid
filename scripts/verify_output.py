import pandas as pd
import os
import logging
import sys

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Output file path
OUTPUT_FILE = r"d:\Lipid Database\data\processed\Coconut_Sugar_Sort.csv"
IMAGE_BASE_DIR = r"d:\Lipid Database\images"

# Set stdout encoding to utf-8 for Windows console compatibility
sys.stdout.reconfigure(encoding='utf-8')

def verify_output():
    logger.info(f"Verifying output file: {OUTPUT_FILE}")
    
    if not os.path.exists(OUTPUT_FILE):
        logger.error("Output file not found!")
        return

    try:
        df = pd.read_csv(OUTPUT_FILE)
        logger.info(f"Successfully loaded CSV. Total rows: {len(df)}")
        
        # Check required columns
        required_cols = ['Aglycone_SMILES', 'Glycan_SMILES', 'Is_Glycoside', 'Sheet_Source']
        missing_cols = [col for col in required_cols if col not in df.columns]
        
        if missing_cols:
            logger.error(f"Missing columns: {missing_cols}")
        else:
            logger.info("All required columns are present.")
            
        # Check for non-empty split results
        # Assuming 'Is_Glycoside' is boolean
        glycosides = df[df['Is_Glycoside'] == True]
        logger.info(f"Found {len(glycosides)} glycosides.")
        
        if not glycosides.empty:
            sample = glycosides.iloc[0]
            logger.info(f"Sample Glycoside from {sample['Sheet_Source']}:")
            logger.info(f"  Precursor: {sample.get('canonical_smiles', 'N/A')}")
            logger.info(f"  Aglycone : {sample.get('Aglycone_SMILES', 'N/A')}")
            logger.info(f"  Glycan   : {sample.get('Glycan_SMILES', 'N/A')}")
            
            # Check for images directory
            sheet_name = sample.get('Sheet_Source')
            if sheet_name:
                sheet_dir = os.path.join(IMAGE_BASE_DIR, sheet_name)
                if os.path.exists(sheet_dir):
                    images = os.listdir(sheet_dir)
                    logger.info(f"Image directory for {sheet_name} exists and contains {len(images)} files.")
                else:
                    logger.warning(f"Image directory for {sheet_name} NOT found.")

    except Exception as e:
        logger.error(f"Error verifying file: {e}")

if __name__ == "__main__":
    verify_output()
