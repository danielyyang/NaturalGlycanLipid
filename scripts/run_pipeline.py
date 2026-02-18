import os
import sys
import argparse
import logging
import pandas as pd
from typing import List, Dict
from tqdm import tqdm

# Ensure we can import from lib
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from lib.aggregator import DataAggregator
from lib.profiler import PropertyProfiler
from lib.visualizer import StructureVisualizer
from lib.researcher import LiteratureResearcher
from lib import sugar_utils

# Configure Logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def run_pipeline(keyword: str, source: str, limit: int):
    """
    Run the full GlycoLipid-Insight pipeline.
    """
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    output_dir = os.path.join(project_root, 'reports')
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize Engines
    aggregator = DataAggregator()
    profiler = PropertyProfiler()
    visualizer = StructureVisualizer()
    researcher = LiteratureResearcher(email="example@email.com")
    
    logger.info(f"--- Starting Pipeline: Source={source}, Keyword='{keyword}', Limit={limit} ---")
    
    # 1. Fetch Data
    raw_data = []
    if source.lower() == 'local':
        # If local, keyword serves as file path
        local_path = keyword
        # Try to resolve relative path if not absolute
        if not os.path.isabs(local_path):
            local_path = os.path.join(project_root, local_path)
            
        if not os.path.exists(local_path):
            # Fallback check standard location
            default_path = os.path.join(project_root, 'data', 'raw', 'coconut.csv.csv')
            if os.path.exists(default_path):
                logger.info(f"Target file not found, using default: {default_path}")
                local_path = default_path
            else:
                logger.error(f"Local file not found: {keyword} or {local_path}")
                return

        # Determine limits (0 or negative means no limit)
        fetch_limit = limit if limit and limit > 0 else None
        raw_data = aggregator.fetch_from_local_file(local_path, limit=fetch_limit)
        
    elif source.lower() == 'pubchem':
        compounds = aggregator.fetch_pubchem_data(keyword, limit=limit)
        raw_data.extend(compounds)
    elif source.lower() == 'coconut':
        compounds = aggregator.fetch_coconut_data(keyword)
        raw_data.extend(compounds)
    elif source.lower() == 'knapsack':
        compounds = aggregator.fetch_knapsack_data(keyword)
        raw_data.extend(compounds)
        
    if not raw_data:
        logger.warning("No data found. Exiting.")
        return

    # 2. Unify Data
    df = aggregator.unify_results(raw_data)
    logger.info(f"Unified {len(df)} records.")
    
    # 3. Sugar Validation & Classification
    valid_glycolipids = []
    non_sugar_compounds = []
    
    logger.info("Validating sugar structures...")
    # Using tqdm for progress bar
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Analyzing Structures"):
        smiles = row.get('SMILES')
        
        # Debugging: Log if SMILES is missing
        if not smiles or not isinstance(smiles, str):
            # Try to recover from other columns if possible or just log
            # logger.warning(f"Row {idx} has invalid SMILES: {smiles}")
            non_sugar_compounds.append(row.to_dict())
            continue
            
        is_sugar, reason = sugar_utils.validate_sugar_structure(smiles)
        
        # Create a mutable copy of the row to add analysis results
        row_dict = row.to_dict()
        row_dict['Is_Sugar'] = is_sugar
        row_dict['Validation_Message'] = reason
        
        if is_sugar:
            valid_glycolipids.append(row_dict)
        else:
            # User Feedback: Even if "non-sugar", user might want them.
            # But we categorize them separatedly.
            non_sugar_compounds.append(row_dict)
            
    df_valid = pd.DataFrame(valid_glycolipids)
    df_invalid = pd.DataFrame(non_sugar_compounds)
    
    logger.info(f"Found {len(df_valid)} valid glycolipids and {len(df_invalid)} non-sugar compounds.")
    
    # 4. Enrichment (Only for Valid Glycolipids)
    if not df_valid.empty:
        logger.info("Enriching valid compounds with PubChem properties...")
        tqdm.pandas(desc="Enriching Properties")
        
        enriched_rows = []
        for idx, row in tqdm(df_valid.iterrows(), total=len(df_valid), desc="Fetching Props"):
            inchikey = row.get('InChIKey')
            if not inchikey and row.get('SMILES'):
                # Try to generate InChIKey if missing? (Skipped for now)
                pass
            
            if inchikey:
                props = profiler.get_pubchem_props(inchikey)
                if props:
                    # Merge properties
                    for k, v in props.items():
                        if k not in row: # Don't overwrite existing
                            row[k] = v
            enriched_rows.append(row)
        
        df_valid = pd.DataFrame(enriched_rows)

    # 5. Export to Excel
    excel_path = os.path.join(output_dir, 'unified_data.xlsx')
    logger.info(f"Saving reports to {excel_path}...")
    
    with pd.ExcelWriter(excel_path, engine='openpyxl') as writer:
        if not df_valid.empty:
            df_valid.to_excel(writer, sheet_name='Valid_Glycolipids', index=False)
        else:
            # Create empty sheet if no valid found to avoid error?
            pass
            
        if not df_invalid.empty:
            df_invalid.to_excel(writer, sheet_name='Non_Sugar_Compounds', index=False)
            
    # 6. Visualization & Image Embedding (Last Step)
    # Determine output folder based on keyword (filename)
    if source.lower() == 'local':
        base_name = os.path.splitext(os.path.basename(local_path))[0]
        image_output_dir = os.path.join(project_root, 'images', base_name)
    else:
        image_output_dir = os.path.join(project_root, 'images', 'batch')

    # Clean existing directory to ensure fresh images
    if os.path.exists(image_output_dir):
        import shutil
        try:
            shutil.rmtree(image_output_dir)
            logger.info(f"Cleaned existing image directory: {image_output_dir}")
        except Exception as e:
            logger.warning(f"Could not clean image directory: {e}")

    if not df_valid.empty:
        logger.info(f"Generating and embedding structures into {image_output_dir}...")
        visualizer.batch_process_from_file(excel_path, image_output_dir)
        
    # Also visualize "invalid" ones if user wants? 
    # Current requirement: "visualize ALL data files". 
    # visualizer.batch_process_from_file processes ALL sheets in the excel.
    # So if we saved 'Non_Sugar_Compounds', it will process that too IF it has SMILES column.
    
    logger.info("Pipeline Complete!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GlycoLipid-Insight Pipeline")
    parser.add_argument("--keyword", type=str, default="Ginsenoside", help="Search keyword or File Path for 'local'")
    parser.add_argument("--source", type=str, default="pubchem", choices=["pubchem", "coconut", "knapsack", "local"], help="Data source")
    parser.add_argument("--limit", type=int, default=0, help="Max results to fetch (0 for all)")
    
    args = parser.parse_args()
    
    run_pipeline(args.keyword, args.source, args.limit)
