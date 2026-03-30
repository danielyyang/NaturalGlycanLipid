import pandas as pd
import requests
import time
import urllib.parse
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

def get_morgan_fp(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    except:
        pass
    return None

def fill_classification(df):
    """
    Fills 'No Classified' rows based on:
    1. InChIKey (first 14 chars) match with a classified row.
    2. Aglycan_SMILE Morgan Fingerprint Tanimoto similarity >= 0.8 with a classified row.
    Returns (df, imputed_cells_set)
    """
    imputed_cells = set()
    classified_df = df[df['chemical_super_class'] != 'No Classified'].copy()
    unclassified_df = df[df['chemical_super_class'] == 'No Classified'].copy()

    if unclassified_df.empty or classified_df.empty:
        return df, imputed_cells
        
    print(f"Attempting to classify {len(unclassified_df)} unclassified rows...")
    
    # Pre-compute InChIKey blocks and fingerprints for classified rows
    classified_df['inchikey_block1'] = classified_df['standard_inchi_key'].astype(str).str[:14]
    
    # Only compute FP for the first Aglycan part
    def extract_primary_aglycan(s):
        if pd.isna(s) or not str(s).strip(): return ""
        return str(s).split('|')[0].strip()
        
    classified_df['primary_aglycan'] = classified_df['Aglycan_SMILE_ALL'].apply(extract_primary_aglycan)
    tqdm.pandas(desc="Computing Classified FPs")
    classified_fps = classified_df['primary_aglycan'].progress_apply(get_morgan_fp)
    classified_df['fp'] = classified_fps
    valid_classified = classified_df.dropna(subset=['fp']).copy()

    changes = 0
    tqdm.pandas(desc="Inferring Classes")
    
    for idx, row in tqdm(unclassified_df.iterrows(), total=len(unclassified_df), desc="Classifying"):
        # Strategy 1: InChIKey Block 1 Match
        ikey = str(row.get('standard_inchi_key', ''))[:14]
        match = classified_df[classified_df['inchikey_block1'] == ikey]
        if not match.empty:
            inferred_class = match.iloc[0]['chemical_super_class']
            df.at[idx, 'chemical_super_class'] = inferred_class
            imputed_cells.add((idx, 'chemical_super_class'))
            changes += 1
            continue
            
        # Strategy 2: Tanimoto Similarity on Primary Aglycan > 0.8
        primary_ag = extract_primary_aglycan(row.get('Aglycan_SMILE_ALL'))
        if primary_ag:
            fp1 = get_morgan_fp(primary_ag)
            if fp1:
                best_sim = 0
                best_class = None
                for c_idx, c_row in valid_classified.iterrows():
                    fp2 = c_row['fp']
                    sim = DataStructs.TanimotoSimilarity(fp1, fp2)
                    if sim > best_sim:
                        best_sim = sim
                        best_class = c_row['chemical_super_class']
                        
                if best_sim >= 0.8 and best_class:
                    df.at[idx, 'chemical_super_class'] = best_class
                    imputed_cells.add((idx, 'chemical_super_class'))
                    changes += 1

    print(f"Successfully inferred classes for {changes} rows.")
    return df, imputed_cells


def fill_np_classifier(df):
    """
    填补空缺的 np_classifier_superclass 列。
    策略与 fill_classification 相同：InChIKey Block-1 匹配 + Aglycan Morgan FP Tanimoto ≥ 0.8。
    Fills empty 'np_classifier_superclass' using the same strategy as fill_classification:
    1. InChIKey Block-1 match with a classified row.
    2. Aglycan Morgan FP Tanimoto similarity >= 0.8.
    Returns (df, imputed_cells_set)
    """
    imputed_cells = set()

    if 'np_classifier_superclass' not in df.columns:
        print("Column 'np_classifier_superclass' not found. Skipping.")
        return df, imputed_cells

    # 标准化空值标记 (Standardize empty markers)
    df['np_classifier_superclass'] = df['np_classifier_superclass'].fillna('No Classified')
    df.loc[df['np_classifier_superclass'].str.strip() == '', 'np_classifier_superclass'] = 'No Classified'

    classified_df = df[df['np_classifier_superclass'] != 'No Classified'].copy()
    unclassified_df = df[df['np_classifier_superclass'] == 'No Classified'].copy()

    if unclassified_df.empty or classified_df.empty:
        return df, imputed_cells

    print(f"Attempting to fill np_classifier_superclass for {len(unclassified_df)} rows...")

    # 预计算 InChIKey Block-1 和分子指纹 (Pre-compute InChIKey blocks and fingerprints)
    classified_df['inchikey_block1'] = classified_df['standard_inchi_key'].astype(str).str[:14]

    def extractPrimaryAglycan(s):
        if pd.isna(s) or not str(s).strip(): return ""
        return str(s).split('|')[0].strip()

    classified_df['primary_aglycan'] = classified_df['Aglycan_SMILE_ALL'].apply(extractPrimaryAglycan)
    tqdm.pandas(desc="Computing NP Classifier FPs")
    classifiedFps = classified_df['primary_aglycan'].progress_apply(get_morgan_fp)
    classified_df['fp'] = classifiedFps
    validClassified = classified_df.dropna(subset=['fp']).copy()

    changes = 0
    for idx, row in tqdm(unclassified_df.iterrows(), total=len(unclassified_df), desc="Filling NP Classifier"):
        # 策略 1: InChIKey Block-1 匹配 (Strategy 1: InChIKey Block-1 Match)
        ikey = str(row.get('standard_inchi_key', ''))[:14]
        match = classified_df[classified_df['inchikey_block1'] == ikey]
        if not match.empty:
            inferredClass = match.iloc[0]['np_classifier_superclass']
            df.at[idx, 'np_classifier_superclass'] = inferredClass
            imputed_cells.add((idx, 'np_classifier_superclass'))
            changes += 1
            continue

        # 策略 2: Tanimoto 相似度 >= 0.8 (Strategy 2: Tanimoto Similarity >= 0.8)
        primaryAg = extractPrimaryAglycan(row.get('Aglycan_SMILE_ALL'))
        if primaryAg:
            fp1 = get_morgan_fp(primaryAg)
            if fp1:
                bestSim = 0
                bestClass = None
                for cIdx, cRow in validClassified.iterrows():
                    fp2 = cRow['fp']
                    sim = DataStructs.TanimotoSimilarity(fp1, fp2)
                    if sim > bestSim:
                        bestSim = sim
                        bestClass = cRow['np_classifier_superclass']

                if bestSim >= 0.8 and bestClass:
                    df.at[idx, 'np_classifier_superclass'] = bestClass
                    imputed_cells.add((idx, 'np_classifier_superclass'))
                    changes += 1

    print(f"Successfully inferred np_classifier_superclass for {changes} rows.")
    return df, imputed_cells

# --- Online Taxonomy Imputation ---
TAXONOMY_DICT = {
    "Arabidopsis thaliana": "Brassicaceae",
    "Oryza sativa": "Poaceae",
    "Zea mays": "Poaceae",
    "Saccharomyces cerevisiae": "Saccharomycetaceae",
    "Homo sapiens": "Hominidae",
    "Mus musculus": "Muridae",
    "Escherichia coli": "Enterobacteriaceae",
    "Bacillus subtilis": "Bacillaceae",
    "Panax ginseng": "Araliaceae",
    "Glycyrrhiza glabra": "Fabaceae",
    "Astragalus membranaceus": "Fabaceae",
    "Staphylococcus aureus": "Staphylococcaceae",
    "Pseudomonas aeruginosa": "Pseudomonadaceae",
}

def get_organism_from_pubchem(iupac_name):
    if not iupac_name or pd.isna(iupac_name):
        return None
    try:
        # 1. Translate Name to CID
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{urllib.parse.quote(str(iupac_name))}/cids/JSON"
        res = requests.get(url, timeout=10)
        if res.status_code == 200:
            cids = res.json().get('IdentifierList', {}).get('CID', [])
            if cids:
                cid = cids[0]
                # 2. Get Taxonomy from PUG View
                tax_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
                tax_res = requests.get(tax_url, timeout=10)
                if tax_res.status_code == 200:
                    data = tax_res.json()
                    record = data.get('Record', {})
                    for sec in record.get('Section', []):
                        if sec.get('TOCHeading') in ('Taxonomy', 'Associated Disorders and Diseases', 'Biomolecular Interactions and Pathways'):
                            for sub_sec in sec.get('Section', []) + [sec]:
                                for info in sub_sec.get('Information', []):
                                    if 'Name' in info and 'Value' in info:
                                        name = info['Name']
                                        if 'Organism' in name or 'Species' in name or 'Taxonomy' in name:
                                            try:
                                                val = info['Value']['StringWithMarkup'][0]['String']
                                                return val
                                            except:
                                                pass
                                        elif sec.get('TOCHeading') == 'Taxonomy':
                                            try:
                                                val = info['Value']['StringWithMarkup'][0]['String']
                                                return val
                                            except:
                                                pass
        time.sleep(0.3)
    except Exception as e:
        pass
    return None

def get_organism_from_coconut(inchikey):
    """
    Level 1: Query native COCONUT database by InChIKey.
    """
    if not inchikey or pd.isna(inchikey):
        return None
    try:
        # COCONUT v1/v2 API endpoints can be volatile but 'search' is often available.
        url = f"https://coconut.naturalproducts.net/api/search/exact?inchikey={urllib.parse.quote(str(inchikey))}"
        res = requests.get(url, timeout=10)
        if res.status_code == 200:
            data = res.json()
            # Handle potential list or object responses
            if isinstance(data, list) and len(data) > 0:
                taxa = data[0].get('taxonomies', [])
                if taxa:
                    for t in taxa:
                        if 'species' in t and t['species']: return t['species']
            elif isinstance(data, dict):
                taxa = data.get('taxonomies', [])
                if taxa:
                    for t in taxa:
                        if 'species' in t and t['species']: return t['species']
        time.sleep(0.2)
    except:
        pass
    return None

def get_organism_from_wikidata(inchikey):
    """
    Level 2: Query Wikidata Knowledge Graph (LOTUS Database records) via SPARQL.
    wdt:P235 = InChIKey, wdt:P703 = found in taxon, wdt:P225 = taxon name.
    """
    if not inchikey or pd.isna(inchikey):
        return None
    try:
        query = f'''
        SELECT ?taxonName WHERE {{
          ?compound wdt:P235 "{inchikey}" .
          ?compound wdt:P703 ?taxon .
          ?taxon wdt:P225 ?taxonName .
        }} LIMIT 1
        '''
        url = "https://query.wikidata.org/sparql"
        params = {'query': query, 'format': 'json'}
        headers = {'User-Agent': 'GlycoLipidInsight-Bot/1.0 (python-requests)'}
        res = requests.get(url, params=params, headers=headers, timeout=10)
        if res.status_code == 200:
            data = res.json()
            bindings = data.get('results', {}).get('bindings', [])
            if bindings:
                return bindings[0]['taxonName']['value']
        time.sleep(0.3)
    except:
        pass
    return None

def get_organism_from_wikipedia(query_name):
    """
    Fallback method: Searches Wikipedia for the compound name and
    attempts to extract the first italicized binomial name (e.g. <i>Panax ginseng</i>)
    from the article's summary.
    """
    if not query_name or pd.isna(query_name):
        return None
    try:
        search_url = f"https://en.wikipedia.org/w/api.php?action=query&list=search&srsearch={urllib.parse.quote(str(query_name))}&utf8=&format=json"
        res = requests.get(search_url, timeout=10)
        if res.status_code == 200:
            data = res.json()
            search_results = data.get('query', {}).get('search', [])
            if search_results:
                page_title = search_results[0]['title']
                
                ex_url = f"https://en.wikipedia.org/w/api.php?action=query&prop=extracts&titles={urllib.parse.quote(page_title)}&format=json&exintro=1"
                ex_res = requests.get(ex_url, timeout=10)
                if ex_res.status_code == 200:
                    pages = ex_res.json().get('query', {}).get('pages', {})
                    if pages:
                        for page_id, info in pages.items():
                            extract = info.get('extract', '')
                            # Extract italicized species names like <i>Panax pseudo-ginseng</i> or <i>M. musculus</i>
                            import re
                            matches = re.findall(r'<i>([A-Z][a-z]+ [a-z]+)</i>', extract)
                            if matches:
                                return matches[0]
        time.sleep(0.3)
    except:
        pass
    return None

def get_family_from_gbif(organism):
    if not organism or pd.isna(organism):
        return None
        
    import re
    # Split by common separators if the dataset provides multiple synonyms
    parts = re.split(r'[|;,]', str(organism))
    
    for part in parts:
        clean_org = part.strip()
        if not clean_org:
            continue
            
        try:
            url = f"https://api.gbif.org/v1/species/match?name={urllib.parse.quote(clean_org)}"
            res = requests.get(url, timeout=10)
            if res.status_code == 200:
                data = res.json()
                family = data.get('family', None)
                if family:
                    return family
            time.sleep(0.1)
        except Exception as e:
            pass
            
        # Fallback to Wikipedia to search for family directly
        try:
            search_url = f"https://en.wikipedia.org/w/api.php?action=query&list=search&srsearch={urllib.parse.quote(clean_org + ' family')}&utf8=&format=json"
            res = requests.get(search_url, timeout=10)
            if res.status_code == 200:
                data = res.json()
                search_results = data.get('query', {}).get('search', [])
                if search_results:
                    page_title = search_results[0]['title']
                    ex_url = f"https://en.wikipedia.org/w/api.php?action=query&prop=extracts&titles={urllib.parse.quote(page_title)}&format=json&exintro=1"
                    ex_res = requests.get(ex_url, timeout=10)
                    if ex_res.status_code == 200:
                        pages = ex_res.json().get('query', {}).get('pages', {})
                        if pages:
                            for page_id, info in pages.items():
                                extract = info.get('extract', '')
                                # Mostly plant families end in -aceae
                                matches = re.findall(r'\b([A-Z][a-z]+aceae)\b', extract)
                                if matches:
                                    return matches[0]
            time.sleep(0.3)
        except:
            pass
            
    return None

def fill_taxonomy_online(df):
    """
    Fills 'organisms' and 'Family' using online databases.
    - organisms: Fetched from PubChem Taxonomy via iupac_name
    - Family: Fetched from GBIF Species match via organisms
    Flags rows that were imputed so they can be colored red in Excel.
    """
    if 'Family' not in df.columns:
        if 'organisms' in df.columns:
            org_idx = df.columns.get_loc('organisms')
            df.insert(org_idx + 1, 'Family', "")
        else:
            df['organisms'] = ""
            df['Family'] = ""
            
    imputed_cells = set()
    changes_org = 0
    changes_fam = 0
    
    print("Applying Online Taxonomy Resolution (PubChem & GBIF)...")
    
    pubchem_cache = {}
    gbif_cache = {}
    
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Fetching Taxonomy from API"):
        org = str(row.get('organisms', '')).strip()
        fam = str(row.get('Family', '')).strip()
        compound_name = str(row.get('name', '')).strip()
        iupac = str(row.get('iupac_name', '')).strip()
        inchikey = str(row.get('standard_inchi_key', '')).strip()
        
        query_name = compound_name if compound_name and compound_name != 'nan' else iupac
        
        # 1. Fill Organism
        if not org or org == 'nan' or org == 'NULL':
            inferred_org = None
            
            # Tier 1: COCONUT exact InChIKey match
            if not inferred_org and inchikey and inchikey != 'nan':
                inferred_org = get_organism_from_coconut(inchikey)
                
            # Tier 2: LOTUS (Wikidata SPARQL) exact InChIKey match
            if not inferred_org and inchikey and inchikey != 'nan':
                inferred_org = get_organism_from_wikidata(inchikey)
            
            # Tier 3: PubChem IUPAC text match
            if not inferred_org and query_name and query_name != 'nan':
                if query_name in pubchem_cache:
                    inferred_org = pubchem_cache[query_name]
                else:
                    inferred_org = get_organism_from_pubchem(query_name)
                    pubchem_cache[query_name] = inferred_org

            # Tier 4: Wikipedia NLP text match
            if not inferred_org and query_name and query_name != 'nan':
                inferred_org = get_organism_from_wikipedia(query_name)
                
            if inferred_org:
                df.at[idx, 'organisms'] = inferred_org
                imputed_cells.add((idx, 'organisms'))
                org = inferred_org # update org for next step
                changes_org += 1
                
        # 2. Fill Family
        if org and org != 'nan' and (not fam or fam == 'nan'):
            if org in TAXONOMY_DICT:
                df.at[idx, 'Family'] = TAXONOMY_DICT[org]
                imputed_cells.add((idx, 'Family'))
                changes_fam += 1
            else:
                if org in gbif_cache:
                    inferred_fam = gbif_cache[org]
                else:
                    inferred_fam = get_family_from_gbif(org)
                    gbif_cache[org] = inferred_fam
                    
                    if inferred_fam:
                        df.at[idx, 'Family'] = inferred_fam
                        imputed_cells.add((idx, 'Family'))
                        changes_fam += 1
                        
    # Replace any remaining blanks with 'Not Result'
    df['organisms'] = df['organisms'].replace(['', 'nan', None], 'NULL').fillna('NULL')
    df['Family'] = df['Family'].replace(['', 'nan', None], 'NULL').fillna('NULL')
    
    # Standardize 'Not Result' flag in imputed cells if we modified them just now?
    # No need to flag 'Not Result' as a successful imputation, leave it as is so it doesn't color green/red incorrectly.

    print(f"Successfully inferred {changes_org} organisms and {changes_fam} Family taxonomy rows via API.")
    return df, imputed_cells
