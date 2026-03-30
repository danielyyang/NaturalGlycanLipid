"""
PubChem Authoritative SMILES Fetcher & Dictionary Patcher
Fetch IsomericSMILES from PubChem REST API for all key sugars,
canonicalize them with RDKit, then patch glycan_reference_library.py.
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from rdkit import Chem
import urllib.request

# All sugars to fetch from PubChem, mapped to PubChem compound names
PUBCHEM_SUGARS = {
    # Hexoses
    "D-Glc": "D-glucose",
    "D-Gal": "D-galactose",
    "D-Man": "D-mannose",
    "D-All": "D-allose",
    "D-Tal": "D-talose",
    "D-Alt": "D-altrose",
    "D-Gul": "D-gulose",
    "D-Ido": "D-idose",
    # L-Hexoses
    "L-Glc": "L-glucose",
    "L-Gal": "L-galactose",
    "L-Man": "L-mannose",
    # Pentoses
    "D-Xyl": "D-xylose",
    "D-Rib": "D-ribose",
    "D-Lyx": "D-lyxose",
    "L-Ara": "L-arabinose",
    "D-Ara": "D-arabinose",
    # Deoxy sugars
    "L-Rha": "L-rhamnose",
    "D-Rha": "D-rhamnose",
    "L-Fuc": "L-fucose",
    "D-Fuc": "D-fucose",
    "D-Qui": "D-quinovose",
    # Amino sugars
    "D-GlcNAc": "N-acetylglucosamine",
    "D-GalNAc": "N-acetylgalactosamine",
    "D-ManNAc": "N-acetylmannosamine",
    "D-GlcN": "D-glucosamine",
    "D-GalN": "D-galactosamine",
    # Uronic acids
    "D-GlcA": "D-glucuronic acid",
    "D-GalA": "D-galacturonic acid",
    "L-IdoA": "L-iduronic acid",
}

results = {}
errors = []

for name, pubchemName in PUBCHEM_SUGARS.items():
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{pubchemName.replace(' ', '%20')}/property/IsomericSMILES,InChI/JSON"
    try:
        req = urllib.request.Request(url)
        req.add_header('Accept', 'application/json')
        with urllib.request.urlopen(req, timeout=15) as resp:
            data = json.loads(resp.read().decode('utf-8'))
        
        props = data["PropertyTable"]["Properties"][0]
        pubchemSmi = props["SMILES"]
        pubchemInchi = props.get("InChI", "")
        
        # Canonicalize with RDKit
        mol = Chem.MolFromSmiles(pubchemSmi)
        if mol is None:
            errors.append(f"{name}: RDKit parse failed for '{pubchemSmi}'")
            continue
        
        canonSmi = Chem.MolToSmiles(mol, isomericSmiles=True)
        rdkitInchi = Chem.MolToInchi(mol)
        
        # Extract CIP from InChI
        results[name] = {
            "pubchem_smiles": pubchemSmi,
            "canonical_smiles": canonSmi,
            "pubchem_inchi": pubchemInchi,
            "rdkit_inchi": rdkitInchi,
        }
        print(f"  {name:12s}: {canonSmi}")
        
    except Exception as e:
        errors.append(f"{name}: {e}")
    
    time.sleep(0.3)

print(f"\nFetched {len(results)}/{len(PUBCHEM_SUGARS)} sugars from PubChem")
if errors:
    print(f"\nErrors ({len(errors)}):")
    for err in errors:
        print(f"  {err}")

# Compare with current dictionary
from lib.glycan_reference_library import RAW_MONOSACCHARIDE_SMILES

print("\n" + "=" * 100)
print("COMPARISON: Current Dictionary vs PubChem Authority")
print("=" * 100)

mismatches = []
for name, info in sorted(results.items()):
    ourKey = (name, "a")
    ourSmi = RAW_MONOSACCHARIDE_SMILES.get(ourKey, None)
    
    if ourSmi is None:
        print(f"{name:12s} MISSING in our dictionary!")
        mismatches.append((name, "MISSING", info["canonical_smiles"]))
        continue
    
    ourMol = Chem.MolFromSmiles(ourSmi)
    if ourMol is None:
        print(f"{name:12s} OUR SMILES IS INVALID: {ourSmi}")
        mismatches.append((name, "INVALID", info["canonical_smiles"]))
        continue
    
    ourInchi = Chem.MolToInchi(ourMol)
    ourCanon = Chem.MolToSmiles(ourMol, isomericSmiles=True)
    
    # Compare InChI (ignoring anomeric C1 differences)
    match = ourInchi == info["rdkit_inchi"]
    
    if not match:
        print(f"{name:12s} MISMATCH!")
        print(f"  Ours:    {ourCanon}")
        print(f"  PubChem: {info['canonical_smiles']}")
        print(f"  Ours InChI:    {ourInchi}")
        print(f"  PubChem InChI: {info['rdkit_inchi']}")
        mismatches.append((name, ourCanon, info["canonical_smiles"]))
    else:
        print(f"{name:12s} OK")

print(f"\n{'='*100}")
print(f"Mismatches: {len(mismatches)} / {len(results)}")

# Output Python code for the correct SMILES
print("\n\n# ============================================================")
print("# PubChem Authoritative SMILES (canonical, for patching)")
print("# ============================================================")
for name, info in sorted(results.items()):
    print(f'    ("{name}", "a"): "{info["canonical_smiles"]}",  # PubChem')
