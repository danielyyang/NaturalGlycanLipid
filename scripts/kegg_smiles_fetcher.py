"""
KEGG Monosaccharide SMILES Fetcher
Batch-fetch MOL files from KEGG REST API, convert to canonical SMILES with RDKit.
Output: authoritative SMILES for rebuilding monosaccharide library.
"""
import sys, os, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from rdkit import Chem
import urllib.request

# All key sugars from KEGG br08082 with their compound IDs
KEGG_SUGARS = {
    # === Hexoses ===
    "D-Glc": "C00031",
    "D-Galf": "C21066",  # Galactofuranose
    "D-Gal": "C00124",
    "L-Gal": "C01825",
    "D-Man": "C00159",
    "D-All": "C01487",
    "L-Alt": "C21032",
    "D-Gul": "C06465",
    "L-Ido": "C21050",
    "D-Tal": "C06467",

    # === Pentoses ===
    "D-Rib-p": "C00121",  # pyranose
    "D-Rib-f": "C21057",  # furanose
    "D-Ara-f": "C21067",  # furanose
    "D-Ara-p": "C00216",  # pyranose
    "L-Ara-f": "C06115",  # furanose
    "L-Ara-p": "C00259",  # pyranose
    "D-Xyl": "C00181",
    "D-Lyx": "C00476",

    # === Deoxy Sugars ===
    "6dL-Alt": "C21033",
    "6dD-Gul": "C22093",
    "6dD-Tal": "C21058",
    "D-Fuc": "C01018",
    "L-Fuc": "C01019",
    "D-Rha": "C01684",
    "L-Rha": "C00507",
    "D-Qui": "C02522",
    "2dD-Glc": "C00586",

    # === 3,6-Dideoxy Sugars ===
    "Oli": "C21054",   # Olivose
    "Tyv": "C21062",   # Tyvelose
    "Asc": "C22244",   # Ascarylose
    "Abe": "C06471",   # Abequose
    "Par": "C21055",   # Paratose
    "Dig": "C21045",   # Digitoxose
    "Col": "C03348",   # Colitose

    # === Amino Sugars ===
    "D-GlcN": "C00329",
    "D-GalN": "C02262",
    "D-ManN": "C03570",
    "D-AllN": "C21038",
    "L-AltN": "C21035",
    "D-GulN": "C21048",
    "L-IdoN": "C21051",
    "D-TalN": "C21060",

    # === N-Acetyl Amino Sugars ===
    "D-GlcNAc": "C00140",
    "D-GalNAc": "C01132",
    "D-ManNAc": "C00645",

    # === Uronic Acids ===
    "D-GlcA": "C00191",
    "D-GalA": "C00333",
    "D-ManA": "C02024",
    "D-GulA": "C21047",
    "L-GulA": "C06477",
    "L-IdoA": "C06472",

    # === Special Sugars ===
    "D-Api": "C21040",    # Apiose
    "The": "C16287",      # Thevetose
    "Cym": "C08234",      # Cymarose

    # === Sialic Acids ===
    "Neu5Ac": "C00270",
    "Neu5Gc": "C03410",

    # === Ketoses ===
    "D-Fru": "C00095",

    # === Special ===
    "KDO": "C21063",      # Ketodeoxyoctonic acid
}

results = {}
errors = []

for name, compId in KEGG_SUGARS.items():
    url = f"https://rest.kegg.jp/get/{compId}/mol"
    try:
        with urllib.request.urlopen(url, timeout=10) as resp:
            molBlock = resp.read().decode('utf-8')

        mol = Chem.MolFromMolBlock(molBlock)
        if mol is None:
            errors.append(f"{name} ({compId}): MolFromMolBlock returned None")
            continue

        smi = Chem.MolToSmiles(mol, isomericSmiles=True)
        results[name] = {
            "kegg_id": compId,
            "smiles": smi,
            "numAtoms": mol.GetNumAtoms(),
        }
        print(f"  {name:12s} ({compId}): {smi}")

    except Exception as e:
        errors.append(f"{name} ({compId}): {e}")

    time.sleep(0.3)  # polite rate limiting

print(f"\nFetched {len(results)}/{len(KEGG_SUGARS)} sugars")
if errors:
    print(f"\nErrors ({len(errors)}):")
    for err in errors:
        print(f"  {err}")

# Now output as Python dict for direct copy-paste
print("\n\n# ============================================================")
print("# KEGG Authoritative SMILES (for glycan_reference_library.py)")
print("# ============================================================")
print("RAW_MONOSACCHARIDE_SMILES = {")
for name, info in sorted(results.items()):
    print(f'    "{name}": "{info["smiles"]}",  # KEGG {info["kegg_id"]}')
print("}")
