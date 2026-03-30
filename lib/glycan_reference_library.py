from rdkit import Chem
from rdkit.Chem import rdqueries

# ==============================================================================
# 1. [DEPRECATED] 旧 SMARTS 库已删除 (Legacy SMARTS Library Removed)
# ==============================================================================
# 历史: 此处曾定义 _MONOSACCHARIDE_SMARTS_DEFS 和 MONOSACCHARIDE_LIBRARY。
# 审计 (2026-03-18): 0/7 SMARTS 无法匹配 PubChem 权威分子 (SubstructMatch 全部 False)。
# 核心匹配引擎 (identify_monosaccharide_v2) 完全依赖 Section 5 的 REFERENCE_MOLS。
# 为向后兼容保留空字典别名。
#
# History: _MONOSACCHARIDE_SMARTS_DEFS and MONOSACCHARIDE_LIBRARY were defined here.
# Audit (2026-03-18): 0/7 SMARTS failed to match PubChem authoritative molecules.
# Core engine uses REFERENCE_MOLS from Section 5 exclusively.
# Empty dict alias preserved for backward compatibility.
# ==============================================================================
MONOSACCHARIDE_LIBRARY = {}  # DEPRECATED — 使用 REFERENCE_MOLS 替代

# ==============================================================================
# 2. Amino Acid Library (Strict Amide Bond Enforced)
# ==============================================================================
# We define SMARTS patterns for Amino Acids that STRICTLY require the amino/carboxylic 
# ends to be attached via an AMIDE bond pattern: [N;!R]-[C;!R](=O)
# If a free amine/carboxyl is detected at the termini of a standalone molecule, we allow it,
# but internally, it must peptide-bond.

# Legend:
# [NX3v3,NX4v4+] = Nitrogen (amine/amide)
# [CX4H] = Alpha Carbon (sp3, exactly 1 H for standard L-AA except Gly)
# [CX3](=O)[OX2,NX3] = Carboxyl or Amide Carbonyl

_AMINO_ACID_SMARTS_DEFS = {
    'Ala': '[NX3,NX4+][CX4H](C)[CX3](=O)[OX2,NX3]',
    'Arg': '[NX3,NX4+][CX4H](CCCNC(=[NH2,NH])[NH2,NH3+])[CX3](=O)[OX2,NX3]',
    'Asn': '[NX3,NX4+][CX4H](CC(=O)N)[CX3](=O)[OX2,NX3]',
    'Asp': '[NX3,NX4+][CX4H](CC(=O)[OH,O-])[CX3](=O)[OX2,NX3]',
    'Cys': '[NX3,NX4+][CX4H](CS)[CX3](=O)[OX2,NX3]',
    'Gln': '[NX3,NX4+][CX4H](CCC(=O)N)[CX3](=O)[OX2,NX3]',
    'Glu': '[NX3,NX4+][CX4H](CCC(=O)[OH,O-])[CX3](=O)[OX2,NX3]',
    'Gly': '[NX3,NX4+][CX4H2][CX3](=O)[OX2,NX3]',
    'His': '[NX3,NX4+][CX4H](Cc1cncn1)[CX3](=O)[OX2,NX3]',
    'Ile': '[NX3,NX4+][CX4H](C(C)CC)[CX3](=O)[OX2,NX3]',
    'Leu': '[NX3,NX4+][CX4H](CC(C)C)[CX3](=O)[OX2,NX3]',
    'Lys': '[NX3,NX4+][CX4H](CCCCN)[CX3](=O)[OX2,NX3]',
    'Met': '[NX3,NX4+][CX4H](CCSC)[CX3](=O)[OX2,NX3]',
    'Phe': '[NX3,NX4+][CX4H](Cc1ccccc1)[CX3](=O)[OX2,NX3]',
    'Pro': '[NX3,NX4+;r5]1[CX4H;r5](CCC1)[CX3](=O)[OX2,NX3]',
    'Ser': '[NX3,NX4+][CX4H](CO)[CX3](=O)[OX2,NX3]',
    'Thr': '[NX3,NX4+][CX4H](C(O)C)[CX3](=O)[OX2,NX3]',
    'Trp': '[NX3,NX4+][CX4H](Cc1c[nH]c2ccccc12)[CX3](=O)[OX2,NX3]',
    'Tyr': '[NX3,NX4+][CX4H](Cc1ccc(O)cc1)[CX3](=O)[OX2,NX3]',
    'Val': '[NX3,NX4+][CX4H](C(C)C)[CX3](=O)[OX2,NX3]',
}

AMINO_ACID_LIBRARY = {k: Chem.MolFromSmarts(v) for k, v in _AMINO_ACID_SMARTS_DEFS.items()}

# A strict sequence definition enforcing standard peptide backbone matching
STRICT_PEPTIDE_BOND_SMARTS = Chem.MolFromSmarts("[NX3][CX3](=O)")

# ==============================================================================
# 3. Nucleobase Library
# ==============================================================================
NUCLEOBASE_LIBRARY = {
    "Adenine": Chem.MolFromSmarts("n1cnc2c1ncnc2"),
    "Guanine": Chem.MolFromSmarts("n1cnc2c1nc(N)[nH]c2=O"),
    "Cytosine": Chem.MolFromSmarts("Nc1ccn([#6])c(=O)n1"), # Attached to sugar C
    "Thymine": Chem.MolFromSmarts("Cc1cn([#6])c(=O)[nH]c1=O"),
    "Uracil": Chem.MolFromSmarts("O=c1ccn([#6])c(=O)[nH]1")
}

# Phosphate rule for Nucleotides
PHOSPHATE_SMARTS = Chem.MolFromSmarts("P(=O)(O)(O)O")

# Reaction SMARTS (SMIRKS) for Modification Stripping
_STRIPPING_SMIRKS_DEFS = {
    # === 1. 无机酸修饰 (Inorganic Esters) ===
    "Sulfated": "[O:1]S(=O)(=O)[O-,OH]",        # 硫酸酯化 (O-Sulfate)
    "Phosphated": "[O:1]P(=O)([O-,OH])[O-,OH]", # 磷酸酯化 (O-Phosphate)
    # === 2. 基础与复杂脂肪酸酯 (Aliphatic Esters) ===
    "Acetylated": "[O:1]C(=O)[CH3]",            # 乙酰化 (O-Acetyl)
    "Formylated": "[O:1]C(=O)[H]",              # 甲酰化 (O-Formyl)
    "Malonylated": "[O:1]C(=O)CC(=O)[OH,O-]",   # 丙二酰化 (O-Malonyl)
    "Succinylated": "[O:1]C(=O)CCC(=O)[OH,O-]", # 琥珀酰化 (O-Succinyl)
    "Lactylated": "[O:1]C(C)C(=O)[OH,O-]",      # 乳酰化 (O-Lactyl)
    "Tigloylated": "[O:1]C(=O)C(=C)C",          # 巴豆酰/当归酰化 (O-Tigloyl/Angeloyl)
    # === 3. 芳香族高亮修饰 (Aromatic Esters - 天然产物重灾区) ===
    "Galloylated": "[O:1]C(=O)c1cc(O)c(O)c(O)c1",       # 没食子酰化 (O-Galloyl)
    "Benzoylated": "[O:1]C(=O)c1ccccc1",                # 苯甲酰化 (O-Benzoyl)
    "p-Coumaroylated": "[O:1]C(=O)/C=C/c1ccc(O)cc1",    # 对香豆酰化 (O-p-Coumaroyl)
    "Caffeoylated": "[O:1]C(=O)/C=C/c1ccc(O)c(O)c1",    # 咖啡酰化 (O-Caffeoyl)
    "Feruloylated": "[O:1]C(=O)/C=C/c1ccc(O)c(OC)c1",   # 阿魏酰化 (O-Feruloyl)
    "Sinapoylated": "[O:1]C(=O)/C=C/c1cc(OC)c(O)c(OC)c1", # 芥子酰化 (O-Sinapoyl)
    # === 4. 醚类修饰 (Ethers) ===
    "Methylated": "[O:1][CH3]",                 # 甲醚化 (O-Methyl, 已修复为安全末端甲基)
    # === 5. 氨基特种修饰 (N-Modifications) ===
    "N-Acetylated": "[N:1]C(=O)[CH3]",          # N-乙酰化 (N-Acetyl)
    "N-Glycolylated": "[N:1]C(=O)CO",           # N-羟乙酰化 (N-Glycolyl)
    "N-Sulfated": "[N:1]S(=O)(=O)[O-,OH]",      # N-硫酸酯化 (N-Sulfate)
    "N-Methylated": "[N:1][CH3]",               # N-甲基化 (N-Methyl)
    "N-Formylated": "[N:1]C(=O)[H]"             # N-甲酰化 (N-Formyl)
}

STRIPPING_REACTIONS = {}
for name, smirks in _STRIPPING_SMIRKS_DEFS.items():
    try:
        rxn = Chem.MolFromSmarts(smirks)
        if rxn:
            STRIPPING_REACTIONS[name] = rxn
    except Exception as e:
        print(f"Failed to load SMARTS for {name}: {e}")


# ==============================================================================
# 5. Comprehensive Monosaccharide SMILES Library (全量单糖 SMILES 库)
# ==============================================================================
# 从 monosaccharide_identifier.py 迁入 — 覆盖 8 大己糖全构型、脱氧糖、二脱氧/三脱氧糖、
# 氨基糖、糖醛酸、呋喃糖、酮糖、支链糖及细菌/海洋特有糖。
# Migrated from monosaccharide_identifier.py — Comprehensive monosaccharide library covering
# all 8 aldohexose stereoisomers, deoxy/dideoxy/trideoxy sugars, amino sugars,
# uronic acids, furanoses, ketoses, branched sugars, and bacterial/marine sugars.
# All SMILES are in closed-ring pyranose/furanose form with full stereochemistry.
# Key format: (sugar_name, anomer_label) -> SMILES
RAW_MONOSACCHARIDE_SMILES = {
    # =====================================================================

    # 1. Hexoses — KEGG Authority (2026-03-17)

    # C1 anomeric = C(O) unspecified -> matches both alpha/beta

    # All C2-C5 chirality verified against KEGG InChI

    # =====================================================================

    # D-Glucose — KEGG C00031

    ("D-Glc", "a"): "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",

    # D-Galactose — KEGG C00124 [C4 epimer of Glc]

    ("D-Gal", "a"): "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O",

    # L-Galactose — KEGG C01825

    ("L-Gal", "a"): "OC[C@@H]1OC(O)[C@@H](O)[C@H](O)[C@@H]1O",

    # D-Mannose — KEGG C00159 [C2 epimer of Glc]

    ("D-Man", "a"): "OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O",

    # D-Allose — KEGG C01487 [C3 epimer of Glc]

    ("D-All", "a"): "OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@@H]1O",

    # D-Talose — KEGG C06467 [C2 epimer of Gal]

    ("D-Tal", "a"): "OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@H]1O",

    # D-Altrose — derived from L-Alt KEGG C21032

    ("D-Alt", "a"): "OC[C@H]1OC(O)[C@@H](O)[C@H](O)[C@@H]1O",

    # D-Gulose — KEGG C06465

    ("D-Gul", "a"): "OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@H]1O",

    # D-Idose — derived from L-Ido KEGG C21050

    ("D-Ido", "a"): "OC[C@H]1OC(O)[C@@H](O)[C@H](O)[C@H]1O",

    # L-Glucose — enantiomer of D-Glc

    ("L-Glc", "a"): "OC[C@@H]1OC(O)[C@@H](O)[C@H](O)[C@H]1O",

    # L-Mannose — enantiomer of D-Man

    ("L-Man", "a"): "OC[C@@H]1OC(O)[C@H](O)[C@H](O)[C@H]1O",

    # L-Altrose — KEGG C21032

    ("L-Alt", "a"): "OC[C@@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O",

    # L-Idose — KEGG C21050

    ("L-Ido", "a"): "OC[C@@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",

    # L-Talose — enantiomer of D-Tal

    ("L-Tal", "a"): "OC[C@@H]1OC(O)[C@H](O)[C@H](O)[C@@H]1O",

    # L-Allose — enantiomer of D-All

    ("L-All", "a"): "OC[C@@H]1OC(O)[C@@H](O)[C@@H](O)[C@H]1O",

    # L-Gulose — enantiomer of D-Gul

    ("L-Gul", "a"): "OC[C@@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O",



    # =====================================================================

    # 2. Deoxy Sugars — KEGG Authority

    # =====================================================================

    # L-Rhamnose — KEGG C00507 (6-deoxy-L-mannose)

    ("L-Rha", "a"): "C[C@@H]1OC(O)[C@H](O)[C@H](O)[C@H]1O",

    # D-Rhamnose — KEGG C01684 (6-deoxy-D-mannose)

    ("D-Rha", "a"): "C[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",

    # L-Fucose — KEGG C01019 (6-deoxy-L-galactose)

    ("L-Fuc", "a"): "C[C@@H]1OC(O)[C@@H](O)[C@H](O)[C@@H]1O",

    # D-Fucose — KEGG C01018 (6-deoxy-D-galactose)

    ("D-Fuc", "a"): "C[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O",

    # D-Quinovose — KEGG C02522 (6-deoxy-D-glucose)

    ("D-Qui", "a"): "C[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",

    # L-Quinovose — enantiomer of D-Qui

    ("L-Qui", "a"): "C[C@@H]1OC(O)[C@@H](O)[C@H](O)[C@H]1O",


    # =====================================================================
    # 3. 二脱氧糖 (Dideoxy Sugars)
    # =====================================================================
    ("D-Oli", "a"): "C[C@@H]1C[C@@H](O)[C@H](O)[C@@H](O)O1",
    ("D-Dig", "a"): "C[C@@H]1C[C@H](O)[C@H](O)[C@@H](O)O1",
    # 修正 (2026-03-27): OMe/OH 位置互换修复, 从 PubChem Oleandrin 碎片提取
    # 修正 (2026-03-27): 手工推导确认真正的 L-Oleandrose 物理构型应表现为 C3,C4,C5=(S,R,S)
    # Fix: Verified L-Oleandrose (2,6-dideoxy-3-O-methyl-L-arabino-hexose) physical chirality is S, R, S
    ("L-Ole", "a"): "C[C@@H]1OC(O)C[C@H](OC)[C@H]1O",
    ("D-Oliose", "a"): "C[C@@H]1C[C@@H](O)[C@@H](O)[C@@H](O)O1",
    ("D-Boi", "a"): "C[C@@H]1C[C@H](O)[C@@H](O)[C@@H](O)O1",
    ("D-Cym", "a"): "C[C@@H]1C[C@H](OC)[C@H](O)[C@@H](O)O1",

    # =====================================================================
    # 4. 三脱氧糖 (Trideoxy Sugars) — 细菌 LPS
    # =====================================================================
    ("D-Abe", "a"): "C[C@@H]1[C@@H](O)C[C@H](O)[C@@H](O)O1",
    ("D-Par", "a"): "C[C@@H]1[C@H](O)C[C@H](O)[C@@H](O)O1",
    ("D-Tyv", "a"): "C[C@@H]1[C@H](O)C[C@@H](O)[C@@H](O)O1",
    ("L-Col", "a"): "C[C@H]1C[C@@H](O)[C@@H](O)[C@@H](O)O1", # Restored 2026-03-27: Safe to include with v10 CIP engine
    ("L-Asc", "a"): "C[C@H]1[C@@H](O)C[C@@H](O)[C@H](O)O1",

    # =====================================================================
    # 5. 氨基糖 (Amino Sugars)
    # =====================================================================
    ("D-GlcNAc", "a"): "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O",  # KEGG C00140

    ("D-GlcNAc", "b"): "O[C@@H]1[C@H](NC(=O)C)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    ("D-GalNAc", "a"): "CC(=O)N[C@H]1C(O)O[C@H](CO)[C@H](O)[C@@H]1O",  # KEGG C01132

    ("D-GalNAc", "b"): "O[C@@H]1[C@H](NC(=O)C)[C@@H](O)[C@@H](O)[C@@H](CO)O1",
    ("D-ManNAc", "a"): "CC(=O)N[C@@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O",  # KEGG C00645

    ("D-ManNAc", "b"): "O[C@@H]1[C@@H](NC(=O)C)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    ("D-GlcN", "a"): "N[C@H]1C(O)O[C@H](CO)[C@@H](O)[C@@H]1O",  # KEGG C00329

    ("D-GalN", "a"): "N[C@H]1C(O)O[C@H](CO)[C@H](O)[C@@H]1O",  # KEGG C02262
    
    # Kanamycin A specific amino sugars:
    ("D-Kanosamine", "a"): "N[C@H]1[C@H](O)[C@@H](CO)OC(O)[C@@H]1O", # 3-amino-3-deoxy-D-glucose
    ("D-6aGlc", "a"): "NC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", # 6-amino-6-deoxy-D-glucose

    ("L-FucNAc", "a"): "C[C@@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@@H]1NC(=O)C",
    ("D-QuiNAc", "a"): "C[C@@H]1[C@H](NC(=O)C)[C@@H](O)[C@H](O)[C@@H](O)O1",
    ("L-RhaNAc", "a"): "C[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1NC(=O)C",
    ("D-Bac", "a"): "C[C@@H]1[C@H](N)C[C@H](N)[C@@H](O)O1",

    # =====================================================================
    # 6. 糖醛酸 (Uronic Acids)
    # =====================================================================
    ("D-GlcA", "a"): "O=C(O)[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",  # KEGG C00191

    ("D-GlcA", "b"): "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](C(=O)O)O1",
    ("D-GalA", "a"): "O=C(O)[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O",  # KEGG C00333

    ("D-GalA", "b"): "O[C@@H]1[C@H](O)[C@@H](O)[C@@H](O)[C@@H](C(=O)O)O1",
    ("D-ManA", "a"): "O[C@H]1[C@@H](O)[C@@H](O)[C@H](O)[C@@H](C(=O)O)O1",
    ("L-IdoA", "a"): "O=C(O)[C@@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",  # KEGG C06472

    ("L-IdoA", "b"): "O[C@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](C(=O)O)O1",
    ("L-GulA", "a"): "O[C@@H]1[C@@H](O)[C@@H](O)[C@H](O)[C@H](C(=O)O)O1",

    # =====================================================================
    # 7. 戊糖 (Pentoses) — 吡喃糖和呋喃糖
    # =====================================================================
    ("D-Xyl", "a"): "OC1OC[C@@H](O)[C@H](O)[C@H]1O",  # KEGG C00181

    ("D-Xyl", "b"): "O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)CO1",
    ("L-Ara", "a"): "OC1OC[C@H](O)[C@H](O)[C@H]1O",  # KEGG C00259

    ("L-Ara", "b"): "O[C@@H]1[C@@H](O)[C@@H](O)[C@@H](O)CO1",
    ("D-Ara", "a"): "O[C@@H]1[C@H](O)[C@H](O)[C@H](O)CO1",
    ("D-Ara", "b"): "O[C@H]1[C@H](O)[C@H](O)[C@H](O)CO1",
    ("D-Rib", "a"): "OC[C@H]1OC(O)[C@H](O)[C@@H]1O",  # KEGG C00121 (furanose)

    ("D-Rib", "b"): "O[C@H]1[C@H](O)[C@H](O)[C@@H](O)CO1",
    ("D-Lyx", "a"): "OC1OC[C@@H](O)[C@H](O)[C@@H]1O",  # KEGG C00476

    ("D-Lyx", "b"): "O[C@H]1[C@H](O)[C@@H](O)[C@H](O)CO1",
    ("L-Lyx", "a"): "O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)CO1",
    # --- 呋喃糖 (Furanoses) ---
    ("D-Ribf", "a"): "OC[C@H]1O[C@H](O)[C@H](O)[C@@H]1O",
    ("D-Ribf", "b"): "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H]1O",
    ("L-Araf", "a"): "OC[C@@H]1O[C@@H](O)[C@H](O)[C@@H]1O",
    ("L-Araf", "b"): "OC[C@@H]1O[C@H](O)[C@H](O)[C@@H]1O",
    ("D-Araf", "a"): "OC[C@H]1O[C@H](O)[C@@H](O)[C@H]1O",
    ("D-Xylf", "a"): "OC[C@@H]1O[C@H](O)[C@H](O)[C@@H]1O",
    ("D-Galf", "a"): "OC[C@H](O)[C@H]1O[C@H](O)[C@H](O)[C@@H]1O",
    ("D-Galf", "b"): "OC[C@H](O)[C@H]1O[C@@H](O)[C@H](O)[C@@H]1O",
    ("D-Glcf", "a"): "OC[C@H](O)[C@@H]1O[C@H](O)[C@H](O)[C@@H]1O",

    # =====================================================================
    # 8. 酮糖 (Ketoses)
    # =====================================================================
    ("D-Fru", "a"): "OC[C@]1(O)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    ("D-Fru", "b"): "OC[C@@]1(O)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    ("D-Frup", "a"): "OC[C@]1(O)OC[C@@H](O)[C@@H](O)[C@@H]1O",
    ("D-Fruf", "a"): "OC[C@H]1O[C@@](O)(CO)[C@@H](O)[C@@H]1O",
    ("L-Sorb", "a"): "OC[C@@]1(O)OC[C@H](O)[C@@H](O)[C@@H]1O",
    ("D-Tag", "a"): "OC[C@]1(O)OC[C@H](O)[C@H](O)[C@@H]1O",
    ("D-Psi", "a"): "OC[C@]1(O)OC[C@H](O)[C@@H](O)[C@H]1O",

    # =====================================================================
    # 9. 支链糖 (Branched-chain Sugars)
    # =====================================================================
    ("D-Api", "a"): "OC[C@@]1(O)COC(O)[C@@H]1O",  # KEGG C21040

    ("D-Api", "b"): "O[C@@H]1[C@@H](O)[C@@](O)(CO)CO1",

    # =====================================================================
    # 10. 庚糖及特殊糖 (Heptoses & Special Sugars)
    # =====================================================================
    ("D-Sed", "a"): "OC[C@H]1O[C@@](O)(CO)[C@@H](O)[C@H](O)[C@@H]1O",
    ("L-D-Hep", "a"): "OC[C@H](O)[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
    ("D-D-Hep", "a"): "OC[C@@H](O)[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
    # L-glycero-D-galacto-heptose — LPS 核心寡糖 (LPS core oligosaccharide)
    ("L-D-galHep", "a"): "OC[C@H](O)[C@H]1O[C@H](O)[C@H](O)[C@H](O)[C@@H]1O",
    # D-glycero-D-galacto-heptose — 植物细胞壁 (Plant cell wall)
    ("D-D-galHep", "a"): "OC[C@@H](O)[C@H]1O[C@H](O)[C@H](O)[C@H](O)[C@@H]1O",
    # D-glycero-D-altro-heptose — 细菌特殊庚糖 (Bacterial special heptose)
    ("D-D-altHep", "a"): "OC[C@@H](O)[C@@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",
    # D-glycero-D-gluco-heptose — 革兰氏阴性菌 LPS (Gram-negative LPS)
    ("D-D-glcHep", "a"): "OC[C@@H](O)[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O",

    # =====================================================================
    # 11. 高级多碳糖 (Advanced Multi-Carbon Sugars)
    # =====================================================================
    # Wiki (N-Acetylneuraminic acid / Sialic acid) — C2 酮糖型季碳
    ("Neu5Ac", "a"): "OC(=O)[C@@]1(O)C[C@H](O)[C@@H](NC(C)=O)[C@@H](O1)[C@H](O)[C@H](O)CO",
    ("Neu5Ac", "b"): "OC(=O)[C@]1(O)C[C@H](O)[C@@H](NC(C)=O)[C@@H](O1)[C@H](O)[C@H](O)CO",
    ("KDO", "a"): "O=C(O)[C@@]1(O)C[C@H](O)[C@@H](O)[C@@H](O1)[C@@H](O)CO",
    ("KDO", "b"): "O=C(O)[C@]1(O)C[C@H](O)[C@@H](O)[C@@H](O1)[C@@H](O)CO",
    ("Neu5Gc", "a"): "O=C(O)[C@@]1(O)C[C@H](O)[C@@H](NC(=O)CO)[C@@H](O1)[C@@H](O)[C@@H](O)CO",

    # =====================================================================
    # 12. Extended Amino Sugars
    # =====================================================================
    ("D-ManN", "a"): "O[C@H]1[C@@H](N)[C@@H](O)[C@H](O)[C@@H](CO)O1",
    ("D-AllN", "a"): "OC[C@H]1O[C@H](O)[C@H](N)[C@H](O)[C@@H]1O",
    ("D-Des", "a"): "C[C@@H]1C[C@@H](O)[C@H](N(C)C)[C@@H](O)O1",
    ("D-MurNAc", "a"): "O[C@H]1[C@H](NC(=O)C)[C@@H](O[C@@H](C)C(=O)O)[C@H](O)[C@@H](CO)O1",
    ("D-GulNAc", "a"): "O[C@H]1[C@H](NC(=O)C)[C@H](O)[C@@H](O)[C@@H](C(=O)O)O1",
    ("D-TalNAc", "a"): "OC[C@H]1O[C@H](O)[C@@H](NC(=O)C)[C@@H](O)[C@@H]1O",
    ("D-AltNAc", "a"): "OC[C@H]1O[C@H](O)[C@@H](NC(=O)C)[C@H](O)[C@@H]1O",
    ("D-AllNAc", "a"): "OC[C@H]1O[C@H](O)[C@H](NC(=O)C)[C@H](O)[C@@H]1O",
    ("D-IdoNAc", "a"): "OC[C@H]1O[C@H](O)[C@@H](NC(=O)C)[C@@H](O)[C@H]1O",
    ("L-IdoNAc", "a"): "OC[C@@H]1O[C@@H](O)[C@H](NC(=O)C)[C@H](O)[C@@H]1O",
    ("D-MurNGc", "a"): "O[C@H]1[C@H](NC(=O)CO)[C@@H](O[C@@H](C)C(=O)O)[C@H](O)[C@@H](CO)O1",
    ("L-FucN", "a"): "C[C@@H]1O[C@H](O)[C@@H](O)[C@H](N)[C@@H]1O",

    # =====================================================================
    # 13. Extended Uronic Acids
    # =====================================================================
    ("D-AllA", "a"): "O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@@H](C(=O)O)O1",
    ("D-TalA", "a"): "OC(=O)[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",
    ("D-AltA", "a"): "O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@@H](C(=O)O)O1",
    ("L-GalA", "a"): "O[C@@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H](C(=O)O)O1",
    ("D-meGlcA", "a"): "O[C@H]1[C@H](O)[C@@H](OC)[C@H](O)[C@@H](C(=O)O)O1",

    # =====================================================================
    # 14. Extended Deoxy Sugars
    # =====================================================================
    ("D-dRib", "a"): "OC[C@H]1C[C@H](O)[C@@H](O)O1",
    ("D-dGlc", "a"): "OC[C@H]1OC(O)[C@H](O)[C@@H](O)C1",
    ("D-dAra", "a"): "O[C@@H]1C[C@H](O)[C@H](O)CO1",
    ("D-dXyl", "a"): "O[C@H]1C[C@@H](O)[C@H](O)CO1",
    ("L-Cla", "a"): "C[C@H]1O[C@H](O)C[C@@](C)(OC)[C@@H]1O",
    ("L-Aco", "a"): "C[C@H]1C[C@@H](N)[C@H](O)[C@@H](O)O1",
    ("D-The", "a"): "C[C@@H]1[C@H](O)[C@@H](OC)[C@H](O)[C@@H](O)O1",
    ("D-6dTalNAc", "a"): "C[C@H]1O[C@H](O)[C@@H](NC(=O)C)[C@@H](O)[C@@H]1O",
    ("D-6dGul", "a"): "C[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    ("D-6dTal", "a"): "C[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",
    ("D-6dAlt", "a"): "C[C@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O",

    # =====================================================================
    # 14b. 心苷糖 (Cardiac Glycoside Sugars) — 2,6-二脱氧/2,3,6-三脱氧糖
    # Cardiac glycoside sugars: 2,6-dideoxy and 2,3,6-trideoxy hexoses
    # 关键: 洋地黄毒苷(Digitoxin)/地高辛(Digoxin)/毛花苷C(Lanatoside C) 等
    # Critical for: Digitoxin, Digoxin, Lanatoside C, Strophanthidin glycosides
    # SMILES 已通过碎片提取+暴力手性匹配确认 (Confirmed by fragment extraction)
    # =====================================================================
    # D-Digitoxose — 2,6-dideoxy-D-ribo-hexose (KEGG C02657)
    # 洋地黄类心苷的核心糖: Digoxin 含 3 个 D-Digitoxose
    # Core sugar of digitalis cardiac glycosides
    ("D-Dtx", "a"): "C[C@H]1OC(O)C[C@H](O)[C@@H]1O",

    # L-Digitoxose — enantiomer
    ("L-Dtx", "a"): "C[C@@H]1OC(O)C[C@@H](O)[C@H]1O",

    # D-Sarmentose — 2,6-dideoxy-3-O-methyl-D-ribo-hexose
    # 夹竹桃苷等的组分 (Component of sarmentogenin glycosides)
    ("D-Sar", "a"): "C[C@H]1OC(O)C[C@H](OC)[C@@H]1O",

    # D-Cymarose — 2,6-dideoxy-3-O-methyl-D-ribo-hexose (different C-3 config)
    # 铃兰毒苷等 (Convallatoxin, etc.)
    ("D-Cma", "a"): "C[C@H]1OC(O)C[C@@H](OC)[C@@H]1O",

    # D-Diginose — 2,6-dideoxy-D-lyxo-hexose
    # C3 差向异构体 (C3 epimer of digitoxose)
    ("D-Din", "a"): "C[C@H]1OC(O)C[C@@H](O)[C@@H]1O",

    # L-Oleandrose — 2,6-dideoxy-3-O-methyl-L-arabino-hexose
    # 夹竹桃糖 (Oleander sugar)
    ("L-Ola", "a"): "C[C@@H]1OC(O)C[C@@H](OC)[C@H]1O",

    # D-Evalose — 2,3,6-trideoxy-D-hexose (no OH on C2 AND C3)
    # 链霉菌次生代谢物 (Streptomyces secondary metabolites)
    ("D-Eva", "a"): "C[C@H]1OC(O)CC[C@@H]1O",

    # L-Evalose — enantiomer
    ("L-Eva", "a"): "C[C@@H]1OC(O)CC[C@H]1O",

    # D-Desosamine — 3-(dimethylamino)-3,4,6-trideoxy-D-xylo-hexopyranose
    # 大环内酯抗生素核心糖: Erythromycin (红霉素)
    # Core sugar of macrolide antibiotics: Erythromycin
    # C3 = -N(CH3)2, C4 = -H (deoxy), C6 = -CH3 (deoxy)
    ("D-Des", "a"): "C[C@@H]1C[C@H](N(C)C)[C@@H](O)C(O)O1",

    # L-Cladinose — 2,6-dideoxy-3-C-methyl-3-O-methyl-L-ribo-hexopyranose
    # Erythromycin 的第二个糖: 含 C3 支链甲基 + C3-OMe
    # Second sugar of Erythromycin: C3 branched methyl + C3-OMe
    # 注意: C3 是季碳 (quaternary carbon), CIP 提取引擎需兼容
    # 修正 (2026-03-27): OMe/OH 位置互换修复, 从 PubChem Erythromycin 碎片提取
    # Fix: OMe/OH position swap corrected, extracted from PubChem Erythromycin fragment
    ("L-Cla", "a"): "CO[C@]1(C)C[C@H](O)O[C@@H](C)[C@@H]1O",

    # =====================================================================
    # 15. Tetroses
    # =====================================================================
    ("D-Thr", "a"): "O[C@H]1[C@H](O)[C@@H](O)CO1",
    ("D-Ery", "a"): "O[C@@H]1[C@H](O)[C@@H](O)CO1",

    # =====================================================================
    # 16. Pentose Ketoses
    # =====================================================================
    ("D-Ribu", "a"): "OC[C@@]1(O)[C@@H](O)[C@@H](O)CO1",
    ("D-Xylu", "a"): "OC[C@]1(O)[C@H](O)[C@@H](O)CO1",

    # =====================================================================
    # 17. Extended Nonoses (9-Carbon)
    # =====================================================================
    ("Kdn", "a"): "O=C(O)[C@@]1(O)C[C@H](O)[C@@H](O)[C@@H](O1)[C@@H](O)[C@@H](O)CO",
    ("Pse", "a"): "O=C(O)[C@@]1(O)C[C@@H](O)[C@@H](NC(=O)C)[C@@H](O1)[C@H](O)[C@H](NC(=O)C)C",
    ("Leg", "a"): "O=C(O)[C@@]1(O)C[C@H](O)[C@@H](NC(=O)C)[C@@H](O1)[C@@H](O)[C@@H](NC(=O)C)C",
    ("Fus", "a"): "O=C(O)[C@@]1(O)C[C@H](O)[C@@H](NC(=O)C)[C@@H](O1)[C@@H](O)[C@@H](O)CO",

    # =====================================================================
    # 18. Extended Octoses (8-Carbon)
    # =====================================================================
    ("Oct", "a"): "O=C(O)[C@]1(O)C[C@@H](O)[C@@H](O)[C@@H](O1)[C@@H](O)CO",

    # =====================================================================
    # 19. Decoses (10-Carbon)
    # =====================================================================
    ("Dha", "a"): "O=C(O)[C@@]1(O)C[C@H](O)[C@@H](O)[C@@H](O1)[C@@H](O)[C@@H](O)[C@@H](O)CO",

    # =====================================================================
    # 20. 抗生素特殊糖 (Antibiotic Special Sugars)
    # =====================================================================
    # L-Streptose — 5-deoxy-3-C-formyl-L-lyxose (分支醛戊糖)
    # Branched aldopentose in Streptomycin. C3 carries -CHO branch.
    # 吡喃型: 5元环 (furanose-like hemiacetal) → 实际为 5C 分支糖
    # PubChem CID 5460942, C6H10O5, MW=162.14
    ("L-Streptose", "a"): "C[C@@H]1O[C@@H](O)[C@H](O)[C@@]1(O)C=O",

    # L-Vancosamine — 3-amino-2,4,6-trideoxy-3-C-methyl-L-lyxo-hexose
    # Vancomycin 的脱氧氨基糖. 含 C3 季碳 (methyl + amino)
    # PubChem CID 189099, C7H15NO3, MW=161.20
    ("L-Vancosamine", "a"): "C[C@@H]1OC(O)C[C@](C)(N)[C@@H]1O",

    # D-Mycosamine — 3-amino-3,6-dideoxy-D-mannose
    # Amphotericin B 的氨基糖 (aminosugar of polyene antifungals)
    # PubChem CID 182095, C6H13NO4, MW=163.17
    ("D-Myc", "a"): "C[C@H]1OC(O)[C@@H](O)[C@@H](N)[C@@H]1O",

    # dRib — 2-deoxy-D-ribofuranose (nucleoside core sugar)
    # 核苷碱基连接的脱氧核糖 (furanose form, 5-membered ring)
    # PubChem CID 5460005, C5H10O4, MW=134.13
    ("dRib", "a"): "OC[C@H]1CC(O)O[C@@H]1O",

    # D-Kanosamine — 3-amino-3-deoxy-D-glucose
    # 卡那霉素 A (Kanamycin A) 的组分氨基糖
    ("D-Kanosamine", "a"): "OC[C@H]1OC(O)[C@H](O)[C@@H](N)[C@@H]1O",

    # D-6aGlc — 6-amino-6-deoxy-D-glucose
    # 卡那霉素 A (Kanamycin A) 的组分氨基糖
    ("D-6aGlc", "a"): "NC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
}


# ==============================================================================
# 铁律 1 (修正): 异头碳手性剥离 + C2-C5 严格保留
# Iron Rule 1 (Revised): Anomeric Relaxation + Strict C2-C5 Chirality
# ==============================================================================
# 设计意图 (Design Intent):
#   异头碳 (C1, 醛糖) 的 α/β 手性是糖苷键连接属性, 不是糖骨架身份标识。
#   在单糖识别阶段强制匹配 C1 手性会导致大量合法糖被误判为 Hex。
#   解决方案: 加载模板时程序化剥离异头碳手性, 保持 C2-C5 绝对严格。
#
#   The anomeric carbon (C1 for aldoses) α/β chirality is a glycosidic bond
#   linkage property, NOT a sugar identity feature. Forcing C1 chirality
#   matching causes massive false Hex results. Solution: programmatically
#   strip anomeric chirality at load time, keep C2-C5 absolutely strict.
# ==============================================================================


def _findAnomericCarbon(mol):
    """定位糖环中的异头碳 (Anomeric Carbon Finder)。
    Locate the anomeric carbon in a sugar ring.

    异头碳特征: 糖环内的碳原子, 同时与环内氧原子和环外氧原子相连。
    对于醛糖是 C1, 对于酮糖是 C2。
    Anomeric carbon: ring carbon bonded to both the ring oxygen
    AND an exocyclic oxygen (OH or OR). C1 for aldoses, C2 for ketoses.

    Returns:
        int or None: 异头碳的原子索引, 未找到则返回 None。
    """
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        # 寻找含 1 个氧的 5 元或 6 元环 (pyranose/furanose)
        # Find 5- or 6-membered ring with exactly 1 oxygen
        if len(ring) not in (5, 6):
            continue
        ringSet = set(ring)
        ringOxygens = [i for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() == 8]
        if len(ringOxygens) != 1:
            continue

        # 在环碳中寻找异头碳: 同时连接环内 O 和环外 O
        # Among ring carbons, find the one bonded to both ring-O and exo-O
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            hasRingO = False
            hasExoO = False
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    if nbr.GetIdx() in ringSet:
                        hasRingO = True
                    else:
                        hasExoO = True
            if hasRingO and hasExoO:
                return idx
    return None


REFERENCE_MOLS = {}
_anomericStripped = 0
for k, smiles in RAW_MONOSACCHARIDE_SMILES.items():
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        continue

    # 异头碳手性剥离 (Anomeric Relaxation):
    # 找到异头碳, 擦除其手性标记, 使 α/β 均可匹配
    # Strip anomeric carbon chirality so both α and β anomers match
    anomericIdx = _findAnomericCarbon(mol)
    if anomericIdx is not None:
        atom = mol.GetAtomWithIdx(anomericIdx)
        atom.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)
        _anomericStripped += 1

    # MolToSmarts: 保留 C2-C5 手性 (@/@@), 放松 O 原子的度约束
    # MolToSmarts: preserves C2-C5 chirality, relaxes O degree constraints
    smartsStr = Chem.MolToSmarts(mol)
    refMol = Chem.MolFromSmarts(smartsStr)
    if refMol is not None:
        REFERENCE_MOLS[k] = refMol


# ==============================================================================
# 铁律 3: 特异性排序遍历顺序 (Iron Rule 3: Specificity-Driven Priority)
# ==============================================================================
# 匹配引擎遍历字典时, 必须按"特异性从高到低"的顺序。
# 罕见/修饰度高的糖先匹配 → 常见基础糖最后匹配。
# 这防止了"先到先得"导致宽松模板劫持精确匹配。
# Matching engine must iterate dictionary in "most specific first" order.
# Rare/highly-modified sugars first → common base sugars last.

# 环大小元数据: 判断参考糖是呋喃(5元环)还是吡喃(6元环)
# Ring size metadata: determine if reference sugar is furanose(5) or pyranose(6)
REFERENCE_RING_SIZE = {}
for k, smiles in RAW_MONOSACCHARIDE_SMILES.items():
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        continue
    ri = mol.GetRingInfo()
    ringSize = 6  # 默认吡喃
    for ring in ri.AtomRings():
        oCount = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if oCount >= 1 and len(ring) in (5, 6):
            ringSize = len(ring)
            break
    REFERENCE_RING_SIZE[k] = ringSize

# 特异性排序: 原子数多 / 约束多的先匹配
# Specificity order: more atoms/constraints = matched FIRST
# 关键原理: 在子结构匹配中, 小分子永远是大分子的子图。
# 所以 D-dGlc (少一个OH) 永远是 D-Glc 的子图 → D-Glc 必须先匹配!
# Key principle: in substructure matching, smaller molecules are always
# subgraphs of larger ones. So D-dGlc is always a subgraph of D-Glc
# → D-Glc MUST be tried first!

# 预计算每个参考糖的重原子数 (用作同组内的二级排序)
REFERENCE_HEAVY_ATOM_COUNT = {}
for k, smiles in RAW_MONOSACCHARIDE_SMILES.items():
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        REFERENCE_HEAVY_ATOM_COUNT[k] = mol.GetNumHeavyAtoms()

def _getSpecificityPriority(key):
    """按子结构匹配的正确特异性返回优先级 (越小越先匹配)。
    Return priority based on correct substructure-matching specificity.
    核心原则: 原子多/约束多的分子先匹配, 因为它们是最特异的。
    脱氧变体 (原子少) 是基础糖的子图 → 必须最后匹配。
    Core: molecules with more atoms/constraints try first (most specific).
    Deoxy variants (fewer atoms) are subgraphs of base sugars → last.
    """
    name = key[0]
    hac = REFERENCE_HEAVY_ATOM_COUNT.get(key, 0)

    # 优先级 0: 大碳糖 (9C+, 最多原子, 最特异) — Neu5Ac, KDO, Kdn 等
    if name in ("Neu5Ac", "Neu5Gc", "KDO", "Kdn", "Pse", "Leg", "Fus", "Oct", "Dha"):
        return (0, -hac)
    # 优先级 1: 氨基糖 (GlcNAc/GalNAc 等带 NAc 基团, 原子多)
    if "NAc" in name or "Mur" in name:
        return (1, -hac)
    # 优先级 2: 糖醛酸 (GlcA/GalA 等带 COOH, 原子多)
    if name.endswith("A") and name not in ("L-Ara", "D-Ara"):
        return (2, -hac)
    # 优先级 3: 氨基糖 (GlcN/GalN 等不带 Ac 但带 N) + 特殊抗生素氨基糖
    if (name.endswith("N") and not name.endswith("NAc") or name == "D-Bac" or name == "D-Des"
            or name in ("L-Vancosamine", "D-Myc", "L-Aco", "D-Kanosamine", "D-6aGlc")):
        return (3, -hac)
    # 优先级 4: 酮糖 (结构有酮基特征)
    if name in ("D-Fru", "D-Frup", "D-Fruf", "L-Sorb", "D-Tag", "D-Psi", "D-Ribu", "D-Xylu"):
        return (4, -hac)
    # 优先级 5: 支链/庚糖 (Api/Sed/Hep 有独特拓扑) + 分支糖
    if name in ("D-Api", "D-Sed", "L-D-Hep", "D-D-Hep", "L-Streptose"):
        return (5, -hac)
    # 优先级 6: 基础己糖 (Glc, Gal, Man 等 — 原子多, 最特异的基础糖)
    if name in ("D-Glc", "L-Glc", "D-Gal", "L-Gal", "D-Man", "L-Man",
                "D-All", "D-Alt", "D-Gul", "D-Ido", "D-Tal"):
        return (6, -hac)
    # 优先级 7: 戊糖 (Xyl, Ara, Rib, Lyx — 原子适中, 在脱氧糖之前!) + 核苷糖
    if name in ("D-Xyl", "L-Ara", "D-Ara", "D-Rib", "D-Lyx", "L-Lyx", "D-Thr", "D-Ery", "dRib"):
        return (7, -hac)
    # 优先级 8: 呋喃糖 (f 后缀)
    if name.endswith("f"):
        return (8, -hac)
    # 优先级 9: 单脱氧糖 (L-Rha, L-Fuc — 6dHex, 原子比基础己糖少 1 个 O)
    if name in ("L-Rha", "D-Rha", "L-Fuc", "D-Fuc", "D-Qui", "L-Qui"):
        return (9, -hac)
    # 优先级 10: 二脱氧/三脱氧糖 (更少原子, 是基础糖的子图)
    if name in ("D-Oli", "D-Dig", "L-Ole", "D-Oliose", "D-Boi", "D-Cym",
                "D-Abe", "D-Par", "D-Tyv", "L-Asc"):
        return (10, -hac)
    # 优先级 11: 扩展脱氧 (D-dGlc, D-dXyl 等 — 原子最少, 最后匹配!)
    if name.startswith(("D-d", "L-d")) or name.startswith("6d") or name in ("D-Cla", "D-Aco", "D-The"):
        return (11, -hac)
    # 优先级 12: 四碳糖 (Thr, Ery — 最小结构)
    return (12, -hac)

# 构建排序后的遍历键列表
SPECIFICITY_ORDER = sorted(REFERENCE_MOLS.keys(), key=_getSpecificityPriority)


# ==============================================================================
# 6. Modification Scan SMARTS (修饰基团扫描 SMARTS)
# ==============================================================================
# 从 modification_scanner.py 迁入 — 用于扫描糖链片段上的常见化学修饰基团。
# Migrated from modification_scanner.py — scan glycan fragments for modifications.
from collections import OrderedDict

MODIFICATION_SMARTS: "OrderedDict[str, str]" = OrderedDict([
    ("NAc",       "[N;!$(N=*);!$(N#*)]-C(=O)[CH3]"),
    ("O-Ac",      "[O;!$(O=*)]-C(=O)[CH3]"),
    # 甲酰化 (O-Formyl) — 专家审查后补充 (Added per expert review)
    ("O-CHO",     "[OX2;!$(O=*)]-[CH1]=O"),                                    # 甲酰酯 (Formyl ester)
    ("Sulfate",   "[O;!$(O=*)]-S(=O)(=O)-[O,OH,$([O-])]"),
    ("Phosphate", "[O;!$(O=*)]-P(=O)(-[O,OH,$([O-])])-[O,OH,$([O-])]"),
    ("O-Me",      "[O;!$(O=*);!$(*C=O)]-[CH3;!$([CH3]C(=O))]"),
    ("COOH",      "[CX3](=O)[OX2H1,$([OX1-])]"),
    ("Acetonide", "[CX4]([CH3])([CH3])([OX2])([OX2])"),
    ("NH2",       "[NX3H2;!$(NC=O)]"),
    # 芳香族保护基/修饰 (Aromatic protecting groups / modifications)
    ("O-Bn",      "[OX2]-[CH2]-c1ccccc1"),                                    # 苄基醚 (Benzyl ether)
    # O-Bz: 注意 — 苯甲酰 SMARTS 是没食子酰的子结构, 需在 scanner 中后处理去重
    # O-Bz: NOTE — benzoyl is substructure of galloyl; dedup handled in modification_scanner
    ("O-Bz",      "[OX2]-C(=O)-c1ccccc1"),                                    # 苯甲酰酯 (Benzoyl ester)
    ("O-Gall",    "[OX2]-C(=O)-c1cc(O)c(O)c(O)c1"),                           # 没食子酰酯 (Galloyl ester)
    ("O-pCou",    "[OX2]-C(=O)-/C=C/-c1ccc(O)cc1"),                           # 对香豆酰 (p-Coumaroyl)
    ("O-Caf",     "[OX2]-C(=O)-/C=C/-c1cc(O)c(O)cc1"),                        # 咖啡酰 (Caffeoyl)
    ("O-Fer",     "[OX2]-C(=O)-/C=C/-c1cc(OC)c(O)cc1"),                       # 阿魏酰 (Feruloyl)
])

