"""
Test using RDKit's topological distance and relative orientation (stereo tags directly)
"""
import sys, os
from rdkit import Chem
from rdkit.Chem import AllChem

SUGAR_SMILES = {
    "Glc": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O",
    "Gal": "OC[C@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O",
    "Man": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    "All": "OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@@H]1O"
}

def get_chiral_tags(mol, ring_atoms):
    # Instead of full CIP R/S which flips, we can literally just ask RDKit for the ChiralTag 
    # of each atom in the ring. The chiral tag (CHI_TETRAHEDRAL_CW or CCW) is relative to 
    # the atom indices! 
    # For a canonical representation, we can orient by finding the O, then C1, then C2...
    # and determining if the exocyclic substituent is 'Up' or 'Down' relative to the ring plane.
    
    # 1. Find O
    o_idx = None
    for idx in ring_atoms:
         if mol.GetAtomWithIdx(idx).GetSymbol() == 'O': o_idx = idx; break
    if o_idx is None: return None
    
    # ... previous approach was too complex.
    pass

# There is a much simpler, highly robust approach for sugars:
# Use SMILES matching but replace the terminal OH with a generic O in the query!

print("Testing SMARTS with recursive environments...")

maltose_smi = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O[C@H]1OC(CO)[C@@H](O)[C@H](O)[C@H]1O"
mol = Chem.MolFromSmiles(maltose_smi)

for name, smi in SUGAR_SMILES.items():
    ref_mol = Chem.MolFromSmiles(smi)
    
    # Convert reference mol to a query mol where any alcohol O can be an ether O
    query = Chem.rdChemReactions.ReactionFromSmarts(
        "[C:1]-[OH:2]>>[C:1]-[O:2]"
    ).RunReactants((ref_mol,))
    # Reactions are complicated. Let's just use string replacement on SMILES to SMARTS.
    
    # A sugar is basically its stereocenters. 
    # If we convert the SMILES strings: "OC[C@H]1O[C@@H](O)" ...
    # Any "O)" can be "O*)"
    # Any "OC[" can be "*OC["
    # Any "]O" can be "]*O"
    pass

# Approach 3: 3D Conformer overlay
def test_3d_overlay():
    mol_3d = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_3d, randomSeed=42)
    # 3D overlay is too slow for 100k database.
    pass

# Approach 4: The definitive solution for CIP inversion: 
# "Virtual Hydrolysis" done correctly
def isolate_sugar_ring(mol, ring_atoms):
    emol = Chem.EditableMol(mol)
    
    # We want to replace any bond from an exocyclic O/N/S to the REST of the molecule with a bond to H.
    # Actually, RDKit lets us use `ReplaceSubstructs` or `FragmentOnBonds`.
    
    bonds_to_cut = []
    
    # Identify C6
    c6_idx = None
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in ring_atoms and nbr.GetSymbol() == 'C':
                for nnbr in nbr.GetNeighbors():
                    if nnbr.GetSymbol() in ('O', 'N', 'S'):
                        c6_idx = nbr.GetIdx()
                        break
                        
    # Find all bonds connecting the sugar to the aglycone/other sugars
    core_atoms = set(ring_atoms)
    if c6_idx is not None: core_atoms.add(c6_idx)
    
    for idx in list(core_atoms):
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in core_atoms and nbr.GetSymbol() in ('O', 'N', 'S'):
                core_atoms.add(nbr.GetIdx())
                
    # Now, any bond from core_atoms to outside should be cut!
    for idx in list(core_atoms):
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in core_atoms:
                bond = mol.GetBondBetweenAtoms(idx, nbr.GetIdx())
                if bond:
                    bonds_to_cut.append(bond.GetIdx())
                    
    if not bonds_to_cut:
        # It's already a monomer
        return Chem.Mol(mol)
        
    # Fragment on these bonds
    # This adds dummy atoms (*) where it cuts.
    frag_mol = Chem.FragmentOnBonds(mol, bonds_to_cut)
    
    # The sugar is the chunk containing the original ring atoms
    frags = Chem.GetMolFrags(frag_mol, asMols=True, frags=None, sanitizeFrags=False)
    
    # Find the fragment with the ring
    for frag in frags:
        # A fragment might not retain the original global indices, 
        # But we can check if it has a ring of size 5 or 6 with exactly one Oxygen
        ri = frag.GetRingInfo()
        valid = False
        for ring in ri.AtomRings():
            atoms = [frag.GetAtomWithIdx(a).GetSymbol() for a in ring]
            if len(ring) in (5,6) and atoms.count('O') == 1:
                valid = True
                break
                
        if valid:
            # We found the sugar fragment!
            # Replace all dummy atoms (*) with Hydrogens (H) mapping to atomic number 1
            efrag = Chem.EditableMol(frag)
            for atom in list(frag.GetAtoms())[::-1]: # reverse to not mess up indices
                 if atom.GetAtomicNum() == 0:
                     # It's a dummy '*'
                     atom_idx = atom.GetIdx()
                     # In an editable mol we can't easily change atomic num.
                     pass 
            # Actually, RDKit can do ReplaceSubstructs
            params = Chem.SmilesParserParams()
            params.sanitize = False
            dummy = Chem.MolFromSmiles("*", params)
            h_atom = Chem.MolFromSmiles("[H]", params)
            
            clean_frag = Chem.ReplaceSubstructs(frag, dummy, h_atom, replaceAll=True)[0]
            Chem.SanitizeMol(clean_frag)
            Chem.AssignStereochemistry(clean_frag, force=True, cleanIt=True)
            return clean_frag
            
    return None

import lib.monosaccharide_identifier as seq

print("Testing Virtual Hydrolysis Approach 4:")
clean_rings = []
ri = mol.GetRingInfo()
for ring in ri.AtomRings():
     atoms = [mol.GetAtomWithIdx(a).GetSymbol() for a in ring]
     if atoms.count('O') == 1:
          clean = isolate_sugar_ring(mol, list(ring))
          if clean:
              clean_rings.append(clean)
              
for i, mol_m in enumerate(clean_rings):
     print(f"Monomer {i} SMILES: {Chem.MolToSmiles(mol_m, isomericSmiles=True)}")
     # Now use our existing RS logic on this CLEAN monomer!
     ring_m = mol_m.GetRingInfo().AtomRings()[0]
     sig = seq.get_rs_signature_core(mol_m, list(ring_m))
     print(f"  RS: {sig}")
     if sig in seq.RS_LIBRARY:
         print(f"  Matched: {seq.RS_LIBRARY[sig]}")
