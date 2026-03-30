"""
Test isolating a sugar ring and hydrolyzing its attachments to restore monomer CIP priorities.
"""
import sys, os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__))))
from rdkit import Chem
from lib.monosaccharide_identifier import get_rs_signature_core

def isolate_and_hydrolyze(mol, ring_atoms):
    """
    Extracts the ring atoms and immediate exocyclic oxygens/carbons,
    capping any broken bonds to restore the monomer state.
    """
    # Create an editable molecule
    emol = Chem.EditableMol(mol)
    
    # We want to keep:
    # 1. Ring atoms
    # 2. Exocyclic atoms (O, N, or C for C6) directly attached to ring atoms
    # 3. For C6, the O attached to it.
    
    atoms_to_keep = set(ring_atoms)
    
    # Identify C6
    c6_idx = None
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in ring_atoms and nbr.GetSymbol() == 'C':
                # Check if this C is attached to an O (typical C6)
                for nnbr in nbr.GetNeighbors():
                    if nnbr.GetSymbol() in ('O', 'N', 'S'):
                        c6_idx = nbr.GetIdx()
                        atoms_to_keep.add(c6_idx)
                        atoms_to_keep.add(nnbr.GetIdx())
                        break
    
    # Add other direct exocyclic heteroatoms
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in ring_atoms and nbr.GetSymbol() in ('O', 'N', 'S'):
                atoms_to_keep.add(nbr.GetIdx())
                
    # Create a submol with just these atoms
    sub_mol = Chem.PathToSubmol(mol, list(atoms_to_keep))
    
    # The submol might still have wildcards/dummy atoms where we cut things.
    # In RDKit, PathToSubmol doesn't add dummy atoms by default, it just leaves valences open.
    # We should sanitize it to add implicit hydrogens where bonds were cut.
    try:
        Chem.SanitizeMol(sub_mol)
    except:
        pass
        
    return sub_mol

maltose = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O[C@H]1OC(CO)[C@@H](O)[C@H](O)[C@H]1O"
mol = Chem.MolFromSmiles(maltose)
Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

ri = mol.GetRingInfo()
for i, ring in enumerate(ri.AtomRings()):
    print(f"\n--- Ring {i} ---")
    sub_mol = isolate_and_hydrolyze(mol, ring)
    print(f"Submol SMILES: {Chem.MolToSmiles(sub_mol)}")
    
    # We need the new indices of the ring atoms in the submol
    # PathToSubmol preserves the order of the requested atoms: No it doesn't always.
    # We can map them by checking atom parity or just passing a mapped ring.
    
    # A better way to get the RS signature of a submol:
    # RDKit's FindMolChiralCenters gives us the (idx, RS) of the submol.
    Chem.AssignStereochemistry(sub_mol, force=True, cleanIt=True)
    centers = Chem.FindMolChiralCenters(sub_mol, includeUnassigned=True)
    print(f"Submol centers: {centers}")
    
