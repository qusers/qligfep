from rdkit import Chem
from rdkit.Chem import AllChem, rdmolfiles, MCS, rdFMCS

def align_molecules(list_smiles, iterations, names):
    '''
        Get maximum common substructure (MCS) from a list of molecules. The first molecule should be the most simplistic of all molecules.
        The function will try to reference all coordinates to this specific molecule.
    '''
    mols = [AllChem.MolFromSmiles(refmol) for refmol in list_smiles]
    mols_with_H = [Chem.AddHs(mol) for mol in mols]
    [AllChem.EmbedMolecule(x,AllChem.ETKDG()) for x in mols_with_H]

    # Find the MCS:
    mcs = rdFMCS.FindMCS(mols_with_H,threshold=0.8,completeRingsOnly=True,ringMatchesRingOnly=True)

    # align everything to the first molecule
    patt = Chem.MolFromSmarts(mcs.smartsString)
    refMol = mols_with_H[0]
    refMatch = refMol.GetSubstructMatch(patt)

    mol_edit = Chem.EditableMol(refMol)
    for atom in refMol.GetAtoms():
        id_atom = atom.GetIdx()
        if id_atom not in refMatch:
            mol_edit.RemoveAtom(id_atom)
    mol_to_allign_to = mol_edit.GetMol()
    
    minimized_rmsd = []
    aligned_minimized_mols = []
    for probeMol in list_smiles:
        aligned_mols = []
        rms_aligned_mols = []
        for i in range(iterations):
            Mol = AllChem.MolFromSmiles(probeMol)
            Mol = Chem.AddHs(Mol)
            AllChem.EmbedMolecule(Mol,AllChem.ETKDG())
            
            # get atoms that match substructure
            smarts_mcs = Mol.GetSubstructMatch(patt)
            rms = AllChem.AlignMol(Mol,mol_to_allign_to,atomMap=list(zip(smarts_mcs,refMatch)))
            rms_aligned_mols.append(rms)
            aligned_mols.append(Mol)

        minimum_rmsd = min(rms_aligned_mols)
        index_lowest_rmsd = rms_aligned_mols.index(minimum_rmsd)
        lowest_rmsd_mol = aligned_mols[index_lowest_rmsd]
        minimized_rmsd.append(minimum_rmsd)
        aligned_minimized_mols.append(lowest_rmsd_mol)

    for compound_number, name in enumerate(names):
        compound = aligned_minimized_mols[compound_number]
        Chem.rdmolfiles.MolToMolFile(compound, name + '.mol')
        
    return minimized_rmsd,aligned_minimized_mols

