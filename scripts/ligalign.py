import glob
import os, sys
import shutil
import pymol
from pymol import cmd

path = os.getcwd()
pymemdyn = path+'/../pymemdyn'
structures = sorted(glob.glob('*/'))

# Returns a class object for each crystal strcutre, with ligand, pdbcode and mode of action attributes.
def YYY(XXXX):
    
    # Retrieves the ligand file.
    os.chdir(XXXX+'1.binding/1.ligprep')
    LigStruct = glob.glob('*.pdb')[0].strip('.pdb').split('_')
    
    # Define a class object for a ligand with its pdb code, its protein pdb code and its mode of action.
    class Ligand:
        def __init__(self, lig, prot):
            self.lig = lig
            self.prot = prot
    
    # Saves the ligand ass class object.
    ligand = Ligand(LigStruct[0], LigStruct[1])
    os.chdir(path)
    
    # Returns the ligand.
    return ligand

# Returns a path to the corresponding pymemdyn protein file for a crystal structure.
def Get_System(XXXX):
    system = '{}/{}/eq/system.pdb'.format(pymemdyn, XXXX)
    return system

# Creates a list of all the ligands.
Ligands = []
for XXXX in structures:
    Ligands.append(YYY(XXXX))
    
#For every ligand...
for YYY in Ligands:
    # ...fetch the corresponding crystal structure.
    cmd.fetch(YYY.prot)
    # (Delete the created .cif file right away.)
    cif = glob.glob('*.cif')
    for file in cif:
        os.remove(file)

    # ...load the corresponding pymemdyn confout.pdb file.
    cmd.load(Get_System(YYY.prot),'system_'+YYY.prot)
    
    # ...align the crystal structure on the pymemdyn confout.pdb file.
    cmd.align(YYY.prot, 'system_'+YYY.prot)
    
    # ...select only the ligand (in its correct coordination) and adds Hydrogens.
    cmd.select('resn '+YYY.lig+' and chain A')
    cmd.h_add('sele')
    cmd.select('resn '+YYY.lig+' and chain A')
    
    # ...and finally save as LIG_PROT_MODE.pdb file.
    os.chdir('{}/{}/1.binding/1.ligprep/'.format(path, YYY.prot))
    cmd.save(YYY.lig+'_'+YYY.prot+'.pdb', 'sele')
    
    # Clear the current workspace.
    cmd.reinitialize()
    os.chdir(path)