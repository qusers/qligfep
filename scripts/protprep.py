import glob
import os, sys
import shutil
import shlex

path = os.getcwd()
pymemdyn = path+'/../pymemdyn'
structures = sorted(glob.glob('*/'))

# Funnction to retrieve the COG of a ligand.
def COG(XXXX):
    os.chdir(XXXX+'1.binding/1.ligprep')
    file = open('COG.txt', 'r')
    xyz = file.readline().strip('\n')
    os.chdir(path)
    return xyz

# Function to retrieve the system.pdb file for a structure.
def pdb(XXXX):
    os.chdir('{}/{}/eq/'.format(pymemdyn, XXXX))
    pdb = os.path.abspath('system.pdb')
    os.chdir(path)
    return pdb

# Neutralized the simulation sphere by (de)protonation of residues
def neutralize(XXXX, xyz, charge):
    os.chdir(XXXX+'1.binding/2.protprep/')
    os.system('python /home/koenekoop/software/qligfep/scripts/neutralize_sphere.py -p protein.pdb -c '+xyz+' -r 25 -T +'+str(charge))

# For every protein structure, preps the protein PyMEMdyn system.pdb file with respect to the COG of the corresponding ligand and generates protein.pdb and water,pdb files.
for XXXX in structures:
    print(XXXX.strip('/'))
    xyz = COG(XXXX)
    protein = pdb(XXXX)

    os.chdir(XXXX+'1.binding/2.protprep')
    os.system('python3 /home/jespers/software/qligfep/protprep.py -p '+protein+' -r 25 -c '+xyz+' --noclean -O gromacs')

#    neutralize(XXXX, xyz, charge[XXXX])
    
    os.chdir(path)