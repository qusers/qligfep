import glob
import os, sys
import shutil
from shutil import copyfile

path = os.getcwd()
structures = sorted(glob.glob('*/'))

for XXXX in structures:
    print(XXXX.strip('/'))
    # Copies the ligand files to the setupQLIE folder (protein simulation) and the solvation folder (water simulation).
    os.chdir(XXXX+'1.binding/1.ligprep')
    ligfiles = glob.glob('*.pdb') + glob.glob('*.prm') + glob.glob('*.lib') + glob.glob('*.sdf')
    for file in ligfiles:
        shutil.copy('./'+file, './../3.setupQLIE/'+file)
        shutil.copy('./'+file, './../../2.solvation/'+file)

    # Copies the protein files to the setupQLIE folder (protein simulation).
    os.chdir('../2.protprep')
    shutil.copy('./protein.pdb', './../3.setupQLIE/protein.pdb')
    shutil.copy('./water.pdb', './../3.setupQLIE/water.pdb')
    shutil.copy('./protPREP.log', './../3.setupQLIE/protPREP.log')
    
    # Sets up the QLIE protein simulation job - and submits the job shell script.
    os.chdir('../3.setupQLIE')
    ligand = glob.glob('*.pdb')[0].split('.pdb')[0]
    os.system('python /home/koenekoop/software/qligfep/QLIE_conv.py -l '+ligand+' -S protein -P CSB -C TETRA -T 298 -r 10 -f OPLS2015')

    os.chdir('LIE_'+ligand)
    os.system('./FEP_submit.sh')
    
    # Sets up the QLIE water simulation job
    os.chdir('../../../2.solvation')
    
    os.system('python /home/koenekoop/software/qligfep/QLIE_conv.py -l '+ligand+' -S water -P CSB -C TETRA -T 298 -r 10 -f OPLS2015 -R 25')
    
    os.chdir('LIE_'+ligand)
    os.system('./FEP_submit.sh')
    
    os.chdir(path)