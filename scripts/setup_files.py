import glob
import os, sys
import shutil
import subprocess

structures = glob.glob('setup/CrystStructs/*.pdb')

for XXXX in structures:
    XXXX = XXXX.strip('.pdb')[-4:]
    os.mkdir('QLIE/{}'.format(XXXX))
    os.mkdir('QLIE/{}/1.binding'.format(XXXX))
    os.mkdir('QLIE/{}/2.solvation'.format(XXXX))
    os.mkdir('QLIE/{}/1.binding/1.ligprep'.format(XXXX))
    os.mkdir('QLIE/{}/1.binding/2.protprep'.format(XXXX))
    os.mkdir('QLIE/{}/1.binding/3.setupQLIE'.format(XXXX))
    os.mkdir('pymemdyn/{}'.format(XXXX))