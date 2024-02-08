import glob
import os, sys
from subprocess import check_output
import shutil
import shlex
sys.path.append('/home/koenekoop/software/qligfep/') # Make sure it can import form parent directory instead of hardcode path
import IO

path = os.getcwd()
structures = sorted(glob.glob('*/'))

def COG(XXXX):
    os.chdir(XXXX+'1.binding/1.ligprep')
    ligand = glob.glob('*.pdb')[0]
    COG = IO.run_command('python /home/koenekoop/software/qligfep/scripts/COG.py ', '-i '+ligand)
    txt = open(r"COG.txt", "w")
    txt.write(COG.decode("utf-8"))
    os.chdir(path)
    return COG

def ligprep(XXXX):
    os.chdir(XXXX+'1.binding/1.ligprep')
    os.system('python /home/koenekoop/software/qligfep/scripts/generate_opls.py')
    exe = '/home/apps/schrodinger2018/run fix_bond_orders.py '
    infile = glob.glob('*.pdb')[0]
    outfile = infile.strip('.pdb')+'.sdf'
    options = infile+' '+outfile
    IO.run_command(exe, options)
    print(XXXX+': prepped')
    os.chdir(path)

for XXXX in structures:
    COG(XXXX)
    ligprep(XXXX)