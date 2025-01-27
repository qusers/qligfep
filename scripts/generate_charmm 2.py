import glob
from subprocess import check_output
import shlex
import os
import shutil

molecules = glob.glob('*.mol2')
# FIX THIS
os.system('module load openbabel/2.4.1')
cgenff = '/home/jespers/software/cgenff/cgenff'

def readprm(prmfile):
    block = 0
    nonbonded = {}
    with open(prmfile) as infile:
        for line in infile:
            line1 = line.split()
            if len(line1) < 1:
                continue
            if line1[0] == 'ATOMS':
                block = 1

            if line1[0] == 'NONBONDED':
                block = 2
            
            if block == 1:
                if line1[0] == 'MASS':
                    nonbonded[line1[2]] = [line]

            if block == 2:
                if len(line1) < 4:
                    continue
                try:
                    nonbonded[line1[0]].append(line)
                except:
                    continue

    return nonbonded

def CGenFF(molecules, nonbonded):
    for molecule in molecules:
        block = 0
        lig = molecule.split('.')[0]
        atomtypes = []
        call = '/home/jespers/software/cgenff/cgenff {} -a > {}.str'.format(molecule, lig)
        os.system(call)
        str_file = '{}.str'.format(lig)
        prm_file = '{}.prm'.format(lig)
        prm_file2 = '{}.merge.prm'.format(lig)
        tpr_file  = '{}.rtf'.format(lig) 
        with open(str_file) as infile,         \
             open(prm_file, 'w') as outfile1,  \
             open(tpr_file, 'w') as outfile2:
            for line in infile:
                if line == 'read rtf card append\n':
                    block = 1
                    continue

                if line == 'read param card flex append\n':
                    block = 2
                    continue

                if block == 1:
                    line_split = line.split()
                    if len(line_split) < 1:
                        continue
                    if line_split[0] == 'ATOM':
                        if line_split[2] not in atomtypes:
                            atomtypes.append(line_split[2])
                    outfile2.write(line)

                if block == 2:
                    outfile1.write(line)
        
        with open(prm_file) as infile,     \
             open(prm_file2, 'w') as outfile:
                
            for line in infile:
                outfile.write(line)
                if len(line.split()) > 0:
                    if line.split()[-1] == 'validation/optimization.':
                        outfile.write('\n')
                        for line in atomtypes:
                            outline = nonbonded[line][0]
                            outfile.write(outline)
                    if line.split()[-1].strip() == 'RETURN':
                        outfile.write('NONBONDED\n')
                        outfile.write('cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5\n') # NOTE: this needs to be dehardcoded!!
                        for line in atomtypes:
                            outfile.write(nonbonded[line][1])

def mol2pdb(molecules):
    for molecule in molecules:
        lig = molecule.split('.')[0]
        call = 'babel -imol2 {} -opdb {}.pdb'.format(molecule, lig)
        os.system(call)

def charmm2Q(molecules):
    for molecule in molecules:
        lig = molecule.split('.')[0]
        shutil.move('{}.merge.prm'.format(lig), '{}.prm'.format(lig))
        call = 'python  /home/jespers/software/setupFEP/charmm2Q.py -l {} -FF CHARMM36'.format(lig)
        os.system(call)

def movefiles(molecules):
    for molecule in molecules:
        lig = molecule.split('.')[0]
        for extension in ['.pdb', '.lib', '.prm']:
            if extension == '.pdb':
                src = lig + extension
                tgt = lig + extension
            else:
                src = 'Q' + lig + extension
                tgt = lig + extension
            
            shutil.move(src, tgt)

nonbonded = readprm('/home/jespers/software/cgenff/par_all36_cgenff.prm')
CGenFF(molecules, nonbonded)
mol2pdb(molecules)
charmm2Q(molecules)
movefiles(molecules)
