import glob
import sys
import os

sys.path.append('/home/jespers/software/setupFEP')

import IO
import settings as s

replicates = 10
replacements = {'PROT_START':'1',
                'PROT_END':'328',
                'Q_ATOMS':'329 330'
               }

def write_qcalc(dcd):
    with open(s.INPUT_DIR+ '/qcalc.inp') as infile,          \
         open('qcalc_tmp.inp', 'w') as outfile:
        for line in infile:
            line = IO.replace(line, replacements)
            if line.rstrip() == 'TRAJECTORIES':
                for trajectory in dcd:
                    trajectory = trajectory + '\n'
                    outfile.write(trajectory)
                continue
                    
            outfile.write(line)

for i in range(1, replicates + 1):
    curdir = os.getcwd()
    os.chdir('FEP1/298/1'.format(i))
    dcd = glob.glob('md*.dcd')
    #dcd = sorted(dcd)
    #print dcd
    write_qcalc(dcd)
    os.system('qcalc < qcalc_tmp.inp > qcalc.log')
    os.chdir(curdir)
