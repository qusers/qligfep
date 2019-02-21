import glob
import sys
import os

sys.path.append('/home/jespers/software/setupFEP')

import IO
import settings as s

replicates = 10
replacements = {'PROT_START':'1',
                'PROT_END':'257',
                'Q_ATOMS':'residue 258'
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

def run_qcalc():
    for i in range(1, replicates + 1):
        curdir = os.getcwd()
        os.chdir('FEP1/298/1'.format(i))
        dcd = glob.glob('md*.dcd')
        #dcd = sorted(dcd)
        #print dcd
        write_qcalc(dcd)
        os.system('qcalc < qcalc_tmp.inp > qcalc.log')
        os.chdir(curdir)

def print_results():
    with open('FEP1/298/1/qcalc.log') as infile:
        for line in infile:
            line=line.split()
            if len(line) == 3:
                try:
                    if abs(float(line[2])) > 1.0 or abs(float(line[1])) > 1.0: 
                        print line[0], line[1], line[2]
                except:
                    continue

run_qcalc()
print_results()
