import argparse
import glob
import os, sys
import shutil
import subprocess
sys.path.append('/home/koenekoop/software/qligfep/')
import IO


path = os.getcwd()
XXXX = sorted(glob.glob('*/'))

for i in XXXX:
    os.chdir(i+'eq/')
    os.system('echo 1 0 | trjconv -pbc mol -center -ur compact -f confout.gro -o confout.pdb')
    pdb = os.path.abspath("confout.pdb")
    
    with open(pdb) as infile, open('system.pdb', 'w') as outfile:
        for line in infile:
            inline = IO.pdb_parse_in(line)
            if inline[0] == 'ATOM  ':
                if inline[4] == 'HOH':
                    inline[4] = 'SOL'
                elif inline[4] == 'NA+':
                    inline[2] = 'SOD'
                    inline[4] = 'SOD'
                outline = IO.pdb_parse_out(inline) + '\n'
                outfile.write(outline)
            else:
                outfile.write(line)
    os.chdir(path)
    