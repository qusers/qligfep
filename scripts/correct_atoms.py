import argparse
import glob
import os, sys
import shutil
import subprocess
sys.path.append('/home/koenekoop/software/qligfep/')
import IO

with open('setup/hoh_ref.pdb', 'r') as infile, open('pymemdyn/hoh.pdb', 'w') as outfile:
    HW = ['H01', 'H02', 'H03', 'H04', 'H05']
    for line in infile:
        line = IO.pdb_parse_in(line)
        if line[2] == 'O' and line[3] != 'B':
            line[2] = 'OW'
            line[3] = ''
            outfile.write(IO.pdb_parse_out(line))
            outfile.write('\n')
        elif line[2] in HW and line[3] != 'B':
            number = line[2][-1]
            if number == '1':
                line[2] = 'HW1'
                line[3] = ''
                outfile.write(IO.pdb_parse_out(line))
                outfile.write('\n')
            else:
                line[2] = 'HW2'
                line[3] = ''
                outfile.write(IO.pdb_parse_out(line))
                outfile.write('\n')
    outfile.close()
    
with open('setup/ions_local_ref.pdb', 'r') as infile, open('pymemdyn/ions_local.pdb', 'w') as outfile:
    for line in infile:
        if line == 'END\n':
            continue
        line = IO.pdb_parse_in(line)
        ion = line[2]
        if line[4] != '{}+'.format(ion):
            line[4] = '{}+'.format(ion)
        if line[-2] != ion:
            line[-2] = ion
        if line[-1] != '':
            line[-1] = ''
        line[5] = ''
        line[6] = 1
        outfile.write(IO.pdb_parse_out(line))
        outfile.write('\n')
    infile.close()
    outfile.close()