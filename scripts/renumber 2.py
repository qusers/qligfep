#!/usr/bin/env python3

# Script to renumber residues in pdb files starting at given index.
# Original by W. Jespers
# Modifications by F.W van der Ent
# 2019-09-18

import sys
installdir = '/home/vanderent/software/qligfep/' 
sys.path.append(installdir)

import IO
import functions
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='input', help='Path to input pdb')
parser.add_argument('start', type=int, help='Starting point for residue renumbering')
parser.add_argument('--overwrite', dest='overwrite', action='store_true')
parser.set_defaults(overwrite=False)

args = parser.parse_args()

pdbfile = args.input
resn_ref = args.start
resn_out = args.start

with open(pdbfile) as infile, open(pdbfile[:-4] + '_renumber.pdb', 'w') as outfile:
    for line in infile:
        if line.startswith('TER'):
            outfile.write(line)
        elif 'SPHERE' in line:
            outfile.write(line)
        elif 'GAP' in line:
            outfile.write(line)
        elif line.startswith('ATOM'):
            line2 = IO.pdb_parse_in(line)
            resn = line2[6]
            if resn != resn_ref:
                resn_out += 1
                resn_ref = resn

            line2[6] = resn_out
            
            outfile.write(IO.pdb_parse_out(line2) + '\n')
