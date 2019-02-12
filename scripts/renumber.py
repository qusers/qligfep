import sys

installdir = '/home/jespers/software/qligfep/' 
sys.path.append(installdir)

import IO
import functions

pdbfile = '6MEO_full.pdb'
resn_ref = 1
resn_out = 1

with open(pdbfile) as infile, open(pdbfile[:-4] + '_renumber.pdb', 'w') as outfile:
    for line in infile:
        if line.startswith('TER'):
            outfile.write(line)
        elif line.startswith('ATOM'):
            line2 = IO.pdb_parse_in(line)
            resn = line2[6]
            if resn != resn_ref:
                resn_out += 1
                resn_ref = resn

            line2[6] = resn_out
            
            outfile.write(IO.pdb_parse_out(line2) + '\n')
