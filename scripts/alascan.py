from subprocess import check_output
import shlex
import math
import re
import sys
import os
import argparse

sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import IO
import settings as s
import functions as f

class Run(object):
    """
    Create FEP files from a common substructure for a given set of
    ligands
    """
    def __init__(self, lig, dist, *args, **kwargs):
        """
        
        """
        self.lig = lig
        self.dist = dist
        self.mapping = {}
        self.prot_coord = {}
        self.include = ('ATOM', 'HETATM')
        self.CYS = []
        self.mutations = []
        self.liglist = []
        self.exclude = ['PRO','GLY']
        if self.lig[0] == '[':
            self.lig = self.lig.strip('[')
            self.lig = self.lig.strip(']')
            self.lig = self.lig.split(',')
        
        if not os.path.exists('protPREP.log'):
            print("Please prepare your protein with protprep.py first")
            exit()
    
    def readlog(self):
        block = 0
        with open('protPREP.log') as infile:
            for line in infile:
                line = line.split()
                if 'mapping' in line:
                    block = 1
                    
                if 'S-S' in line:
                    block = 2 
                    
                if block == 1 and len(line) == 4:
                    try:
                        self.mapping[int(line[0])] = int(line[1])
                    except:
                        continue
                    
                if block == 2 and len(line) == 4:
                    try:
                        self.CYS.append(int(line[0]))
                        self.CYS.append(int(line[1]))
                    except:
                        continue
                        
    def get_mutations(self):
        test = []
        with open('protein.pdb') as infile:
            for line in infile:
                if line.startswith(self.include):
                    line = IO.pdb_parse_in(line)
                    if int(line[6]) in self.CYS:
                        continue

                    if line[4] in self.exclude:
                        continue
                        
                    coord = (line[8],line[9],line[10])
                    self.prot_coord[line[1]] = [coord, line[6], line[4]]

        for lig in self.lig:
            self.liglist.append("'{}'".format(lig))
            with open(lig + '.pdb') as infile:
                for line in infile:
                    if line.startswith(self.include):
                        line = IO.pdb_parse_in(line)
                        coord1 = (line[8],line[9],line[10])
                        for at in self.prot_coord:
                            coord2 = self.prot_coord[at][0]
                            if f.euclidian_overlap(coord1, 
                                                   coord2, 
                                                   float(self.dist)) \
                            == True:
                                res = self.prot_coord[at][1]
                                if self.prot_coord[at][2] != 'ALA':
                                    mutation = self.prot_coord[at][2] + str(self.mapping[res]) + 'A'
                                else:
                                    mutation = self.prot_coord[at][2] + str(self.mapping[res]) + 'G'
                                if mutation not in self.mutations:
                                    self.mutations.append(mutation)
                                    
    def write_script(self):
        with open(s.INPUT_DIR + '/alaSCAN.py') as infile, \
             open('write_alascan.py', 'w') as outfile:
            for line in infile:
                if line.strip() == 'MUTATIONS':
                    for mutation in self.mutations:
                        outfile.write("{:15s}'{}',\n".format(' ',mutation))
                    continue
                if line.strip() == 'SYSTEMS':
                    ligs = ','.join(self.liglist)
                    outfile.write("systems=[{}]\n".format(ligs))
                    continue
                if line.strip() == 'QLIGFEP':
                    outfile.write("qligfep='{}'\n".format(os.path.dirname(os.path.dirname(__file__))))
                    continue
                outfile.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='ala-scan',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Generate QresFEP inputscript for Ala scan protocol == ')

    
    parser.add_argument('-l', '--ligand',
                        nargs='*',
                        dest = "lig",
                        required = True,
                        help = "name of the ligand, or [ligand1, ligand2]")
    
    parser.add_argument('-d', '--distance',
                        dest = "dist",
                        required = True,
                        help = "distance of residues to select")

    args = parser.parse_args()
    run = Run(lig = args.lig,
              dist = args.dist
             )
    
    run.readlog()
    run.get_mutations()
    run.write_script()
