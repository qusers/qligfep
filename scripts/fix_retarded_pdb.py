import glob
import os, sys
from subprocess import check_output
import shutil
import shlex
sys.path.append('/home/koenekoop/software/qligfep/') # Make sure it can import form parent directory instead of hardcode path
import IO
import argparse

def pdb_parse_out(line):
    """
    Takes a list and parses it into a pdb writeable line
    """
    if len(line[2]) <= 3: 
        line = '{:6s}{:5d}  {:3s}{:1s}{:4s}{:1s}{:4d}{:1s}   '\
               '{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format(*line)

    elif len(line[2]) == 4: 
        line = '{:6s}{:5d}  {:4s}{:s}{:4s}{:1s}{:4d}{:1s}   '\
               '{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format(*line)
    return line

class Run(object):
    """
    Create FEP files from a common substructure for a given set of
    ligands
    """
    def __init__(self, ipdb, *args, **kwargs):
        self.ipdb = ipdb

    def Fix(self, include='ATOM,HETATM'):

        # Defines the ATOM and HETATM start of pdb files
        include = ('ATOM','HETATM')

        # Reads all the ligands in your directory and saves them in a list
        ligand = self.ipdb

        # Opens the ligand file to read the lines from as infile, and opens a temporary pdb file to write the corrections to as outfile
        with open(ligand) as infile, \
            open("tmp.pdb", "w") as outfile:
            # loops over all the lines of your pdb file
            for line in infile:
                # Makes sure to only do this for ATOM and HETATM lines
                if line.startswith(include):
                    # Converts the pdb lines to a list with the function pdb_parse_in() from IO.py in qligfep
                    line = IO.pdb_parse_in(line)
                    if len(line[2]) == 4:
                        line[2] = line[2][1:]+line[2][0]
                        line[3] =''
                    # converts the line list back to pdb line
                    outline = pdb_parse_out(line)
                    # Writes the outline to the temporary pdb file
                    outfile.write(outline + "\n")
            infile.close()
            outfile.close()
    
        # Opens the temporary pdb file and rewrites the original pdb file
        with open("tmp.pdb") as infile, \
            open(ligand, "w") as outfile:
                for line in infile:
                    outfile.write(line)
        infile.close()
        outfile.close()

        # removes the temporary pdb file
        os.remove("tmp.pdb")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='fix_atomnames',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == realigns the atom name column for Q of a given .pdb == ')


    parser.add_argument('-i', '--ipdb',
                        dest = "ipdb",
                        required = True,
                        help = "name of pdbfile")

    args = parser.parse_args()
    run = Run(ipdb = args.ipdb)
            
