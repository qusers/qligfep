import glob
import os, sys
from subprocess import check_output
import shutil
import shlex
sys.path.append('/home/koenekoop/software/qligfep/') # Make sure it can import form parent directory instead of hardcode path
import IO

# Defines the ATOM and HETATM start of pdb files
include = ('ATOM','HETATM')

# Reads all the ligands in your directory and saves them in a list
ligands = glob.glob("*.pdb")

# Loops over all the ligand files in your ligands list
for ligand in ligands:
    # Opens the ligand file to read the lines from as infile, and opens a temporary pdb file to write the corrections to as outfile
    with open(ligand) as infile, \
        open("tmp.pdb", "w") as outfile:
        # loops over all the lines of your pdb file
        for line in infile:
            # Makes sure to only do this for ATOM and HETATM lines
            if line.startswith(include):
                # Converts the pdb lines to a list with the function pdb_parse_in() from IO.py in qligfep
                line = IO.pdb_parse_in(line)
                # fixes the retarded atom name
                line[2] = line[2][:1]+str(line[1])
                # converts the line list back to pdb line
                outline = IO.pdb_parse_out(line)
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
            