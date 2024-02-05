import glob
import statistics
# import pandas
import argparse
import os
import shutil
import re
import functions as f
import settings as s
import IO
import pymol
from pymol import cmd

class Run(object):
    """
    """
    def __init__(self, FEP, INPUTS, FF, convert_files, *args, **kwargs):
        self.FEP = FEP.strip('/')
        self.INPUTS = INPUTS.strip('/')
        self.dihedralangle1 = []
        self.dihedralangle2 = []
        self.rootdir = os.getcwd() + "/"
        self.FF = FF
        self.convert_files = convert_files

    def replace(self, string, replacements):
        pattern = re.compile(r'\b(' + '|'.join(replacements.keys()) + r')\b')
        replaced_string = pattern.sub(lambda x: replacements[x.group()], string)
        return replaced_string

    def qprep(self, writedir):
        os.chdir(writedir)
        # Somehow Q is very annoying with this < > input style so had to implement
        # another function that just calls os.system instead of using the preferred
        # subprocess module....
        IO.run_command('/home/wjespers/software/Q/bin/qprep', ' <convert_re2pdb_temp.inp > convert_re2pdb_temp.out', string = True)
        os.chdir(self.rootdir)
        
    def make_pdb_files(self):
        replacements = {}
        repfolders = sorted(glob.glob(self.FEP + '/*/'))
        directories = [repfolders[0] + entry for entry in os.listdir(repfolders[0]) if os.path.isdir(repfolders[0] + entry)]
        directories = sorted(directories)
        inputdirectory = self.rootdir + self.INPUTS 
        top_file = glob.glob(inputdirectory + '/*.top')[0]
        libfile_file = glob.glob(inputdirectory + '/*.lib')[0]
        replacements['TOP_FILE_PATH'] = top_file
        replacements['LIB_FILE_PATH'] = libfile_file
        replacements["FF_FILE_PATH"] = self.rootdir + 'FF/' + self.FF + '.lib'
        for folder in directories:
            qfepfilesrep = sorted(glob.glob(folder + '/*.re'))
            qfepfilesrep = [file_ for file_ in qfepfilesrep if "md"in file_]
            qfepfilesrep = sorted(qfepfilesrep)
            for file in qfepfilesrep:
                if self.convert_files == True:
                    shutil.copy("INPUTS/convert_re2pdb.inp", self.rootdir + folder)
                    full_path = self.rootdir + file
                    replacements['RE_INPUT'] = file.split("/")[-1]
                    replacements['PDB_OUPUT'] = file.split("/")[-1].strip('.re') + ".pdb"
                    run.replace_placeholders(replacements, self.rootdir + folder + "/convert_re2pdb.inp",  self.rootdir + folder + "/convert_re2pdb_temp.inp")
                    writedirectory = self.rootdir + folder
                    run.qprep(writedirectory)
                filepath_created_pdb = file.strip('.re') + ".pdb"
                cmd.load(filepath_created_pdb,"ethylbenzene")
        no_of_states = cmd.count_states("ethylbenzene")
        print(no_of_states)
        for state_num in range(1, no_of_states + 1):
            # first dihedral
            print(cmd.get_dihedral(atom1="resi 111 and name N",atom2="resi 111 and name CA",atom3="resi 111 and name CB",atom4="resi 111 and name CG1",state=state_num))

            # second dihedral
            #print(cmd.get_dihedral(atom1="resi 111 and name N",atom2="resi 111 and name CA",atom3="resi 111 and name CB",atom4="resi 111 and name CG2",state=state_num))
    
    def replace_placeholders(self, replacements, convertion_Q_file_path, convertion_Q_output):
        with open(convertion_Q_file_path) as infile, open(convertion_Q_output, 'w') as outfile:
            for line in infile:
                line = run.replace(line, replacements)
                outfile.write(line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='retreive_angles',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Analyse FEP == ')
    
    parser.add_argument('-F', '--FEP',
                        dest = "FEP",
                        required = True,
                        help = "name of FEP directory (FEP_$)")

    parser.add_argument('-I', '--INPUTS',
                    dest = "INPUTS",
                    required = True,
                    help = "top and lib file location")

    parser.add_argument('-FF', '--forcefield',
                        dest = "FF",
                        required = True,
                        choices = ['OPLS2005', 'OPLS2015', 'AMBER14sb', 'CHARMM36', 'CHARMM22', 'CHARMM_TEST'],
                        help = "Forcefield to be used")

    parser.add_argument('-convert', '--convert_files',
                        dest = "convert_files",
                        required = False,
                        action='store_true',
                        help = "Convert .re to .pdb files")

    args = parser.parse_args()
    run = Run(FEP = args.FEP,
              INPUTS = args.INPUTS,
              FF = args.FF,
              convert_files = args.convert_files)

    run.make_pdb_files()
    
