import glob
import statistics
import math
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
from pymol import stored
import numpy as np
import time

class Run(object):
    """
    """
    def __init__(self, FEP, INPUTS, FF, convert_files, align_pdbs, *args, **kwargs):
        self.FEP = FEP.strip('/')
        self.INPUTS = INPUTS.strip('/')
        self.dihedralangle1 = []
        # self.dihedralangle2 = []
        self.rootdir = os.getcwd() + "/"
        self.FF = FF
        self.convert_files = convert_files
        self.align_pdbs = align_pdbs
        self.residue_list = ["108", "109", "110", "111", "112", "113"]
        self.fixed_point_residue = "107"
        self.lib_for_df_solute = {}
        self.lib_for_df_shell = {}
        self.dG_BAR = {}
        self.md_file_length = ""

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
        angles = []
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
            self.md_file_length = len(qfepfilesrep)
            #print(len(qfepfilesrep))
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
        # cmd.load("reference.pdb", "ethylbenzene")
        no_of_states = cmd.count_states("ethylbenzene")
        # cmd.intra_fit("(name ca)", no_of_states)
        print(no_of_states)
        # cmd.intra_fit()
        for state_num in range(1, no_of_states):
            # first dihedral
            angles.append(str(round(cmd.get_dihedral(atom1="resi 111 and name N",atom2="resi 111 and name CA",atom3="resi 111 and name CB",atom4="resi 111 and name CG1",state=state_num),2)))
            
            # second dihedral
            # print(cmd.get_dihedral(atom1="resi 111 and name N",atom2="resi 111 and name CA",atom3="resi 111 and name CB",atom4="resi 111 and name CG2",state=state_num))
        angles_per_replicate = list(run.chunks(angles, self.md_file_length))
        print(self.FEP.split("/")[-1])
        print("angles")
        angles_tabular_form = list(zip(*angles_per_replicate))
        for angles_per_rep in angles_tabular_form:
            print(', '.join(angles_per_rep))

    def replace_placeholders(self, replacements, convertion_Q_file_path, convertion_Q_output):
        with open(convertion_Q_file_path) as infile, open(convertion_Q_output, 'w') as outfile:
            for line in infile:
                line = run.replace(line, replacements)
                outfile.write(line)

    def align_to_reference_benzene(self):
        if self.align_pdbs == True:
        
            # cmd.reinitialize()
            # cmd.load("reference.pdb","reference")
            repfolders = sorted(glob.glob(self.FEP + '/*/'))
            directories = [repfolders[0] + entry for entry in os.listdir(repfolders[0]) if os.path.isdir(repfolders[0] + entry)]
            #start = time.time()
            for folder in directories:    
                cmd.reinitialize()
                cmd.load("reference.pdb","reference")
                qfepfilesrep = sorted(glob.glob(folder + '/md*.pdb'))
                for md_pdb_file in qfepfilesrep:
                    #cmd.reinitialize()
                    #cmd.load("reference.pdb","reference")
                    # print(md_pdb_file.strip(".pdb") + "_2.pdb")
                    # print(md_pdb_file)
                    # print(md_pdb_file.split("/")[-1].split(".")[0])
                    cmd.load(md_pdb_file, md_pdb_file.split("/")[-1].split(".")[0])
                    #cmd.align(md_pdb_file.split("/")[-1].split(".")[0], "reference")
                    #no_of_states = cmd.count_states()
                    #outfile = md_pdb_file.strip(".pdb") + "_2.pdb"
                    # print(outfile)
                    #cmd.save(outfile, md_pdb_file.split("/")[-1].split(".")[0], format="pdb")   
                cmd.alignto("reference")
                for aligned_pdbfile in qfepfilesrep:
                    outfile = aligned_pdbfile.replace(".pdb", "_2.pdb")
                    # print(outfile)
                    cmd.save(outfile, aligned_pdbfile.split("/")[-1].split(".")[0], format="pdb")
            #end = time.time()
            #print(f"Time taken to run the code was {end-start} seconds")
        else:
            pass

    def save_aligned_trajectory_pdb_in_single_state(self):
        repfolders = sorted(glob.glob(self.FEP + '/*/'))
        directories = [repfolders[0] + entry for entry in os.listdir(repfolders[0]) if os.path.isdir(repfolders[0] + entry)]
        for folder in directories:
            cmd.reinitialize()
            qfepfilesrep = sorted(glob.glob(folder + '/md*_2.pdb'))
            for md_pdb_file in qfepfilesrep:
                cmd.load(md_pdb_file, "trajectory")
            cmd.save(folder + "/trajectory.pdb", state=0, format="pdb")

    def find_RMSD_ligand(self):
        repfolders = sorted(glob.glob(self.FEP + '/*/'))
        directories = [repfolders[0] + entry for entry in os.listdir(repfolders[0]) if os.path.isdir(repfolders[0] + entry)]
        ligand_starting_coordinates = []
        RMSD_values = []
        directory_ligand_in_list = self.FEP.split("/")[:-1]
        directory_ligand = "/".join(directory_ligand_in_list)
        with open(directory_ligand + "/inputfiles/complexnotexcluded.pdb") as starting_structure:
            for line in starting_structure:
                if "LIG" in line:
                    if "C" in line:
                        x_coordinate = float(list(filter(None,line.split(" ")))[5])
                        y_coordinate = float(list(filter(None,line.split(" ")))[6])
                        z_coordinate = float(list(filter(None,line.split(" ")))[7])
                        coordinates = (x_coordinate, y_coordinate, z_coordinate)
                        ligand_starting_coordinates.append(coordinates)
        for folder in directories:
            qfepfilesrep = sorted(glob.glob(folder + '/md*_2.pdb'))
            for aligned_pdb_file in qfepfilesrep:
                md_coordinates = []
                with open(aligned_pdb_file) as aligned_pdbfile:
                    for line in aligned_pdbfile:
                        if "LIG" in line:
                            if "C" in line:
                                x_coordinate_md = float(list(filter(None,line.split(" ")))[5])
                                y_coordinate_md = float(list(filter(None,line.split(" ")))[6])
                                z_coordinate_md = float(list(filter(None,line.split(" ")))[7])
                                coordinates_md = (x_coordinate_md, y_coordinate_md, z_coordinate_md)
                                md_coordinates.append(coordinates_md)
                #print(md_coordinates)
                #print(str((run.rmsd(ligand_starting_coordinates, md_coordinates))))
                RMSD_values.append(str(round((run.rmsd(ligand_starting_coordinates, md_coordinates)),2)))
        #self.md_file_length = 112   
        rmsd_per_replicate = list(run.chunks(RMSD_values, self.md_file_length))
        #print(len(rmsd_per_replicate))
        rmsd_tabular_form = list(zip(*rmsd_per_replicate))
        print("RMSD")
        for rmsd_per_rep in rmsd_tabular_form:
            print(', '.join(rmsd_per_rep))

    def rmsd(self, allcoordsA, allcoordsB):
        deviation = sum(run.squared_distance(atomA, atomB) for 
                        (atomA, atomB) in zip(allcoordsA, allcoordsB))
        return math.sqrt(deviation / float(len(allcoordsA)))

    def squared_distance(self, coordsA, coordsB):
        sqrdist = sum( (a-b)**2 for a, b in zip(coordsA, coordsB) )
        return sqrdist

    def return_restraints(self):
        repfolders = sorted(glob.glob(self.FEP + '/*/'))
        directories = [repfolders[0] + entry for entry in os.listdir(repfolders[0]) if os.path.isdir(repfolders[0] + entry)]
        for folder in directories:
            logfilesrep = sorted(glob.glob(folder + '/md*.log'))
            solute_rep_restraints = []
            shell_rep_restraints = []
            for file in logfilesrep:
                #print(file)
                with open(file) as logfile:
                    for line in logfile:
                        if line[:10] == "restraints":
                            solute_restraint = list(filter(None,line.split(" ")))[-1]
                            #print(solute_restraint)
                            shell_restraints = list(filter(None,line.split(" ")))[-2]
                            #print(shell_restraints)
                            solute_rep_restraints.append(float(solute_restraint))
                            shell_rep_restraints.append(float(shell_restraints))
                            break
                if solute_rep_restraints != []:
                    self.lib_for_df_solute[file.split("/")[-2]] = solute_rep_restraints
                    print_statement = ""
                if shell_rep_restraints != []:
                    self.lib_for_df_shell[file.split("/")[-2]] = shell_rep_restraints
                    print_statement = ""
                else:
                    self.lib_for_df_shell[file.split("/")[-2]] = 'nan'
                    self.lib_for_df_solute[file.split("/")[-2]] = 'nan'
                    print_statement = 'Could not retrieve restraints for:' + folder
            if print_statement == "":
                continue
            else:
                print(print_statement)

    def return_dG_BAR_method(self):
        repfolders = sorted(glob.glob(self.FEP + '/*/'))
        directories = [repfolders[0] + entry for entry in os.listdir(repfolders[0]) if os.path.isdir(repfolders[0] + entry)]
        for folder in directories:    
            qfepfilesrep = sorted(glob.glob(folder + '/qfep*.out'))
            for file in qfepfilesrep:
                # print (qfepfilesrep)
                with open(file, "r") as energies:
                    for line in energies:
                        pass
                    dGBAR = line
                    dGBAR = dGBAR.strip('\n').split(" ")[-1]
                    try:
                        dGBAR = float(dGBAR)
                        self.dG_BAR[file.split("/")[-2]] = dGBAR
                    except ValueError:
                        self.dG_BAR[file.split("/")[-2]] = "nan"

    def print_results(self):
        averages_restraints_solute = []
        averages_restraints_shell = []
        dG_BAR_values = []
        for rep_key in self.lib_for_df_solute:
            if self.lib_for_df_solute[rep_key] == [] or "nan" in self.lib_for_df_solute[rep_key] or self.dG_BAR[rep_key] == "nan":
                averages_restraints.append('nan')
                dG_BAR_values.append('nan')
            else:
                average_solute_restraint = statistics.mean(self.lib_for_df_solute[rep_key])
                averages_restraints_solute.append(str(round(average_solute_restraint,3)))
                average_shell_restraint = statistics.mean(self.lib_for_df_shell[rep_key])
                averages_restraints_shell.append(str(round(average_shell_restraint,3)))
                dG_BAR_values.append(str(self.dG_BAR[rep_key]))
        print("dG and restraints")
        print(', '.join(dG_BAR_values))
        print(', '.join(averages_restraints_solute))
        print(', '.join(averages_restraints_shell))

    # Function to calculate angles relative to xy-plane and z-axis
    def calculate_angles(self, line_vector):
        xy_angle_rad = np.arctan2(line_vector[1], line_vector[0])
        xy_angle_deg = np.degrees(xy_angle_rad)

        z_angle_rad = np.arccos(line_vector[2] / np.linalg.norm(line_vector))
        z_angle_deg = np.degrees(z_angle_rad)
        return xy_angle_deg, z_angle_deg

    def check_alpha_helix(self):
        residue_list = ["108", "109", "110", "111", "112", "113"]
        fixed_point_residue = "107"
        tilted_alfahelix_angles_xy = []
        tilted_alfahelix_angles_z = []
        # coordinates_alphahelix_Ca_1 = []
        # coordinates_fixed_point_1 = ""

        # with open("reference.pdb") as pdbfile:
        #     for line in pdbfile:
        #         line_in_list = list(filter(None,line.split(" ")))
        #         # print(line_in_list[4])
        #         try:
        #             if line_in_list[4] in residue_list:
        #                 # print(line)
        #                 if "CA" in line_in_list:
        #                     # print([float(line_in_list[5]),float(line_in_list[6]),float(line_in_list[7])])
        #                     coordinates_alphahelix_Ca_1.append([float(line_in_list[5]),float(line_in_list[6]),float(line_in_list[7])])
        #             if line_in_list[4] == fixed_point_residue and "CA" in line_in_list:
        #                 coordinates_fixed_point_1 = [float(line_in_list[5]),float(line_in_list[6]),float(line_in_list[7])]
        #         except IndexError:
        #             continue
        # # Calculate lines passing through centroids of both sets
        # centroid_1 = np.mean(coordinates_alphahelix_Ca_1, axis=0)
        # direction_1 = centroid_1 - coordinates_fixed_point_1
        # # Calculate angles for Alpha-Helix 1
        # xy_angle_1, z_angle_1 = run.calculate_angles(direction_1)
        # print(180+xy_angle_1)
        # print(180-z_angle_1)

        repfolders = sorted(glob.glob(self.FEP + '/*/'))
        directories = [repfolders[0] + entry for entry in os.listdir(repfolders[0]) if os.path.isdir(repfolders[0] + entry)]
        for folder in directories:    
            qfepfilesrep = sorted(glob.glob(folder + '/md*_2.pdb'))
            for md_file in qfepfilesrep:
                #print(md_file)
                coordinates_alphahelix_Ca_2 = []
                coordinates_fixed_point_2 = ""
                with open(md_file) as pdbfile:
                    for line in pdbfile:
                        line_in_list = list(filter(None,line.split(" ")))
                        # print(line_in_list[4])
                        try:
                            if line_in_list[4] in residue_list:
                                # print(line)
                                if "CA" in line_in_list:
                                    # print([float(line_in_list[5]),float(line_in_list[6]),float(line_in_list[7])])
                                    coordinates_alphahelix_Ca_2.append([float(line_in_list[5]),float(line_in_list[6]),float(line_in_list[7])])
                            if line_in_list[4] == fixed_point_residue and "CA" in line_in_list:
                                coordinates_fixed_point_2 = [float(line_in_list[5]),float(line_in_list[6]),float(line_in_list[7])]
                        except IndexError:
                            continue

                # Calculate lines passing through centroids of both sets
                centroid_2 = np.mean(coordinates_alphahelix_Ca_2, axis=0)
                direction_2 = centroid_2 - coordinates_fixed_point_2
                # Calculate angles for Alpha-Helix 1
                xy_angle_2, z_angle_2 = run.calculate_angles(direction_2)
                xy_angle = 180+xy_angle_2
                z_angle = 180-z_angle_2
                tilted_alfahelix_angles_xy.append(str(round(xy_angle,2)))
                tilted_alfahelix_angles_z.append(str(round(z_angle,2)))
                # print(180+xy_angle_2 180-z_angle_2)
        #self.md_file_length = 112
        xy_angle_per_replicate = list(run.chunks(tilted_alfahelix_angles_xy, self.md_file_length))
        xy_angle_tabular_form = list(zip(*xy_angle_per_replicate))
        print("helix xy plane tilt")
        for xy_angle_per_rep in xy_angle_tabular_form:
            print(', '.join(xy_angle_per_rep))
        
        z_angle_per_replicate = list(run.chunks(tilted_alfahelix_angles_z, self.md_file_length))
        z_angle_tabular_form = list(zip(*z_angle_per_replicate))
        print("helix z plane tilt")
        for z_angle_per_rep in z_angle_tabular_form:
            print(', '.join(z_angle_per_rep))

    def chunks(self, lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

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

    parser.add_argument('-align', '--align_pdbs',
                    dest = "align_pdbs",
                    required = False,
                    action='store_true',
                    help = "align trajectory pdb files")

    args = parser.parse_args()
    run = Run(FEP = args.FEP,
              INPUTS = args.INPUTS,
              FF = args.FF,
              convert_files = args.convert_files,
              align_pdbs = args.align_pdbs)

    run.make_pdb_files()
    run.align_to_reference_benzene()
    run.find_RMSD_ligand()
    run.check_alpha_helix()
    run.return_dG_BAR_method()
    run.return_restraints()
    run.print_results()
    run.save_aligned_trajectory_pdb_in_single_state()
