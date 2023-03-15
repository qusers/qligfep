import argparse
import re
import glob
import os
import shutil
import stat
import shlex
from subprocess import check_output
from rdkit import Chem 
from rdkit.Chem import rdmolfiles

import functions as f
import settings as s
import IO

class Run(object):
    """
    Create dual topology FEP files based on two ligands
    """
    def __init__(self, 
                 lig1, 
                 FF, 
                 system, 
                 cluster, 
                 sphereradius, 
                 cysbond, 
                 start, 
                 temperature, 
                 replicates,
                 sampling,
                 *args, 
                 **kwargs):
        
        self.lig = lig
        self.FF = FF
        self.system = system
        self.rootdir = os.getcwd()
        self.cluster = cluster
        self.sphereradius = sphereradius
        self.cysbond = cysbond
        self.start = start
        self.include = ['ATOM', 'HETATM']
        self.temperature = temperature
        self.replicates = replicates
        self.sampling = sampling
        self.final_mol_number = 1
        # needs to be user defined
        self.waters_to_perturb = 26
    

    def replace(self, string, replacements):
        pattern = re.compile(r'\b(' + '|'.join(replacements.keys()) + r')\b')
        replaced_string = pattern.sub(lambda x: replacements[x.group()], string)
        return replaced_string


    def makedir(self):
        directory1 = self.rootdir + '/FEP_' + self.lig
        if not os.path.exists(directory):
            os.makedirs(directory)
            
        if not os.path.exists(directory + '/inputfiles'):
            os.makedirs(directory + '/inputfiles')

        return directory

    def retreive_lig_coordinates():
        oxygencoords = []
        
        ligand = rdmolfiles.MolFromPDBFile(self.lig + '.pdb')

        for atom in .GetAtoms():
            atomindex = atom.GetIdx()
            if atomindex % 2:
                continue
            else:
                coordinates = ligand.GetConformer().GetAtomPosition(atomindex)
                oxygencoords.append((format(round(coordinates.x + 0.001,3), '.3f'),
                format(round(coordinates.y + 0.001,3), '.3f'),
                format(round(coordinates.z + 0.001,3), '.3f')))
        return oxygencoords


    def write_waters_to_pdb(coordinates_list, solvent_waters = False):

        with open('oxygens.pdb', "a+") as oxygenfile:
            
            for watnumber in range(water_amount):

                iteration_oxygen_coords = oxygencoords[watnumber]
                length_xcoord = len(str(iteration_oxygen_coords[0]))
                length_ycoord = len(str(iteration_oxygen_coords[1]))
                length_zcoord = len(str(iteration_oxygen_coords[2]))
                xspaces = 12 - length_xcoord
                yspaces = 8  - length_ycoord
                zspaces = 8  - length_zcoord

                new_molnumber = self.final_mol_number + 1
                molspaces = 6 - len(str(new_molnumber))

                oxygenfile.write('HETATM    1  O1  WAT')
                oxygenfile.write(molspaces*' ' + str(new_molnumber))
                oxygenfile.write(xspaces*' ' + iteration_oxygen_coords[0] +
                                yspaces*' ' + iteration_oxygen_coords[1]+
                                zspaces*' ' + iteration_oxygen_coords[2])
                oxygenfile.write('  0.00  0.00           O  \n')

                self.final_mol_number = new_molnumber

    # load qprepfile with waters solvated
    def get_solvent_oxygens():
        axescoords = []
        trigonalcoords = []
        diagonalcoords = []

        with open('complexnotexcluded.pdb') as complexx:
            for line in complexx:
                splitted = line.strip('\n').split(' ')
                line_in_list = list(filter(None, splitted))
                try:
                    # get waters only
                    if line_in_list[3] == 'HOH':
                        coordinates = line_in_list[5:]
                        if coordinates[0].strip('-') == coordinates[1].strip('-') \
                        or coordinates[1].strip('-') == coordinates[2].strip('-') \
                        or coordinates[0].strip('-') == coordinates[2].strip('-'):
                            if coordinates[0] == '0.000' \
                            or coordinates[1] == '0.000' \
                            or coordinates[2] == '0.000':
                                if coordinates[1] == '0.000' \
                                and coordinates[2] == '0.000':
                                    axescoords.append(coordinates)
                            if coordinates[0] == coordinates[1] \
                            and coordinates[2] == '0.000':
                                diagonalcoords.append(coordinates)
                            elif coordinates[0].strip('-') == coordinates[1].strip('-') == coordinates[2].strip('-'):
                                trigonalcoords.append(coordinates)
                                
                except IndexError:
                    continue

        # coordinates on axes
        highest_axis_axes = axescoords[-2][0]
        pos_x_axis_water = [highest_axis_axes, '0.000', '0.000']
        neg_x_axis_water = ['-' + highest_axis_axes, '0.000', '0.000']
        pos_y_axis_water = ['0.000', highest_axis_axes, '0.000']
        neg_y_axis_water = ['0.000', '-' + highest_axis_axes, '0.000']
        pos_z_axis_water = ['0.000','0.000', highest_axis_axes]
        neg_z_axis_water = ['0.000','0.000', '-' + highest_axis_axes]

        # coordinates on trigonals
        highest_trigonalcoords = trigonalcoords[-1][0]
        pos_pos_pos_trigonal = [highest_trigonalcoords, highest_trigonalcoords, highest_trigonalcoords]
        neg_neg_neg_trigonal = ['-' + highest_trigonalcoords, '-' + highest_trigonalcoords, '-' + highest_trigonalcoords]
        pos_neg_pos_trigonal = [highest_trigonalcoords, '-' + highest_trigonalcoords, highest_trigonalcoords]
        neg_pos_neg_trigonal = ['-' + highest_trigonalcoords, highest_trigonalcoords, '-' + highest_trigonalcoords]
        neg_neg_pos_trigonal = ['-' + highest_trigonalcoords, '-' + highest_trigonalcoords, highest_trigonalcoords]
        pos_pos_neg_trigonal = [highest_trigonalcoords, highest_trigonalcoords, '-' + highest_trigonalcoords]
        pos_neg_neg_trigonal = [highest_trigonalcoords, '-' + highest_trigonalcoords, '-' + highest_trigonalcoords]
        neg_pos_pos_trigonal = ['-' + highest_trigonalcoords, highest_trigonalcoords, highest_trigonalcoords]

        # coordinates on diagonal
        highest_diagonalcoords = diagonalcoords[-1][0] 
        pos_x_axis_pos_y_axis_water = [highest_diagonalcoords, highest_diagonalcoords, '0.000']
        neg_x_axis_neg_y_axis_water = [ '-' + highest_diagonalcoords, '-' + highest_diagonalcoords, '0.000']
        pos_x_pos_z_axis_water = [highest_diagonalcoords,'0.000', highest_diagonalcoords]
        neg_x_neg_z_axis_water = [ '-' + highest_diagonalcoords,'0.000', '-' + highest_diagonalcoords]
        pos_y_axis_pos_z_axis_water = ['0.000', highest_diagonalcoords, highest_diagonalcoords]
        neg_y_axis_neg_z_axis_water = ['0.000', '-' + highest_diagonalcoords, '-' + highest_diagonalcoords]
        pos_x_axis_neg_y_axis_water = [highest_diagonalcoords,  '-' + highest_diagonalcoords, '0.000']
        neg_x_axis_pos_y_axis_water = [ '-' + highest_diagonalcoords, highest_diagonalcoords, '0.000']
        pos_x_neg_z_axis_water = [highest_diagonalcoords,'0.000', '-' + highest_diagonalcoords]
        neg_x_pos_z_axis_water = [ '-' + highest_diagonalcoords,'0.000', highest_diagonalcoords]
        pos_y_axis_neg_z_axis_water = ['0.000', highest_diagonalcoords, '-' + highest_diagonalcoords]
        neg_y_axis_pos_z_axis_water = ['0.000',  '-' + highest_diagonalcoords, highest_diagonalcoords]

        all_solvent_waters = [pos_x_axis_water,
                      neg_x_axis_water,
                      pos_y_axis_water,
                      neg_y_axis_water,
                      pos_z_axis_water,
                      neg_z_axis_water,
                      pos_pos_pos_trigonal,
                      neg_neg_neg_trigonal,
                      pos_neg_pos_trigonal,
                      neg_pos_neg_trigonal,
                      neg_neg_pos_trigonal,
                      pos_pos_neg_trigonal,
                      pos_neg_neg_trigonal,
                      neg_pos_pos_trigonal,
                      pos_x_axis_pos_y_axis_water,
                      neg_x_axis_neg_y_axis_water,
                      pos_x_pos_z_axis_water,
                      neg_x_neg_z_axis_water,
                      pos_y_axis_pos_z_axis_water,
                      neg_y_axis_neg_z_axis_water,
                      pos_x_axis_neg_y_axis_water,
                      neg_x_axis_pos_y_axis_water,
                      pos_x_neg_z_axis_water,
                      neg_x_pos_z_axis_water,
                      pos_y_axis_neg_z_axis_water,
                      neg_y_axis_pos_z_axis_water]
                      
        return all_solvent_waters

    def remove_waters(solvent_waters):
        waterfile_lines = []

        with open('complexnotexcluded.pdb') as complexx:
            for line in complexx:
                splitted = line.strip('\n').split(' ')
                line_in_list = list(filter(None, splitted))
                try:
                    # get waters only
                    if line_in_list[3] == 'HOH':
                        coordinates = line_in_list[5:]
                        for solvent_water in solvent_waters:
                            if solvent_water[0]==coordinates[0] \
                            and solvent_water[1]==coordinates[1] \
                            and solvent_water[2]==coordinates[2] in line:
                                waternumber = line_in_list[4]
                        if waternumber == line_in_list[4]:
                            continue
                        else:
                            waterfile_lines.append(line)

                except IndexError:
                    continue

        # need to add sphere radius here
        with open('waterfile.pdb', 'w') as waterfile:
            waterfile.write('22 sphere\n')
            for line in waterfile_lines:
                waterfile.write(line)


    def prepare_complex():
        replacements = {}
        center = f.COG(self.lig + '.pdb')
        center = '{:} {:} {:}'.format(center[0], center[1], center[2])
        qprep_in = s.ROOT_DIR + '/INPUTS/qprep.inp'
        qprep_out = writedir + '/qprep.inp'
        replacements['FF_LIB'] = s.ROOT_DIR + '/FF/' + self.FF + '.lib'
        replacements['LIG1']   = self.lig + '.lib'
        replacements['LIG2']   = 'templates/template_water_lib.lib'
        replacements['LIGPRM'] = self.FF + '_' + self.lig + 'water_merged.prm'
        replacements['LIGPDB'] = self.lig + '.pdb'
        replacements['CENTER'] = center
        replacements['SPHERE'] = self.sphereradius
        if self.system =='vacuum':
            replacements['solvate'] = '!solvate'
        if self.system == 'water':
            replacements['SOLVENT']  = '1 HOH'         

        with open(qprep_in) as infile, open(qprep_out, 'w') as outfile:
            for line in infile:
                line = run.replace(line, replacements)
                if line == '!addbond at1 at2 y\n' and self.cysbond != None:
                    cysbond = self.cysbond.split(',')
                    for cys in cysbond:
                        at1 = cys.split('_')[0]
                        at2 = cys.split('_')[1]
                        outfile.write('addbond ' + at1 + ' ' + at2 + ' y \n' )
                    continue
                outfile.write(line)

    def qprep_initial(self, writedir):
        os.chdir(writedir)
        cluster_options = getattr(s, self.cluster)
        qprep = cluster_options['QPREP']
        options = ' < qprep.inp > qprep.out'
        # Somehow Q is very annoying with this < > input style so had to implement
        # another function that just calls os.system instead of using the preferred
        # subprocess module....
        IO.run_command(qprep, options, string = True)
        os.chdir('../../')


    def change_prm(self, replacements, writedir):

        # find ligand file and (stored in list. Thats why the [0] is there)
        lig_prm_file = glob.glob(self.lig + '.prm')[0]
        # for now only amber is possible
        prm_file = s.FF_DIR + '/' + self.FF + '.prm'
        prm_merged = {'vdw':[],
                       'bonds':[],
                       'angle':[],
                       'torsion':[],
                       'improper':[]
                       }
        
        # need to have a directory with all templates such as template_water.prm
        for file in ['template_water.prm', lig_prm_file]:
            with open(file) as infile:
                block = 0
                for line in infile:
                    if line == '[atom_types]\n':
                        block = 1
                        continue
                    elif line == '[bonds]\n':
                        block = 2
                        continue
                    elif line == '[angles]\n':
                        block = 3
                        continue
                    elif line == '[torsions]\n':
                        block = 4
                        continue
                    if line == '[impropers]\n':
                        block = 5
                        continue

                    if block == 1:
                        prm_merged['vdw'].append(line)

                    elif block == 2:
                        prm_merged['bonds'].append(line)                

                    elif block == 3:
                        prm_merged['angle'].append(line)                

                    elif block == 4:
                        prm_merged['torsion'].append(line)                

                    elif block == 5:
                        prm_merged['improper'].append(line)
                        
        with open(prm_file) as infile, open(writedir + 
                                            '/' + 
                                            self.FF + 
                                            '_' + 
                                            self.lig +
                                            'water_merged.prm', 'w') as outfile:
            for line in infile:
                block = 0
                outfile.write(line)
                if len(line) > 1:
                    if line == "! Ligand vdW parameters\n":
                        block = 1
                    if line == "! Ligand bond parameters\n":
                        block = 2     
                    if line == "! Ligand angle parameters\n":
                        block = 3
                    if line == "! Ligand torsion parameters\n":
                        block = 4
                    if line == "! Ligand improper parameters\n":
                        block = 5
                # Read the parameters in from file and store them
                if block == 1: 
                    for line in prm_merged['vdw']:
                        outfile.write(line)

                if block == 2:
                    for line in prm_merged['bonds']:
                        outfile.write(line)

                if block == 3:
                    for line in prm_merged['angle']:
                        outfile.write(line)

                if block == 4:
                    for line in prm_merged['torsion']:
                        outfile.write(line)

                if block == 5:
                    for line in prm_merged['improper']:
                        outfile.write(line)                  
        
        #AND return the vdW list for the FEP file
        FEP_vdw = []
        for line in prm_merged['vdw']:
            if len(line) > 1 and line[0] != '!' and line[0:1]:
                line = line.split()
                line2 = "{:10}{:10}{:10}{:10}{:10}{:10}{:10}{:10}".format(line[0],
                                                                      line[1],
                                                                      line[3],
                                                                      str(0),
                                                                      str(0),
                                                                      line[4],
                                                                      line[5],
                                                                      line[6]
                                                                      )
                FEP_vdw.append(line2)
        FEP_vdw.append('DUM       0.0000    0.0000    0         0         0.0000    0.0000    1.0080')
        return FEP_vdw

    def read_files(self):
        charges = []
        atomtypes = []
        merged_molsize = 0
        
        with open(self.lig1 + '.lib') as infile:
            block = 0
            for line in infile:
                line = line.split()
                if len(line) > 0:
                    if line[0] == '[atoms]':
                        block = 1
                        continue
                    if line[0] == '[bonds]':
                        block = 2
                        
                if block == 1 and len(line) > 0:
                    #construct for FEP file
                    merged_molsize = molsize + 1
                    charges.append([merged_molsize, line[3], '0.000'])
                    atomtypes.append([merged_molsize, line[2], 'DUM'])
                
                if block == 2:
                    break
                    
            ligmolsize = merged_molsize

        # maybe add template forcefield water charges?
        for water in range(self.waters_to_perturb):
            merged_molsize = molsize + 1
            charges.append([merged_molsize, '-0.834', '0.000'])
            charges.append([merged_molsize, '0.417', '0.000'])
            charges.append([merged_molsize, '0.417', '0.000'])

        for water in range(self.waters_to_perturb):
            merged_molsize = molsize + 1
            charges.append([merged_molsize,'0.000', '-0.834'])
            charges.append([merged_molsize, '0.000', '0.417'])
            charges.append([merged_molsize, '0.000', '0.417'])

        return ligmolsize, charges, atomtypes, merged_molsize

    def write_FEP1_file(self, change_charges, change_vdw, FEP_vdw):      
        
        with open(writedir + '/FEP1.fep', 'w') as outfile:
            total_atoms = len(change_charges)
            outfile.write('!info: ' + self.lig1 + '_water_discharge -->\n')
            outfile.write('[FEP]\n')
            outfile.write('states 2\n')
            outfile.write('softcore_use_max_potential on\n\n')

            # defining the atom order taken user given offset into account
            outfile.write('[atoms]\n')
            for i in range(1, change_charges + 1):
                outfile.write("{:5}{:5}\n".format(str(i),
                                                    str(i)))
            outfile.write('\n\n')
            
            # changing charges
            outfile.write('[change_charges]\n')

            for line in change_charges:
                outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                        line[1],
                                                        line[2]))
            outfile.write('\n\n')
            
            #add the Q atomtypes
            outfile.write('[atom_types]\n')
            for line in FEP_vdw:
                outfile.write(line + '\n')
            outfile.write('\n\n')
            
            outfile.write('[softcore]\n')
            # ADD softcore
            for i in range(1, change_charges + 1):
                outfile.write("{:<5}{:>10}{:>10}\n".format(str(i),'0', '0'))
            
            outfile.write('\n\n')
            
            # changing atom types
            ### for FEP 3 keep original ###
            outfile.write('[change_atoms]\n')
            for line in change_vdw:
                outfile.write('{:<5}{:>10}{:>10}\n'.format(line[0],
                                                        line[1],
                                                        line[1]))











# if __name__ == "__main__":
#     args = parseargs()
#     run = Run(lig1 = args.lig1,
#               FF= args.FF,
#               system = args.system,
#               cluster = args.cluster,
#               sphereradius = args.sphereradius,
#               cysbond = args.cysbond,
#               start = args.start,
#               temperature = args.temperature,
#               replicates = args.replicates,
#               sampling = args.sampling
#              )

#     writedir = run.makedir()
#     inputdir = writedir + '/inputfiles'
#     a = run.read_files()
#     changes_for_libfiles = a[0][1]
#     changes_for_prmfiles = a[0][1]
#     change_charges       = a[1][0]
#     change_vdw           = a[1][1]
#     changes_for_pdbfiles = a[0][0]
#     lig_size1, lig_size2 = a[2][0], a[2][1]

#     # Write the merged files
#     run.change_lib(changes_for_libfiles, inputdir)
#     FEP_vdw = run.change_prm(changes_for_prmfiles, inputdir)
#     run.write_FEP_file(change_charges, change_vdw, FEP_vdw, inputdir, lig_size1, lig_size2)
#     run.merge_pdbs(inputdir)
#     if args.system == 'protein':
#         run.write_water_pdb(inputdir)
#     lambdas = run.get_lambdas(args.windows, args.sampling)
#     overlapping_atoms = run.overlapping_atoms(writedir)
    
#     # Handling the correct offset here
#     if args.start == '0.5':
#         file_list = run.write_MD_05(lambdas, inputdir, lig_size1, lig_size2)
#         run.write_runfile(inputdir, file_list)    
        
#     if args.start == '1':
#         file_list = run.write_MD_1(lambdas, inputdir, lig_size1, lig_size2, overlapping_atoms)
#         run.write_runfile(inputdir, file_list)    
    
#     run.write_submitfile(writedir)
#     run.write_qfep(inputdir, args.windows, lambdas)
#     run.write_qprep(inputdir)
#     run.qprep(inputdir)