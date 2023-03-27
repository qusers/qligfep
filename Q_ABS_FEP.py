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
                 lig, 
                 FF,  
                 cluster, 
                 sphereradius, 
                 cysbond, 
                 temperature, 
                 replicates,
                 sampling1,
                 sampling2,
                 sampling3,
                 sampling4,
                 sampling5,
                 *args, 
                 **kwargs):
        
        self.lig = lig
        self.FF = FF
        self.rootdir = os.getcwd()
        self.cluster = cluster
        self.sphereradius = sphereradius
        self.cysbond = cysbond
        self.include = ['ATOM', 'HETATM']
        self.temperature = temperature
        self.replicates = replicates
        self.sampling1 = sampling1
        self.sampling2 = sampling2
        self.sampling3 = sampling3
        self.sampling4 = sampling4
        self.sampling5 = sampling5
        # needs to be user defined
        self.waters_to_perturb = 3
    

    def replace(self, string, replacements):
        pattern = re.compile(r'\b(' + '|'.join(replacements.keys()) + r')\b')
        replaced_string = pattern.sub(lambda x: replacements[x.group()], string)
        return replaced_string


    def makedir(self):
        directory = self.rootdir + '/FEP_' + self.lig

        if not os.path.exists(directory):
            os.makedirs(directory)
            
        if not os.path.exists(directory + '/inputfiles1'):
            os.makedirs(directory + '/inputfiles1')
        
        if not os.path.exists(directory + '/inputfiles2'):
            os.makedirs(directory + '/inputfiles2')

        if not os.path.exists(directory + '/inputfiles3'):
            os.makedirs(directory + '/inputfiles3')

        if not os.path.exists(directory + '/inputfiles4'):
            os.makedirs(directory + '/inputfiles4')

        if not os.path.exists(directory + '/inputfiles5'):
            os.makedirs(directory + '/inputfiles5')
        
        if not os.path.exists(directory + '/prep_input'):
            os.makedirs(directory + '/prep_input')

        return directory

    def retreive_lig_coordinates(self, prepdir):
        oxygencoords = []
        atomnumbers = []
        
        ligand = rdmolfiles.MolFromPDBFile(prepdir + self.lig + '.pdb')

        for atom in ligand.GetAtoms():
            atomindex = atom.GetIdx()
            
            if atomindex % 2:
                continue
            else:
                atomnumbers.append(atomindex+1)
                coordinates = ligand.GetConformer().GetAtomPosition(atomindex)
                oxygencoords.append([format(round(coordinates.x + 0.001,3), '.3f'),
                format(round(coordinates.y + 0.001,3), '.3f'),
                format(round(coordinates.z + 0.001,3), '.3f')])
        return oxygencoords, atomnumbers


    def write_waters_to_pdb(self, coordinates_list, input1folder):

        with open(input1folder + self.lig + '.pdb', "a+") as oxygenfile:
            oxygenfile.write('\n')
            for watnumber in range(len(coordinates_list)):

                iteration_oxygen_coords = coordinates_list[watnumber]
                length_xcoord = len(str(iteration_oxygen_coords[0]))
                length_ycoord = len(str(iteration_oxygen_coords[1]))
                length_zcoord = len(str(iteration_oxygen_coords[2]))
                xspaces = 12 - length_xcoord
                yspaces = 8  - length_ycoord
                zspaces = 8  - length_zcoord

                molnumber = watnumber + 2
                molspaces = 6 - len(str(molnumber))

                oxygenfile.write('HETATM    1  O1  WAT')
                oxygenfile.write(molspaces*' ' + str(molnumber))
                oxygenfile.write(xspaces*' ' + iteration_oxygen_coords[0] +
                                yspaces*' ' + iteration_oxygen_coords[1]+
                                zspaces*' ' + iteration_oxygen_coords[2])
                oxygenfile.write('  0.00  0.00           O  \n')


    # load qprepfile with waters solvated
    def get_solvent_oxygens(self, prepdir):
        axescoords = []
        trigonalcoords = []
        diagonalcoords = []

        with open(prepdir + 'prep_complex.pdb') as complexx:
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
                      
        return all_solvent_waters[:self.waters_to_perturb]

    def remove_waters(self, solvent_waters, prepdir, input1dir):
        waterfile_lines = []
        waternumber = ''

        with open(prepdir + 'prep_complex.pdb') as complexx:
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

        # need to add path?
        with open(input1dir + 'waterfile.pdb', 'w') as waterfile:
            waterfile.write(self.sphereradius + ' sphere\n')
            for line in waterfile_lines:
                waterfile.write(line)


    def prepare_complex(self, system, write_dir):
        replacements = {}
        center = f.COG(self.lig + '.lib')
        center = '{:} {:} {:}'.format(center[0], center[1], center[2])
        qprep_in = s.ROOT_DIR + '/INPUTS/qprep.inp'
        qprep_out = write_dir + '/qprep.inp'
        replacements['FF_LIB'] = s.ROOT_DIR + '/FF/' + self.FF + '.lib'
        replacements['LIG1']   = self.lig + '.lib'
        replacements['LIG2']   = s.ROOT_DIR + '/INPUTS/template_water_lib.lib'
        replacements['LIGPRM'] = self.FF + '_' + self.lig + '_water_merged.prm'
        replacements['LIGPDB'] = self.lig + '.pdb'
        replacements['CENTER'] = center
        replacements['SPHERE'] = self.sphereradius
        if system =='vacuum':
            replacements['solvate'] = '!solvate'

        if system == 'prepare':
            replacements['SOLVENT']  = '1 HOH'
            replacements['complexnotexcluded.pdb'] = 'prep_complex.pdb'
        if system == 'water':
            replacements['SOLVENT']  = '4 waterfile.pdb' 

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


    def qprep(self, writedir):
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
        water_prm_file = s.ROOT_DIR + '/INPUTS/template_water.prm'
        prm_merged = {'vdw':[],
                       'bonds':[],
                       'angle':[],
                       'torsion':[],
                       'improper':[]
                       }

        # need to have a directory with all templates such as template_water.prm
        for file in [water_prm_file, lig_prm_file]:
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
                                            '/prep_input/' + 
                                            self.FF + 
                                            '_' + 
                                            self.lig +
                                            '_water_merged.prm', 'w') as outfile:
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
        shutil.copy(writedir + '/prep_input/' + self.FF + '_' + self.lig + '_water_merged.prm', writedir + '/inputfiles1/' + self.FF + '_' + self.lig + '_water_merged.prm')
        shutil.copy(writedir + '/prep_input/' + self.FF + '_' + self.lig + '_water_merged.prm', writedir + '/inputfiles5/' + self.FF + '_' + self.lig + '_water_merged.prm')
        return FEP_vdw

    def read_files(self):
        charges = []
        atomtypes = []
        merged_molsize = 0

        with open(self.rootdir + '/' + self.lig + '.lib') as infile:
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
                    merged_molsize += 1
                    charges.append([merged_molsize, line[3]])
                    atomtypes.append([merged_molsize, line[2]])

                if block == 2:
                    break

            ligmolsize = merged_molsize

        ### Path needs to be adjusted ###
        ### use s.ROOT_DIR, see how its done in QligFEP ###
        for water in range(self.waters_to_perturb * 2):
            with open(s.ROOT_DIR + '/INPUTS/template_water_lib.lib') as waterfile:
                block = 0
                for line in waterfile:
                    line = line.split()
                    if len(line) > 0:
                        if line[0] == '[atoms]':
                            block = 1
                            continue
                        if line[0] == '[bonds]':
                            block = 2

                    if block == 1 and len(line) > 0:
                        #construct for FEP file
                        merged_molsize += 1
                        charges.append([merged_molsize, line[3]])
                        atomtypes.append([merged_molsize, line[2]])

                    if block == 2:
                        break

        return ligmolsize, charges, atomtypes, merged_molsize

    def write_FEP1_file(self, change_charges, change_vdw, FEP_vdw, writedir, ligsize):
        all_atoms_in_FEP_file = []

        with open(writedir + '/FEP1.fep', 'w') as outfile:
            total_atoms = len(change_charges)
            outfile.write('!info: ' + self.lig + '_water_discharge -->\n')
            outfile.write('[FEP]\n')
            outfile.write('states 2\n')
            outfile.write('softcore_use_max_potential on\n\n')

            # defining the atom order taken user given offset into account
            outfile.write('[atoms]\n')
            for i in range(1, len(change_charges) + 1):
                outfile.write("{:5}{:5}\n".format(str(i),
                                                    str(i)))
                all_atoms_in_FEP_file.append(i)
            outfile.write('\n\n')

            # changing charges
            outfile.write('[change_charges]\n')

            for atomnumber, line in enumerate(change_charges, start =1):
                if atomnumber <= ligsize:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            line[1],
                                                            '0.000'))

                # waters from solvent
                elif atomnumber <= ligsize + self.waters_to_perturb * 3:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            '0.000',
                                                            '0.000'))

                # solute waters
                else:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            line[1],
                                                            '0.000'))

            outfile.write('\n\n')

            #add the Q atomtypes
            outfile.write('[atom_types]\n')
            for line in FEP_vdw:
                outfile.write(line + '\n')
            outfile.write('\n\n')

            outfile.write('[softcore]\n')

            ### same as above can be done for other FEP files
            # ADD softcore
            for i in range(1, len(change_charges) + 1):
                outfile.write("{:<5}{:>10}{:>10}\n".format(str(i),'0', '0'))

            outfile.write('\n\n')

            # changing atom types
            ### for FEP 3 keep original ###
            outfile.write('[change_atoms]\n')
            for atomnumber, line in enumerate(change_vdw, start =1):
                if atomnumber <= ligsize:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            line[1],
                                                            line[1]))

                # waters from solvent
                elif atomnumber <= ligsize + self.waters_to_perturb * 3:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            'DUM',
                                                            'DUM'))

                # solute waters
                else:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            line[1],
                                                            line[1]))
        return all_atoms_in_FEP_file

    def write_FEP2_file(self, change_charges, change_vdw, FEP_vdw, writedir, ligsize):

        with open(writedir + '/FEP2.fep', 'w') as outfile:

            # change total atoms for len(change_charges) in range for loop
            total_atoms = len(change_charges)
            outfile.write('!info: ' + self.lig + '_add_softcore -->\n')
            outfile.write('[FEP]\n')
            outfile.write('states 2\n')
            outfile.write('softcore_use_max_potential on\n\n')

            # defining the atom order taken user given offset into account
            outfile.write('[atoms]\n')
            for i in range(1, len(change_charges) + 1):
                outfile.write("{:5}{:5}\n".format(str(i),
                                                    str(i)))
            outfile.write('\n\n')

            # changing charges
            outfile.write('[change_charges]\n')

            for i in range(1, len(change_charges) + 1):
                outfile.write("{:<5}{:>10}{:>10}\n".format(str(i),
                                                        '0.000', 
                                                        '0.000'))

            outfile.write('\n\n')

            #add the Q atomtypes
            outfile.write('[atom_types]\n')
            for line in FEP_vdw:
                outfile.write(line + '\n')
            outfile.write('\n\n')

            outfile.write('[softcore]\n')

            ### same as above can be done for other FEP files
            # ADD softcore
            for i in range(1, len(change_charges) + 1):
                outfile.write("{:<5}{:>10}{:>10}\n".format(str(i),'0', '20'))

            outfile.write('\n\n')

            # changing atom types
            ### for FEP 3 keep original ###
            outfile.write('[change_atoms]\n')

            for atomnumber, line in enumerate(change_vdw, start =1):
                if atomnumber <= ligsize:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            line[1],
                                                            line[1]))

                # waters from solvent
                elif atomnumber <= ligsize + self.waters_to_perturb * 3:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            'DUM',
                                                            'DUM'))

                # solute waters
                else:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            line[1],
                                                            line[1]))


    def write_FEP3_file(self, change_charges, change_vdw, FEP_vdw, writedir, ligsize):

        with open(writedir + '/FEP3.fep', 'w') as outfile:
            total_atoms = len(change_charges)
            outfile.write('!info: ' + self.lig + '_remove_vdW_and_softcore -->\n')
            outfile.write('[FEP]\n')
            outfile.write('states 2\n')
            outfile.write('softcore_use_max_potential on\n\n')

            # defining the atom order taken user given offset into account
            outfile.write('[atoms]\n')
            for i in range(1, len(change_charges) + 1):
                outfile.write("{:5}{:5}\n".format(str(i),
                                                    str(i)))
            outfile.write('\n\n')

            # changing charges
            outfile.write('[change_charges]\n')

            for i in range(1, len(change_charges) + 1):
                outfile.write("{:<5}{:>10}{:>10}\n".format(str(i),
                                                        '0.000', 
                                                        '0.000'))

            outfile.write('\n\n')

            #add the Q atomtypes
            outfile.write('[atom_types]\n')
            for line in FEP_vdw:
                outfile.write(line + '\n')
            outfile.write('\n\n')

            outfile.write('[softcore]\n')


            ### remove softcore from solvent waters
            # ADD softcore
            for atomnumber, line in enumerate(change_charges, start =1):
                if atomnumber <= ligsize:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            '20',
                                                            '20'))

                # waters from solvent
                elif atomnumber <= ligsize + self.waters_to_perturb * 3:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            '20',
                                                            '20'))

                # solute waters
                else:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            '20',
                                                            '0'))


            outfile.write('\n\n')

            # changing atom types
            ### for FEP 3 keep original ###
            outfile.write('[change_atoms]\n')

            for atomnumber, line in enumerate(change_vdw, start =1):
                if atomnumber <= ligsize:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            line[1],
                                                            line[1]))

                # waters from solvent
                elif atomnumber <= ligsize + self.waters_to_perturb * 3:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            'DUM',
                                                            'DUM'))

                # solute waters
                else:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            line[1],
                                                            'DUM'))

    def write_FEP4_file(self, change_charges, change_vdw, FEP_vdw, writedir, ligsize):

        with open(writedir + '/FEP4.fep', 'w') as outfile:
            total_atoms = len(change_charges)
            outfile.write('!info: ' + self.lig + '_water_RBFE -->\n')
            outfile.write('[FEP]\n')
            outfile.write('states 2\n')
            outfile.write('softcore_use_max_potential on\n\n')

            # defining the atom order taken user given offset into account
            outfile.write('[atoms]\n')
            for i in range(1, len(change_charges) + 1):
                outfile.write("{:5}{:5}\n".format(str(i),
                                                    str(i)))
            outfile.write('\n\n')

            # changing charges
            outfile.write('[change_charges]\n')

            for atomnumber, line in enumerate(change_charges, start =1):
                if atomnumber <= ligsize:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            '0.000',
                                                            '0.000'))

                # waters from solvent
                elif atomnumber <= ligsize + self.waters_to_perturb * 3:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            '0.000',
                                                            line[1]))

                # solute waters
                else:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            '0.000',
                                                            '0.000'))

            outfile.write('\n\n')

            #add the Q atomtypes
            outfile.write('[atom_types]\n')
            for line in FEP_vdw:
                outfile.write(line + '\n')
            outfile.write('\n\n')

            outfile.write('[softcore]\n')


            ### remove softcore from solvate waters and ligand
            # ADD softcore
            for atomnumber, line in enumerate(change_charges, start =1):
                if atomnumber <= ligsize:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            '20',
                                                            '0'))

                # waters from solvent
                elif atomnumber <= ligsize + self.waters_to_perturb * 3:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            '20',
                                                            '0'))

                # solute waters
                else:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            '0',
                                                            '0'))


            outfile.write('\n\n')

            # changing atom types
            ### for FEP 3 keep original ###
            outfile.write('[change_atoms]\n')

            for atomnumber, line in enumerate(change_vdw, start =1):
                if atomnumber <= ligsize:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            line[1],
                                                            'DUM'))

                # waters from solvent
                elif atomnumber <= ligsize + self.waters_to_perturb * 3:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            'DUM',
                                                            line[1]))

                # solute waters
                else:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            'DUM',
                                                            'DUM'))

    # need other qprep inputfile
    def write_FEP5_file(self, change_charges, change_vdw, FEP_vdw, writedir, ligsize):

        with open(writedir + '/FEP5.fep', 'w') as outfile:
            total_atoms = len(change_charges)
            outfile.write('!info: ' + self.lig + '_water_RBFE -->\n')
            outfile.write('[FEP]\n')
            outfile.write('states 2\n')
            outfile.write('softcore_use_max_potential on\n\n')

            # defining the atom order taken user given offset into account
            outfile.write('[atoms]\n')
            for i in range(1, len(change_charges) + 1):
                if i <= ligsize:
                    outfile.write("{:5}{:5}\n".format(str(i),
                                                        str(i)))
            outfile.write('\n\n')

            # changing charges
            outfile.write('[change_charges]\n')

            for atomnumber, line in enumerate(change_charges, start =1):
                if atomnumber <= ligsize:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            line[1],
                                                            '0.000'))

            outfile.write('\n\n')

            #add the Q atomtypes
            outfile.write('[atom_types]\n')
            for line in FEP_vdw:
                outfile.write(line + '\n')
            outfile.write('\n\n')

            outfile.write('[softcore]\n')


            ### remove softcore from solvate waters and ligand
            # ADD softcore
            for i in range(1, len(change_charges) + 1):
                if i <= ligsize:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(str(i),'0', '0'))

            outfile.write('\n\n')

            # changing atom types
            ### for FEP 3 keep original ###
            outfile.write('[change_atoms]\n')

            for atomnumber, line in enumerate(change_vdw, start =1):
                if atomnumber <= ligsize:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            line[1],
                                                            'DUM'))
    
    def copyfiles(self):
        shutil.copy(self.lig + '.lib', prepinputdir + self.lig + '.lib')
        shutil.copy(self.lig + '.pdb', prepinputdir + self.lig + '.pdb')
        shutil.copy(self.lig + '.lib', inputdir1 + self.lig + '.lib')
        shutil.copy(self.lig + '.pdb', inputdir1 + self.lig + '.pdb')
        shutil.copy(self.lig + '.lib', inputdir5 + self.lig + '.lib')
        shutil.copy(self.lig + '.pdb', inputdir5 + self.lig + '.pdb')

    def copyfiles2(self):
        shutil.copy(inputdir1 + 'dualtop.top', inputdir2 + 'dualtop.top')
        shutil.copy(inputdir1 + 'dualtop.top', inputdir3 + 'dualtop.top')
        shutil.copy(inputdir1 + 'dualtop.top', inputdir4 + 'dualtop.top')


    def write_MD(self, lambdas, writedir, total_atoms, atomnumbers_lig, ligmolsize, step):
        totallambda = len(lambdas)
        replacements = {}

        lig_water_oxygens = total_atoms[ligmolsize: ligmolsize + self.waters_to_perturb*3:][::3]
        overlapping_atoms = list(zip(atomnumbers_lig, lig_water_oxygens))
        waters_begin = total_atoms[ligmolsize:][::3]
        waters_end = total_atoms[ligmolsize:][::-3][::-1]
        water_indexes_seq_restr = list(zip(waters_begin, waters_end))
        lig_indexes_seq_restr = [('1', ligmolsize)]

        if step == 'step1' or step == 'step2' or step == 'step3':
            only_solvent_waters = water_indexes_seq_restr[int(len(water_indexes_seq_restr)/2):]
            all_sequence_restraints = lig_indexes_seq_restr + only_solvent_waters

        if step == 'step4':
            # get the last index of solvate waters
            all_sequence_restraints = [('1', water_indexes_seq_restr[:int(len(water_indexes_seq_restr)/2)][-1][-1])]

        if step == 'step5':
            all_sequence_restraints = lig_indexes_seq_restr

        replacements['SPHERE']          =   self.sphereradius    
        replacements['EQ_LAMBDA']       =   '1.000 0.000'

        if step == 'step1' or step == 'step5':
            for eq_file_in in sorted(glob.glob(s.ROOT_DIR + '/INPUTS/eq*.inp')):
                eq_file = eq_file_in.split('/')[-1:][0]
                eq_file_out = writedir + '/' + eq_file

                with open(eq_file_in) as infile, open(eq_file_out, 'w') as outfile:
                    for line in infile:
                        line = run.replace(line, replacements)
                        # print(line)
                        if 'WATER_RESTRAINT' in line:
                            for seq in all_sequence_restraints:
                                outfile.write('{:<7}{:<8} 1.0 0 1\n'.format(seq[0], seq[1]))
                        elif 'ATOM_START_LIG1' in line:
                            force = line[25:]
                            for seq in all_sequence_restraints:
                                outfile.write('{:<7}{:<8}{forceinput}'.format(seq[0], seq[1], forceinput = force))
                        elif 'distance_restraints' in line:
                            outfile.write(line)
                            for indexxx in overlapping_atoms:
                                outfile.write('{:d} {:d} 0.0 0.2 0.5 0\n'.format(indexxx[0], indexxx[1]))
                        else:
                            outfile.write(line)

        file_in = s.INPUT_DIR + '/md_1000_0000.inp'
        file_out = writedir + '/md_1000_0000.inp' 
        with open(file_in) as infile, open(file_out, 'w') as outfile:
            for line in infile:
                line = run.replace(line, replacements)
                if line == '[distance_restraints]\n':
                    outfile.write(line)
                    for line in overlapping_atoms:
                        outfile.write('{:d} {:d} 0.0 0.2 0.5 0\n'.format(line[0], line[1]))
                elif 'WATER_RESTRAINT' in line:
                    for seq in all_sequence_restraints:
                        outfile.write('{:<7}{:<8} 1.0 0 1\n'.format(seq[0], seq[1]))
                else:
                    outfile.write(line)

        filenr = 0

        for l in lambdas:
            if l == '1.000':
                filename_N = 'md_1000_0000'
                continue
            else:
                step_n = totallambda - filenr - 2

                lambda1 = l
                lambda2 = lambdas[step_n]
                filename = 'md_' + lambda1.replace('.', '') + '_' + lambda2.replace('.', '')
                replacements['FLOAT_LAMBDA1']   =   lambda1 
                replacements['FLOAT_LAMBDA2']   =   lambda2
                replacements['FILE']          =   filename
                replacements['FILE_N'] = filename_N

                # Move to functio
                pattern = re.compile(r'\b(' + '|'.join(replacements.keys()) + r')\b')
                file_in = s.INPUT_DIR + '/md_XXXX_XXXX.inp'
                file_out = writedir + '/' + filename + '.inp'

                with open(file_in) as infile, open(file_out, 'w') as outfile:
                    for line in infile:
                        line = pattern.sub(lambda x: replacements[x.group()], line)
                        if 'WATER_RESTRAINT' in line:
                            for seq in all_sequence_restraints:
                                outfile.write('{:<7}{:<8} 1.0 0 1\n'.format(seq[0], seq[1]))
                        elif 'ATOM_START_LIG1' in line:
                            force = line[25:]
                            for seq in all_sequence_restraints:
                                outfile.write('{:<7}{:<8}{forceinput}'.format(seq[0], seq[1], forceinput = force))
                        elif 'distance_restraints' in line:
                            outfile.write(line)
                            for indexxx in overlapping_atoms:
                                outfile.write('{:d} {:d} 0.0 0.2 0.5 0\n'.format(indexxx[0], indexxx[1]))
                        else:
                            outfile.write(line)
                filename_N = filename
                filenr += 1

    def get_lambdas(self, windows, sampling):
        # Constructing the lambda partition scheme
        windows = int(windows)
        step = int(windows/2)
        lambdas = []
        lmbda_1 = []
        lmbda_2 = []
        k_dic = {'sigmoidal':-1.1, 
                 'linear':1000,
                 'exponential':-1.1,
                 'reverse_exponential':1.1
                }
        k = k_dic[sampling]

        if sampling == 'sigmoidal': 
            for i in range(0, step + 1):
                lmbda1 = '{:.3f}'.format(0.5 * (f.sigmoid(float(i)/float(step), k) + 1))
                lmbda2 = '{:.3f}'.format(0.5 * (-f.sigmoid(float(i)/float(step), k) + 1))
                lmbda_1.append(lmbda1)
                lmbda_2.append(lmbda2)

            lmbda_2 = lmbda_2[1:]

            for i in reversed(lmbda_2):
                lambdas.append(i)

            for i in lmbda_1:
                lambdas.append(i)

        else:
            for i in range(0, windows + 1):
                lmbda = '{:.3f}'.format(f.sigmoid(float(i)/float(windows), k))
                lambdas.append(lmbda)
        
        lambdas = lambdas[::-1]
        return lambdas

    def write_submitfile(self, writedir):
        replacements = {}
        replacements['TEMP_VAR']    = self.temperature
        replacements['RUN_VAR']     = self.replicates
        replacements['RUNFILE']     = 'run' + self.cluster + '.sh'
        submit_in = s.ROOT_DIR + '/INPUTS/FEP_submit.sh'
        submit_out = writedir + ('/FEP_submit.sh')
        with open(submit_in) as infile, open (submit_out, 'w') as outfile:
            for line in infile:
                line = run.replace(line, replacements)
                if 'inputfiles=$workdir/inputfiles' in line:
                    continue
                elif 'submitfile=$inputfiles/' in line:
                    submitline = line
                    outfile.write('inputefiles_fep1=$workdir/inputfiles1\n')
                    outfile.write('inputefiles_fep2=$workdir/inputfiles2\n')
                    outfile.write('inputefiles_fep3=$workdir/inputfiles3\n')
                    outfile.write('inputefiles_fep4=$workdir/inputfiles4\n')
                    outfile.write('inputefiles_fep5=$workdir/inputfiles5\n')
                    outfile.write('for inputfiles in "$inputefiles_fep1" "$inputefiles_fep2" "$inputefiles_fep3" "$inputefiles_fep4" "$inputefiles_fep5"; do\n')
                    outfile.write('submitfile=$inputfiles/runALICE.sh')
                elif 'sed -i s#inputfiles=.*#inputfiles="$inputfiles"#g $submitfile' in line:
                    outfile.write(line)
                    outfile.write('done\n\n')
                    outfile.write('for inputfiles in "$inputefiles_fep1" "$inputefiles_fep2" "$inputefiles_fep3" "$inputefiles_fep4" "$inputefiles_fep5"; do\n')
                    outfile.write(submitline)
                else:
                    outfile.write(line)
            outfile.write('wait\n')
            outfile.write('done')

        try:
            st = os.stat(submit_out)
            os.chmod(submit_out, st.st_mode | stat.S_IEXEC)
        
        except:
            print("WARNING: Could not change permission for " + submit_out)


    def write_runfile(self, writedir, step):
        ntasks = getattr(s, self.cluster)['NTASKS']
        src = s.INPUT_DIR + '/run_abs_solv.sh'
        tgt = writedir + '/run' + self.cluster + '.sh'
        EQ_files = sorted(glob.glob(writedir + '/eq*.inp'))
        MD_files = reversed(sorted(glob.glob(writedir + '/md*.inp')))
        replacements = getattr(s, self.cluster)

        if step == 'step1':
            replacements['FEPS']='FEP1.fep'
            replacements['RUNNUMBER'] = 'cp $inputfiles/eq*.inp .\nsed -i s/SEED_VAR/"$[1 + $[RANDOM % 9999]]"/ eq1.inp\n'
            replacements['FEPDIRR']='FEP1'

        if step == 'step2':
            replacements['FEPS']='FEP2.fep'
            replacements['RUNNUMBER'] = 'lastfep=FEP1\ncp $workdir/$lastfep/$temperature/$run/$finalMDrestart $rundir/eq5.re\n'
            replacements['FEPDIRR']='FEP2'

        if step == 'step3':
            replacements['FEPS']='FEP3.fep'
            replacements['RUNNUMBER'] = 'lastfep=FEP2\ncp $workdir/$lastfep/$temperature/$run/$finalMDrestart $rundir/eq5.re\n'
            replacements['FEPDIRR']='FEP3'

        if step == 'step4':
            replacements['FEPS']='FEP4.fep'
            replacements['RUNNUMBER'] = 'lastfep=FEP3\ncp $workdir/$lastfep/$temperature/$run/$finalMDrestart $rundir/eq5.re\n'
            replacements['FEPDIRR']='FEP4'

        if step == 'step5':
            replacements['FEPS']='FEP5.fep'
            replacements['RUNNUMBER'] = 'cp $inputfiles/eq*.inp .\nsed -i s/SEED_VAR/"$[1 + $[RANDOM % 9999]]"/ eq1.inp\n'
            replacements['FEPDIRR']='FEP5'
        
        with open(src) as infile, open(tgt, 'w') as outfile:
            for line in infile:
                if line.strip() == '#SBATCH -A ACCOUNT':
                    try:
                        replacements['ACCOUNT']
                        
                    except:
                        line = ''
                outline = IO.replace(line, replacements)
                outfile.write(outline)
                
                if line.strip() == '#EQ_FILES':
                    if step == 'step1' or step == 'step5':
                        for line in EQ_files:
                            file_base = line.split('/')[-1][:-4]
                            outline = 'time mpirun -np {} $qdyn {}.inp' \
                                    ' > {}.log\n'.format(ntasks,
                                                        file_base,
                                                        file_base)
                            outfile.write(outline)
                        
                if line.strip() == '#RUN_FILES':
                    for line in MD_files:
                        file_base = line.split('/')[-1][:-4]
                        outline = 'time mpirun -np {} $qdyn {}.inp'  \
                                    ' > {}.log\n'.format(ntasks,
                                                        file_base,
                                                        file_base)
                        outfile.write(outline)
                            

    def write_qfep(self, inputdir, windows, lambdas):
        qfep_in = s.ROOT_DIR + '/INPUTS/qfep.inp' 
        qfep_out = inputdir + 'qfep.inp'
        i = 0
        total_l = len(lambdas)
        
        # TO DO: multiple files will need to be written out for temperature range
        kT = f.kT(float(self.temperature))
        replacements = {}
        replacements['kT']=kT
        replacements['WINDOWS']=windows
        replacements['TOTAL_L']=str(total_l)
        with open(qfep_in) as infile, open(qfep_out, 'w') as outfile:
            for line in infile:
                line = run.replace(line, replacements)
                outfile.write(line)
            
            if line == '!ENERGY_FILES\n':
                for i in range(0, total_l):
                    j = -(i + 1)
                    lambda1 = lambdas[i]
                    lambda2 = lambdas[j]
                    filename = 'md_' +                          \
                                lambda1.replace('.', '') +      \
                                '_' +                           \
                                lambda2.replace('.', '') +      \
                                '.en\n'
                                
                    outfile.write(filename)
                    
def parseargs(args: list[str] = []) -> argparse.Namespace:
    """Return a Namespace after parsing an argument string.

    If args is not provided, defaults to args from command line.
    """
    parser = argparse.ArgumentParser(
        prog='QligFEP',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Generate FEP files for dual topology ligand FEP == ')

    parser.add_argument('-lig', '--ligand',
                        dest = "lig",
                        required = True,
                        help = "name of ligand 1")

    parser.add_argument('-FF', '--forcefield',
                        dest = "FF",
                        required = False,
                        choices = ['OPLS2005', 'OPLS2015', 'AMBER14sb', 'CHARMM36', 'CHARMM22', 'CHARMM_TEST'],
                        default = 'AMBER14sb',
                        help = "Forcefield to be used")

    parser.add_argument('-c', '--cluster',
                        dest = "cluster",
                        required = True,
                        help = "cluster you want to submit to, cluster specific parameters added to settings"
                       )

    parser.add_argument('-r', '--sphereradius',
                        dest = "sphereradius",
                        required = False,
                        default = '15',
                        help = "size of the simulation sphere"
                       )

    parser.add_argument('-b', '--cysbond',
                        dest = "cysbond",
                        default = None,
                        help = "Function to add cysbonds given as at1_at2,at3_at4,(...)"
                       )

    parser.add_argument('-T', '--temperature',
                        dest = "temperature",
                        default = '298',
                        help = "Temperature(s), mutliple tempereratures given as 'T1,T2,...,TN'"
                       )

    parser.add_argument('-R', '--replicates',
                        dest = "replicates",
                        default = '10',
                        help = "How many repeats should be run"
                       )

    parser.add_argument('-S1', '--sampling1',
                        dest = "sampling1",
                        default = 'linear',
                        choices = ['linear', 'sigmoidal', 'exponential', 'reverse_exponential'],
                        help = "Lambda spacing type to be used for removal of coulombic interactions"
                       )

    parser.add_argument('-S2', '--sampling2',
                        dest = "sampling2",
                        default = 'linear',
                        choices = ['linear', 'sigmoidal', 'exponential', 'reverse_exponential'],
                        help = "Lambda spacing type to be used for addition of softcore"
                       )

    parser.add_argument('-S3', '--sampling3',
                        dest = "sampling3",
                        default = 'linear',
                        choices = ['linear', 'sigmoidal', 'exponential', 'reverse_exponential'],
                        help = "Lambda spacing type to be used to remove solvent waters vdW and softcore"
                       )

    parser.add_argument('-S4', '--sampling4',
                        dest = "sampling4",
                        default = 'linear',
                        choices = ['linear', 'sigmoidal', 'exponential', 'reverse_exponential'],
                        help = "Lambda spacing type to be used to perform RBFE with waters"
                       )

    parser.add_argument('-S5', '--sampling5',
                        dest = "sampling5",
                        default = 'linear',
                        choices = ['linear', 'sigmoidal', 'exponential', 'reverse_exponential'],
                        help = "Lambda spacing type to be used"
                       )

    parser.add_argument('-w1', '--windows1',
                        dest = "windows1",
                        default = '50',
                        help = "Total number of windows that will be run for removal of coulombic interactions"
                       )

    parser.add_argument('-w2', '--windows2',
                        dest = "windows2",
                        default = '50',
                        help = "Total number of windows that will be run for addition of softcore"
                       )

    parser.add_argument('-w3', '--windows3',
                        dest = "windows3",
                        default = '50',
                        help = "Total number of windows that will be run to remove solvent waters vdW and softcore"
                       )

    parser.add_argument('-w4', '--windows4',
                        dest = "windows4",
                        default = '50',
                        help = "Total number of windows that will be run to perform RBFE with waters"
                       )

    parser.add_argument('-w5', '--windows5',
                        dest = "windows5",
                        default = '50',
                        help = "Total number of windows that will be run to calculate vacuum annihilation"
                       )

    if args:
        return parser.parse_args(args)
    else:
        return parser.parse_args()


if __name__ == "__main__":
    args = parseargs()
    run = Run(lig = args.lig,
              FF= args.FF,
              cluster = args.cluster,
              sphereradius = args.sphereradius,
              cysbond = args.cysbond,
              temperature = args.temperature,
              replicates = args.replicates,
              sampling1 = args.sampling1,
              sampling2 = args.sampling2,
              sampling3 = args.sampling3,
              sampling4 = args.sampling4,
              sampling5 = args.sampling5
              )

    writedir = run.makedir()
    prepinputdir = writedir + '/prep_input/'
    inputdir1 = writedir + '/inputfiles1/'
    inputdir2 = writedir + '/inputfiles2/'
    inputdir3 = writedir + '/inputfiles3/'
    inputdir4 = writedir + '/inputfiles4/'
    inputdir5 = writedir + '/inputfiles5/'
    # somthing like inputdir + '/FEP1' etc
    # save all FEPs in seperate inputfiules

    ligmolsize, charges, atomtypes, merged_molsize = run.read_files()
    FEP_vdw = run.change_prm(charges, writedir)

    # maybe place these copys somewhere at the end of a function
    run.copyfiles()

    run.prepare_complex('prepare', prepinputdir)
    run.qprep(prepinputdir)

    # make new initialize file
    solvent_oxygens = run.get_solvent_oxygens(prepinputdir)
    run.remove_waters(solvent_oxygens, prepinputdir, inputdir1)
    lig_oxygens, overlap_atomnumbers_lig = run.retreive_lig_coordinates(inputdir1)
    all_waters =  lig_oxygens + solvent_oxygens
    run.write_waters_to_pdb(all_waters, inputdir1)

    run.prepare_complex('water', inputdir1)
    run.qprep(inputdir1)
    run.copyfiles2()

    run.prepare_complex('vacuum', inputdir5)
    run.qprep(inputdir5)

    # write FEP files
    all_atoms_FEP = run.write_FEP1_file(charges, atomtypes, FEP_vdw, inputdir1, ligmolsize)
    run.write_FEP2_file(charges, atomtypes, FEP_vdw, inputdir2, ligmolsize)
    run.write_FEP3_file(charges, atomtypes, FEP_vdw, inputdir3, ligmolsize)
    run.write_FEP4_file(charges, atomtypes, FEP_vdw, inputdir4, ligmolsize)
    run.write_FEP5_file(charges, atomtypes, FEP_vdw, inputdir5, ligmolsize)

    # get lambdas
    lambdas1 = run.get_lambdas(args.windows1, args.sampling1)
    lambdas2 = run.get_lambdas(args.windows2, args.sampling2)
    lambdas3 = run.get_lambdas(args.windows3, args.sampling3)
    lambdas4 = run.get_lambdas(args.windows4, args.sampling4)
    lambdas5 = run.get_lambdas(args.windows5, args.sampling5)

    # write MD files
    run.write_MD(lambdas1, inputdir1, all_atoms_FEP, overlap_atomnumbers_lig, ligmolsize, 'step1')
    run.write_MD(lambdas2, inputdir2, all_atoms_FEP, overlap_atomnumbers_lig, ligmolsize, 'step2')
    run.write_MD(lambdas3, inputdir3, all_atoms_FEP, overlap_atomnumbers_lig, ligmolsize, 'step3')
    run.write_MD(lambdas4, inputdir4, all_atoms_FEP, overlap_atomnumbers_lig, ligmolsize, 'step4')
    run.write_MD(lambdas5, inputdir5, all_atoms_FEP, overlap_atomnumbers_lig, ligmolsize, 'step5')

    run.write_runfile(inputdir1, 'step1')
    run.write_runfile(inputdir2, 'step2')
    run.write_runfile(inputdir3, 'step3')
    run.write_runfile(inputdir4, 'step4')
    run.write_runfile(inputdir5, 'step5')

    run.write_submitfile(writedir)

    run.write_qfep(inputdir1, args.windows1, lambdas1)
    run.write_qfep(inputdir2, args.windows2, lambdas2)
    run.write_qfep(inputdir3, args.windows3, lambdas3)
    run.write_qfep(inputdir4, args.windows4, lambdas4)
    run.write_qfep(inputdir5, args.windows5, lambdas5)
