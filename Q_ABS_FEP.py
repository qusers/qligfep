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
                 sampling,
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
        self.sampling = sampling
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
        
        ligand = rdmolfiles.MolFromPDBFile(prepdir + self.lig + '.pdb')

        for atom in ligand.GetAtoms():
            atomindex = atom.GetIdx()
            if atomindex % 2:
                continue
            else:
                coordinates = ligand.GetConformer().GetAtomPosition(atomindex)
                oxygencoords.append([format(round(coordinates.x + 0.001,3), '.3f'),
                format(round(coordinates.y + 0.001,3), '.3f'),
                format(round(coordinates.z + 0.001,3), '.3f')])
        return oxygencoords


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
                                                            line[1],
                                                            '0.000'))

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
                                                            line[1],
                                                            line[1]))

                # solute waters
                else:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            'DUM',
                                                            'DUM'))


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
                                                            line[1],
                                                            line[1]))

                # solute waters
                else:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            'DUM',
                                                            'DUM'))


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
                                                            '0'))

                # solute waters
                else:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            '20',
                                                            '20'))


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
                                                            line[1],
                                                            'DUM'))

                # solute waters
                else:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            'DUM',
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
                                                            '0.000'))

                # solute waters
                else:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            '0.000',
                                                            line[1]))

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
                                                            '0',
                                                            '0'))

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
                                                            'DUM'))

                # waters from solvent
                elif atomnumber <= ligsize + self.waters_to_perturb * 3:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            'DUM',
                                                            'DUM'))

                # solute waters
                else:
                    outfile.write("{:<5}{:>10}{:>10}\n".format(line[0],
                                                            'DUM',
                                                            line[1]))

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

    parser.add_argument('-S', '--sampling',
                        dest = "sampling",
                        default = 'linear',
                        choices = ['linear', 'sigmoidal', 'exponential', 'reverse_exponential'],
                        help = "Lambda spacing type to be used"
                       )
                        ### need 5 windows? ###
    parser.add_argument('-w', '--windows',
                        dest = "windows",
                        default = '50',
                        help = "Total number of windows that will be run"
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
              sampling = args.sampling
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
    lig_oxygens = run.retreive_lig_coordinates(inputdir1)
    all_waters = solvent_oxygens + lig_oxygens
    run.write_waters_to_pdb(all_waters, inputdir1)

    run.prepare_complex('water', inputdir1)
    run.qprep(inputdir1)
    run.copyfiles2()

    run.prepare_complex('vacuum', inputdir5)
    run.qprep(inputdir5)

    # write FEP files
    run.write_FEP1_file(charges, atomtypes, FEP_vdw, inputdir1, ligmolsize)
    run.write_FEP2_file(charges, atomtypes, FEP_vdw, inputdir2, ligmolsize)
    run.write_FEP3_file(charges, atomtypes, FEP_vdw, inputdir3, ligmolsize)
    run.write_FEP4_file(charges, atomtypes, FEP_vdw, inputdir4, ligmolsize)
    run.write_FEP5_file(charges, atomtypes, FEP_vdw, inputdir5, ligmolsize)

    # lambdas = run.get_lambdas(args.windows, args.sampling)

    # file_list = run.write_MD_1(lambdas, inputdir, lig_size1, lig_size2, overlapping_atoms)
    # run.write_runfile(inputdir, file_list)    
    
    # run.write_submitfile(writedir)
    # run.write_qfep(inputdir, args.windows, lambdas)
    # run.write_qprep(inputdir)
    # run.qprep(inputdir)