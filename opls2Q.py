import settings as s
from subprocess import check_output
import shlex
import math
import re
import sys
import os
import argparse

import IO

class Run(object):
    """
    Create FEP files from a common substructure for a given set of
    ligands
    """
    def __init__(self, lig, FF, output, merge, *args, **kwargs):
        """
        The init method is a kind of constructor, called when an instance
        of the class is created. The method serves to initialize what you
        want to do with the object.
        """
        supported_FFs = ['OPLS2005', 'OPLS2015']
        programs = ['Q']#, 'GROMACS']
        self.lig = lig
        self.FF = FF
        self.merge = merge
        self.ff_list = []
        
        if output in programs:
            self.output = output
        else:
            print('Conversion to ' + ','.join(programs)) + ' currently supported'
            sys.exit()
            
        if self.FF in supported_FFs:
            self.FF = FF
        else:
            print('Forcefields ' + ','.join(supported_FFs)) + ' currently supported'
            sys.exit()
        
    def vdw_calc(self, sig, eps):
        # calculate LJ types
        sig = float(sig)
        eps = float(eps)
        Avdw1 = math.sqrt(4*eps*(sig**12))
        Avdw2 = math.sqrt(4*eps*(sig**12))
        Bvdw1 = math.sqrt(4*eps*(sig**6))
        Avdw3 = (math.sqrt(4*eps*(sig**12)))/math.sqrt(2)
        Bvdw23 = (math.sqrt(4*eps*(sig**6)))/math.sqrt(2)

        # returns unparsed floats, consider changing
        return[Avdw1, Avdw2, Bvdw1, Avdw3, Bvdw23]

    def get_mass(self, atom):
        mass_dic = {"H":"1.0080",
                    "C":"12.0110",
                    "N":"14.0070",
                    "O":"15.9994",
                    "F":"19.00",
                    "P":"30.97",
                    "S":"32.0600",
                    "Cl":"35.0000",
                    "Br":"79.90",
                    "I":"126.90", 
                    "DUM": "0.00"
                   }
        at = re.findall('\d+|\D+', atom)
        mass = mass_dic[at[0]]
        return mass

    def bond_calc(self, k):
        k = float(k)*2
        return k

    # Double function, but this might be needed for different forcefields
    def angle_calc(self, k):
        k = float(k)*2
        return k


    # FIX
    def torsion_calc(self, tors):
        cnt = 0
        tors_Q = []
        for line in tors:
            cnt = cnt + 1
            if cnt == 1 or cnt == 3:
                Vn = (float(line))/2
                Vn_Q = [Vn, -cnt, '0.000', '1.000']
                tors_Q.append(Vn_Q)
            if cnt == 2:
                Vn = (float(line))/2
                Vn_Q = [Vn, -cnt, '180.000', '1.000']
                tors_Q.append(Vn_Q)
            if cnt == 4:
                Vn = (float(line))/2
                Vn_Q = [Vn, cnt, '180.000', '1.000']
                tors_Q.append(Vn_Q)

            #add some debugging here (if n tors > 4)
            else:
                continue
        return tors_Q

    def improper_calc(self, imp_V2):
        imp_V2 = float(imp_V2)/2

        return[imp_V2]



    # This here should later input the prefererd forcefield 
    # defined by user and generate outputfile

    def get_charge_groups(self, charges, bonds):
        dic = {}
        charge_group = []
        charge_group_list = []
        charge_group_dic = {}

        # Generating bond dictionary
        for line in bonds:
            if line[0] in dic:
                dic[line[0]].append(line[1])
            else:
                dic[line[0]] = [line[1]]

        # Now first we do a naive check for X-H charge group partners
        for key in dic:
            for line in dic[key]:
                charge = float(charges[key]) + float(charges[line])
                if charge.is_integer() == True:
                    charge_group_list.append([key, line])
                    charge_group_dic[key] = line
                    charge = 0

        charge = 0
        # Now for larger groups
        for line in bonds:
            if line[0] in charge_group_dic or line[1] in charge_group_dic:
                continue

            else:
                charge_a1 = float(charges[line[0]])
                charge_a2 = float(charges[line[1]])

                if line[0] in charge_group:
                    charge = round((charge + charge_a2), 5)
                    charge_group.append(line[1])

                else:
                    charge = charge + charge_a1 + charge_a2
                    charge_group = [line[0], line[1]]

                if charge.is_integer() == True:
                    charge_group_list.append(charge_group)
                    charge_group = []
                    charge = 0

        return charge_group_list                        

    def get_parameters(self):
        # Add other parameter generators later
        v = 'OPLS 14'
        v = v.split(' ')
        # Running command line tool has been moved to IO, change function!
        if v[0] == 'OPLS':
            ffld_serv = s.SCHROD_DIR + '/utilities/ffld_server'
            options = ' -ipdb '+ self.lig + '.pdb -print_parameters -version ' + v[1]
            args = shlex.split(ffld_serv + options)
            out = check_output(args,universal_newlines=True)
        
        elif v[0] == 'amber':
            print('not supported yet')
            sys.exit()

        with open(self.lig+'.log', 'w') as outfile:
            # added tag??
            outfile.write('Forcefield: ' + v[0] + ' ' + v[1] + '\n')
            for line in out:
                outfile.write(line)

    def read_log(self):
        opls_list = []    
        # Read output from ffld_server ## FIX FOR AMBER
        with open(self.lig + '.log', 'r') as infile:
            for line in infile:
                if len(line) > 2:
                    opls_list.append(line.split())

        self.ff_list = opls_list


    def convert_toQ(self):
        ff_list = self.ff_list
        halogen_cnt = 0

        charge_sum = 0
        charge_list = []
        charge_dic = {}
        vdw_list = []
        bond_list = []
        bonded = []
        charge_group_list = []
        angle_list = []
        torsion_list = []
        improper_list = []

        block = 0

        # READ BLOCKS
        for line in ff_list:
            if len(line) > 1:
                if line[0] == "OPLSAA":
                    block = 1
                if line[0] == "Stretch":
                    block = 2     
                if line[0] == "Bending":
                    block = 3
                if line[0] == "proper":
                    block = 4
                if line[0] == "improper":
                    block = 5

            # Create lists 
            if len(line[0]) <= 4 and line[0] != 'atom':
                if block == 1:
                    charge = [line[0], float(line[4])]
                    charge_dic[line[0]] = line[4]
                    charge_sum = charge_sum + float(line[4])
                    vdw = [line[0], run.vdw_calc(line[5], line[6]) + [run.get_mass(line[0])]]
                    charge_list.append(charge)
                    vdw_list.append(vdw)

                if block == 2:
                    bond = [line[0:2], [run.bond_calc(line[2]), line[3]]]
                    bond_list.append(bond)
                    bonded.append(line[0:2])

                if block == 3:
                    angle = [line[0:3], [run.angle_calc(line[3]), line[4]]]
                    angle_list.append(angle)

                if block == 4:
                    torsion = [line[0:4], run.torsion_calc(line[4:8])] 
                    torsion_list.append(torsion)

                if block == 5:
                    improper = [line[0:4], run.improper_calc(line[4]) + ['180.000']]
                    improper_list.append(improper)
        
        charge_group_list = run.get_charge_groups(charge_dic, bonded)

        self.Q_FF =[charge_list, 
                    charge_sum, 
                    vdw_list, 
                    bond_list, 
                    angle_list, 
                    torsion_list, 
                    improper_list,
                    charge_group_list]

    def write_lib_Q(self):
        # this is just for readability, not necessary
        parm = self.Q_FF
        at_len = len(parm[0])
        charge = parm[1]
        charges = parm[0]
        bonds = parm[3]
        improper = parm[6]
        charge_groups = parm[7]

        with open(self.lig + '.lib', 'w') as outfile:
            outfile.write("{LIG}     ! atoms no %4d   total charge %.3f \n\n" % (at_len,
                                                                               charge
                                                                              )
                         )
            outfile.write("[info] \n SYBYLtype RESIDUE \n\n")

            # atom and charge block:
            outfile.write("[atoms] \n")
            cnt = 0
            for line in charges:
                cnt = cnt + 1
                outfile.write('{:>4d}   {:10}X{:11}{: .4f}\n'.format(cnt, 
                                                                     line[0], 
                                                                     line[0], 
                                                                     line[1]
                                                                    )
                             )

            outfile.write("\n[bonds]\n")
            for line in bonds:
                outfile.write('{:10}{}\n'.format(line[0][0], line[0][1]))

            outfile.write("\n[impropers]\n")
            for line in improper:
                outfile.write('{:10}{:10}{:10}{}\n'.format(line[0][0], 
                                                           line[0][1], 
                                                           line[0][2], 
                                                           line[0][3]))

            outfile.write("\n[charge_groups]\n")
            for i in charges:
                if i[0][0] != 'H':
                    outfile.write('{}'.format(i[0]))
                    for j in bonds:
                        if j[0][0]==i[0] and j[0][1][0] == 'H':
                            outfile.write(' {} '.format(j[0][1]))
                        if j[0][1][0] == i[0] and j[0][1][0] =='H':
                            outfile.write(" H%i" % j[0])
                    outfile.write("\n")

    def write_prm_Q(self):
        parm = self.Q_FF
        vdw = parm[2]
        bond = parm[3]
        angle = parm[4]
        torsion = parm[5]
        improper = parm[6]

        if self.FF == 'OPLS2015' and self.merge == True:
            prm_file = os.path.join(s.FF_DIR, 'OPLS2015.prm')
            prm_file_out = self.FF + '_' + self.lig + '.prm'
            
        elif self.FF == 'OPLS2005' and self.merge == True:
            prm_file = os.path.join(s.FF_DIR, 'OPLS2005.prm')
            prm_file_out = self.FF + '_' + self.lig + '.prm'
            
        elif self.merge == False:
            prm_file = os.path.join(s.FF_DIR, 'NOMERGE.prm')
            prm_file_out = self.lig + '.prm'
        
        with open(prm_file) as infile, open(prm_file_out, 'w') as outfile:
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

                # Create lists
                if block == 1: 
                    for line in vdw:
                        outfile.write("""X{:6}{: 8.2f}{: 10.2f}{: 10.2f}{: 10.2f}{: 10.2f}{:>10s}\n""".format(line[0], 
                                                                          line[1][0],
                                                                          line[1][1],
                                                                          line[1][2],
                                                                          line[1][3],
                                                                          line[1][4],
                                                                          line[1][5]))

                if block == 2:
                    for line in bond:
                        outfile.write('X{:10}X{:10}{:10.1f}{:>10.5s}\n'.format(line[0][0], 
                                                                     line[0][1],
                                                                     line[1][0],
                                                                     line[1][1]))

                if block == 3:
                    for line in angle:
                        outfile.write("""X{:10}X{:10}X{:10}{: 8.2f}{:>12.7}\n""".format(line[0][0],
                                                                   line[0][1],
                                                                   line[0][2],
                                                                   line[1][0],
                                                                   line[1][1]))

                if block == 4:
                    for line in torsion:
                        for line2 in line[1]:
                            outfile.write("""X{:10}X{:10}X{:10}X{:10}{:<10.3f}{: d}.000{:>10s}{:>10s}\n""".format(line[0][0],
                                                                          line[0][1],
                                                                          line[0][2],
                                                                          line[0][3],
                                                                          line2[0],
                                                                          line2[1],
                                                                          line2[2],
                                                                          line2[3]))

                if block == 5:
                    for line in improper:
                        outfile.write("""X{:10}X{:10}X{:10}X{:10}{:10.3f}{:>10s}\n""".format(line[0][0],
                                                      line[0][1],
                                                      line[0][2],
                                                      line[0][3],
                                                      line[1][0],
                                                      line[1][1]))

    def write_itp_GROMACS(self):
        return True

    def rename_pdb(self, include = ('ATOM', 'HETATM')):
        pdb_in = self.lig + '.pdb'
        pdb_out = self.lig + '_out.pdb'
        index = -1
        atomnames = self.Q_FF[0]
        with open(pdb_in) as infile, open(pdb_out, 'w') as outfile:
            for line in infile:
                line = IO.pdb_parse_in(line)
                if line[0].strip() in include:
                    line[4] = 'LIG'
                    index += 1
                    line[2] = atomnames[index][0]
                    line[6] = 1
                    outline = IO.pdb_parse_out(line)
                    outfile.write(outline + '\n')
                    
                print(len(line))
                    
        os.rename(pdb_out, pdb_in)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='lig_prep',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Generate parameter files for ligands. == ')

    
    parser.add_argument('-l', '--ligand',
                        dest = "lig",
                        required = True,
                        help = "name of the ligand")
    
    parser.add_argument('-FF', '--forcefield',
                    dest = "FF",
                    required = True,
                    help = "forcefield to use, OPLS2005 or OPLS2015")
    
    parser.add_argument('-o', '--output',
                    dest = "output",
                    required = True,
                    help = "generate parameter files for Q or GROMACS")

    parser.add_argument('-m', '--merge',
                        dest = "merge",
                        action = 'store_false',
                        default = True,        
                        help = "Use this flag if you do not want the ligand prms to be merged")
    
    args = parser.parse_args()
    run = Run(lig = args.lig,
              FF = args.FF,
              output = args.output,
              merge = args.merge
             )

    run.get_parameters()
    run.read_log()
    run.convert_toQ()
    run.write_lib_Q()
    run.write_prm_Q()
    run.rename_pdb()
