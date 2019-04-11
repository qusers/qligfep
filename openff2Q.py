import settings as s
from subprocess import check_output
import shlex
import math
import re
import sys
import os
import argparse

# QligFEP modules
import IO

# openFF modules
from openforcefield.topology import Molecule, Topology
from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.utils import get_data_filename

class Run(object):
    """
    Create FEP files from a common substructure for a given set of
    ligands
    """
    def __init__(self, lig, FF, merge, *args, **kwargs):
        """
        The init method is a kind of constructor, called when an instance
        of the class is created. The method serves to initialize what you
        want to do with the object.
        """
        self.lig = lig
        self.FF = FF
        self.merge = merge
        self.ff_list = []
        self.mapping = {}
        self.total_charge = 0
        masses =   {"H"     : "1.0080",
                    "C"     : "12.0110",
                    "N"     : "14.0070",
                    "O"     : "15.9994",
                    "F"     : "19.00",
                    "P"     : "30.97",
                    "S"     : "32.0600",
                    "Cl"    : "35.0000",
                    "Br"    : "79.90",
                    "I"     : "126.90", 
                    "DUM"   : "0.00"
                   }
    
    def openff(self):
        # Load the molecule (for now mol2, until charges are saved on sdf)
        molecule = Molecule.from_file(self.lig + '.mol2')
        topology = Topology.from_molecules([molecule])

        # Label using the smirnoff99Frosst force field
        self.forcefield = ForceField('smirnoff99Frosst.offxml')
        self.parameters = self.forcefield.label_molecules(topology)[0]
    
    def read_mol2(self):
        """
            This is basically to get the charge, will later be deprecated when charges are
            transferable in openff
        """
        with open(self.lig + '.mol2') as infile:
            cnt = -1
            for line in infile:
                line = line.split()
                if len(line) == 9:
                    cnt += 1
                    self.mapping[cnt] = [line[0],               # at idex
                                         line[1],               # atname
                                         line[5].split('.')[0], # attype
                                         line[8]                # charge
                                        ]
                    self.total_charge += float(line[8]) 
        
        if self.total_charge != 0.0:
            print('WARNING: residual charge {} check your mol2 file!'.format(self.total_charge))
            
    def write_lib_Q(self):
        with open(self.lig + '.lib', 'w') as outfile:
            outfile.write('{}    ! atoms no {}   total charge {} \n\n'.format('{LIG}',
                                                                                 len(self.mapping),
                                                                                 self.total_charge)
                         )
            
            outfile.write("[info] \n SYBYLtype RESIDUE \n\n")

            #atom and charge block:
            outfile.write("[atoms] \n")
            for i, at in enumerate(self.mapping):
                outfile.write('{:>4s}   {:10}X{:11}{:>10s}\n'.format(self.mapping[at][0], 
                                                                     self.mapping[at][1], 
                                                                     self.mapping[at][1], 
                                                                     self.mapping[at][3]
                                                                    )
                             )
                
            # bonded block
            outfile.write("\n[bonds]\n")
            for i, bond in enumerate(self.parameters['Bonds']):
                ai = self.mapping[bond[0]][1]
                aj = self.mapping[bond[1]][1]
                outfile.write('{:10s}{:}\n'.format(ai,
                                                   aj))
            
            # improper block
            outfile.write("\n[impropers]\n")
            for i, torsion in enumerate(self.parameters['ImproperTorsions']):
                ai = self.mapping[torsion[0]][1]
                aj = self.mapping[torsion[1]][1]
                ak = self.mapping[torsion[2]][1]
                al = self.mapping[torsion[3]][1]
                outfile.write('{:10}{:10}{:10}{}\n'.format(ai, 
                                                           aj, 
                                                           ak, 
                                                           al))

            outfile.write("\n[charge_groups]")
            for i, atom in enumerate(self.mapping):
                if self.mapping[atom][2] != 'H':
                    outfile.write('\n{}'.format(self.mapping[atom][1]))
                for j, bond in enumerate(self.parameters['Bonds']):
                    if bond[0] == i:
                        if self.mapping[bond[1]][2] == 'H':
                            outfile.write(' {}'.format(self.mapping[bond[1]][1]))
            
    def write_prm_Q(self):
        if self.FF == 'AMBER14sb' and self.merge == True:
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
                    ## STUCK HERE
                    for i, vdw in enumerate(self.parameters['vdW']):
                        id = vdw.id
                        print(i, self.forcefield.get_handler('vdW').parameters[i])
                        #outfile.write("""X{:6}{: 8.2f}{: 10.2f}{: 10.2f}{: 10.2f}{: 10.2f}{:>10s}\n""".format(line[0], 
                        #                                                  line[1][0],
                        #                                                  line[1][1],
                        #                                                  line[1][2],
                        #                                                  line[1][3],
                        #                                                  line[1][4],
                        #                                                  line[1][5]))

                #if block == 2:
                #    for line in bond:
                #        outfile.write('X{:10}X{:10}{:10.1f}{:>10.5s}\n'.format(line[0][0], 
                #                                                     line[0][1],
                #                                                     line[1][0],
                #                                                     line[1][1]))

                #if block == 3:
                #    for line in angle:
                #        outfile.write("""X{:10}X{:10}X{:10}{: 8.2f}{:>12.7}\n""".format(line[0][0],
                #                                                   line[0][1],
                #                                                   line[0][2],
                #                                                   line[1][0],
                #                                                   line[1][1]))

                #if block == 4:
                #    for line in torsion:
                #        for line2 in line[1]:
                #            outfile.write("""X{:10}X{:10}X{:10}X{:10}{:<10.3f}{: d}.000{:>10s}{:>10s}\n""".format(line[0][0],
                #                                                          line[0][1],
                #                                                          line[0][2],
                #                                                          line[0][3],
                #                                                          line2[0],
                #                                                          line2[1],
                #                                                          line2[2],
                #                                                          line2[3]))

                #if block == 5:
                #    for line in improper:
                #        outfile.write("""X{:10}X{:10}X{:10}X{:10}{:10.3f}{:>10s}\n""".format(line[0][0],
                #                                      line[0][1],
                #                                      line[0][2],
                #                                      line[0][3],
                #                                      line[1][0],
                #                                      line[1][1]))    
    
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
                        choices = ['AMBER14sb'],
                        help = "forcefield to use")

    parser.add_argument('-m', '--merge',
                        dest = "merge",
                        action = 'store_false',
                        default = True,        
                        help = "Use this flag if you do not want the ligand prms to be merged")
    
    args = parser.parse_args()
    run = Run(lig = args.lig,
              FF = args.FF,
              merge = args.merge
             )

    run.openff()
    run.read_mol2()
    #run.convert_to_Q()
    run.write_lib_Q()
    run.write_prm_Q()
    #run.rename_pdb()
