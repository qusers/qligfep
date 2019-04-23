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
        self.masses =   {"H"     : "1.0080",
                         "C"     : "12.0110",
                         "N"     : "14.0070",
                         "O"     : "15.9994",
                         "F"     : "19.0000",
                         "P"     : "30.9700",
                         "S"     : "32.0600",
                         "Cl"    : "35.0000",
                         "Br"    : "79.9000",
                         "I"     : "126.90", 
                        "DUM"   : "0.0000"
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
                                         line[8],               # charge
                                         line[2],               # X coordinate
                                         line[3],               # Y coordinte
                                         line[4]                # Z coordinate
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
                outfile.write('{:>4s}   {:10}{:11}{:>10s}\n'.format(self.mapping[at][0], 
                                                                    self.mapping[at][1], 
                                                                    self.mapping[at][1].lower(), 
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

            #outfile.write("\n[charge_groups]")
            #for i, atom in enumerate(self.mapping):
            #    if self.mapping[atom][2] != 'H':
            #        outfile.write('\n{}'.format(self.mapping[atom][1]))
            #    for j, bond in enumerate(self.parameters['Bonds']):
            #        if bond[0] == i:
            #            if self.mapping[bond[1]][2] == 'H':
            #                outfile.write(' {}'.format(self.mapping[bond[1]][1]))
            
    def write_prm_Q(self):
        if self.FF == 'AMBER14sb' and self.merge == True:
            prm_file = os.path.join(s.FF_DIR, 'AMBER14sb.prm')
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

                if block == 1:
                    for (atom_indices, parameter) in self.parameters['vdW'].items():
                        ai = atom_indices[0]
                        ai_name = self.mapping[ai][1].lower()
                        # This is a bit hacky, check how to get the float out directly
                        epsilon = float('{}'.format(parameter.epsilon).split()[0])
                        epsilon23 = epsilon/2
                        # TO DO: CHECK IF THIS IS CORRECT!!
                        Rmin = float('{}'.format(parameter.sigma).split()[0])/2
                        mass = self.masses[self.mapping[ai][2]]
                        outfile.write("""{:6}{: 8.3f}{: 10.3f}{: 10.3f}{: 10.3f}{: 10.3f}{:>10s}\n""".format(ai_name, 
                                                                          Rmin,
                                                                          0.00,
                                                                          epsilon,
                                                                          Rmin,
                                                                          epsilon23,
                                                                          mass))

                if block == 2:
                    for (atom_indices, parameter) in self.parameters['Bonds'].items():
                        ai      = atom_indices[0]
                        ai_name = self.mapping[ai][1].lower()
                        aj      = atom_indices[1]
                        aj_name = self.mapping[aj][1].lower()
                        fc      = float('{}'.format(parameter.k).split()[0])
                        l       = float('{}'.format(parameter.length).split()[0])
                        outfile.write('{:10}{:10}{:10.1f}{:>10.3f}\n'.format(ai_name, 
                                                                               aj_name,
                                                                               fc,
                                                                               l))

                if block == 3:
                    for (atom_indices, parameter) in self.parameters['Angles'].items():
                        ai      = atom_indices[0]
                        ai_name = self.mapping[ai][1].lower()
                        aj      = atom_indices[1]
                        aj_name = self.mapping[aj][1].lower()
                        ak      = atom_indices[2]
                        ak_name = self.mapping[ak][1].lower()
                        fc      = float('{}'.format(parameter.k).split()[0])
                        angle   = float('{}'.format(parameter.angle).split()[0])
                        
                        outfile.write("""{:10}{:10}{:10}{: 8.2f}{:>12.3f}\n""".format(ai_name,
                                                                                         aj_name,
                                                                                         ak_name,
                                                                                         fc,
                                                                                         angle))

                if block == 4:
                    for (atom_indices, parameter) in self.parameters['ProperTorsions'].items():
                        forces = []
                        ai      = atom_indices[0]
                        ai_name = self.mapping[ai][1].lower()
                        aj      = atom_indices[1]
                        aj_name = self.mapping[aj][1].lower()
                        ak      = atom_indices[2]
                        ak_name = self.mapping[ak][1].lower()
                        al      = atom_indices[3]
                        al_name = self.mapping[al][1].lower()
                        max_phase = len(parameter.phase)
                        
                        # Now check if there are multiple minima
                        for i in range(0, max_phase):
                            fc      = float('{}'.format(parameter.k[i]).split()[0])
                            phase   = float('{}'.format(parameter.phase[i]).split()[0])
                            paths   = int(parameter.idivf[i])
                            
                            if i != max_phase-1 and max_phase > 1:
                                minimum  = float(parameter.periodicity[i])*-1
                                
                            else:
                                minimum  = float(parameter.periodicity[i])
                            
                            force = (fc,minimum,phase,paths)
                            forces.append(force)

                        for force in forces:
                            outfile.write("""{:10}{:10}{:10}{:10}{:>10.3f}{:>10.3f}{:>10.3f}{:>5d}\n""".format(ai_name,
                                                                                                                   aj_name,
                                                                                                                   ak_name,
                                                                                                                   al_name,
                                                                                                                   force[0],
                                                                                                                   force[1],
                                                                                                                   force[2],
                                                                                                                   force[3]))

                if block == 5:
                    for (atom_indices, parameter) in self.parameters['ImproperTorsions'].items():
                        ai      = atom_indices[0]
                        ai_name = self.mapping[ai][1].lower()
                        aj      = atom_indices[1]
                        aj_name = self.mapping[aj][1].lower()
                        ak      = atom_indices[2]
                        ak_name = self.mapping[ak][1].lower()
                        al      = atom_indices[3]
                        al_name = self.mapping[al][1].lower()                      
                        fc      = float('{}'.format(parameter.k[0]).split()[0])
                        phase   = float('{}'.format(parameter.phase[0]).split()[0])
                        outfile.write("""{:10}{:10}{:10}{:10}{:10.3f}{:10.3f}\n""".format(ai_name,
                                                                                              aj_name,
                                                                                              ak_name,
                                                                                              al_name,
                                                                                              fc,
                                                                                              phase))    
                        
    def write_PDB(self):
        with open(self.lig + '.pdb', 'w') as outfile:
            for atom in self.mapping:
                ai      = atom + 1
                ai_name = self.mapping[atom][1]
                a_el    = self.mapping[atom][2]
                ax      = float(self.mapping[atom][4])
                ay      = float(self.mapping[atom][5])
                az      = float(self.mapping[atom][6])
                at_entry = ['HETATM',   #  0 ATOM/HETATM
                            ai,         #  1 ATOM serial number
                            ai_name,    #  2 ATOM name
                            '',         #  3 Alternate location indicator
                            'LIG',      #  4 Residue name
                            '',         #  5 Chain identifier
                            1,          #  6 Residue sequence number
                            '',         #  7 Code for insertion of residue
                            ax,         #  8 Orthogonal coordinates for X
                            ay,         #  9 Orthogonal coordinates for Y
                            az,         # 10 Orthogonal coordinates for Z
                            0.0,        # 11 Occupancy
                            0.0,        # 12 Temperature factor
                            a_el,       # 13 Element symbol
                            ''          # 14 Charge on atom
                           ]
                outfile.write(IO.pdb_parse_out(at_entry) + '\n')
        
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
    run.write_lib_Q()
    run.write_prm_Q()
    run.write_PDB()
