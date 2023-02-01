import settings as s
from subprocess import check_output
import shlex
import math
import re
import sys
import os
import argparse
from copy import deepcopy
import numpy as np
from simtk import unit
import openmm

# QligFEP modules
import IO

# openFF modules
from openff.toolkit import Molecule, ForceField, Topology
from openff.interchange import Interchange
#from openforcefield.utils import get_data_filename

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
        self.charges_list_magnitude = []
        self.FF = FF
        self.merge = merge
        self.ff_list = []
        self.mapping = {}
        self.total_charge = 0
        self.molecule = Molecule.from_file(self.lig + '.sdf')
        self.topology = Topology.from_molecules([self.molecule])
        self.forcefield = ForceField('openff-1.0.0.offxml')
        self.parameters = self.forcefield.label_molecules(self.topology)[0]
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
    
    def charges_and_mapping(self):

        ### Assign partial charges and return in form of dictionary (can only load sdf, mol or smi) ###
        # sage = ForceField(self.forcefield)
        # molecule = Molecule.from_file(self.lig + '.sdf')
        self.partial_charges = Interchange.from_smirnoff(self.forcefield, [self.molecule])['Electrostatics'].charges

        # loop through dictionary and safe charges in list
        for charge in self.partial_charges.values():
            self.charges_list_magnitude.append(charge._magnitude)

        # get total charge and round
        self.total_charge = round(sum(self.charges_list_magnitude),10)

        if self.total_charge != 0.0:
            print('WARNING: residual charge {} check your mol2 file!'.format(self.total_charge))

        # loop through sdf lig file and safe relevant information in mapping dictionary
        with open(self.lig + '.sdf') as infile:
            cnt = -1
            # loop through sdf file and retreive information
            for line in infile:
                line = line.split()
                if len(line) == 16:
                    cnt += 1
                    atomindex = cnt + 1
                    charge = self.charges_list_magnitude[cnt]
                    self.mapping[cnt] = [   str(atomindex),                         # at idex
                                            line[3] + str(atomindex),               # atname
                                            line[3],                                # attype
                                            str(charge),                            # charge
                                            line[0],                                # X coordinate
                                            line[1],                                # Y coordinte
                                            line[2]                                 # Z coordinate
                                        ]
            
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

            # This was commented out?

            #outfile.write("\n[charge_groups]")
            #for i, atom in enumerate(self.mapping):
            #    if self.mapping[atom][2] != 'H':
            #        outfile.write('\n{}'.format(self.mapping[atom][1]))
            #    for j, bond in enumerate(self.parameters['Bonds']):
            #        if bond[0] == i:
            #            if self.mapping[bond[1]][2] == 'H':
            #                outfile.write(' {}'.format(self.mapping[bond[1]][1]))
            
    def write_prm_Q(self):

        # dont know what happens here
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
                        # TO DO: CHECK IF THIS IS CORRECT!
                        Rmin = '{}'.format(parameter.rmin_half)
                        Rmin = Rmin.split()[0]
                        Rmin = float(Rmin)
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

    def report_missing_parameters(self):
        """
        Analyze a molecule using a provided ForceField, generating a report of any
        chemical groups in the molecule that are lacking parameters.

        Parameters
        ----------
        molecule : an openforcefield.topology.FrozenMolecule
            The molecule to analyze
        forcefield : an openforcefield.typing.engine.smirnoff.ForceField
            The ForceField object to use

        Returns
        -------
        missing_parameters : dict[tagname: list[dict[tagged_smiles:string, image:PIL.Image, atom indices:list[int]]]]
            A hierarchical dictionary, with first level keys indicating ForceField tag 
            names (eg. "Bonds"), and first-level values which are lists of dictionaries. 
            Each dictionary in this list reflects one missing parameter, and contains the 
            following key:value pairs :
            * "image": PIL.Image 
                * shows a 2D drawing, highlighting the feature that could not be parametrized
            * "tagged_smiles": string
                * SMILES of the whole molecule, tagging the atom indices which could not be
                  parametrized
            * "atom_indices": tuple(int)
                * The indices of atoms which could not be parametrized

        """
        highlight_color = (0.75, 0.75, 0.75)

        # Make deepcopies of both inputs, since we may modify them in this function
        forcefield = deepcopy(self.forcefield)

        ###### lig.sdf file ######
        molecule = deepcopy(self.molecule)

        # Set partial charges to placeholder values so that we can skip AM1-BCC 
        # during parameterization

        # the partial charges have to be added to the molecule somehow with an array like thing multiplied with the elementary charge
        # They have to be the same format. It works for the molecule.partial_charges = (np.zeros(molecule.n_atoms) + 0.1) * unit.elementary_charge
        # molecule.n_atoms just gives 45 (number of atoms). So that should not be 
        molecule.partial_charges = np.array(self.charges_list_magnitude) * unit.elementary_charge

        # Prepare dictionary to catch parameterization failure info
        success = False
        missing_params = {}

        while not success:
            # Try to parameterize the system, catching the exception if there is one.
            try:

                ###### error occurs, Molecule Brc1cccc(Nc2nc(OCC3CCCCC3)c3nc[nH]c3n2)c1 has a net charge of 4.5 ######
                    #### why do you calculate charges when you have them already? ####
                    #### why molecule.to_topology() when we have that already in the openff function? ####
                    #### charge_from_molecules (List[openff.toolkit.molecule.Molecule], optional) – If specified, partial #
                    #    charges will be taken from the given molecules instead of being determined by the force field. #####

                forcefield.create_openmm_system(molecule.to_topology(), charge_from_molecules=[molecule])
                success = True

            ###### NameError: name 'UnassignedValenceParameterException' is not defined ######
            except UnassignedValenceParameterException as e:      
                success = False

                # Ensure that there is a list initialized for missing parameters 
                # under this tagname
                handler_tagname = e.handler_class._TAGNAME
                if handler_tagname not in missing_params:
                    missing_params[handler_tagname] = []

                # Create a shortcut to the topology atom tuples attached to 
                # the parametrization error
                top_atom_tuples =  e.unassigned_topology_atom_tuples

                # Make a summary of the missing parameters from this attempt and add it to 
                # the missing_params dict
                rdmol = molecule.to_rdkit()
                for top_atom_tuple in top_atom_tuples:
                    orig_atom_indices = [i.topology_atom_index for i in top_atom_tuple]
                    # Make a copy of the input RDMol so that we don't modify the original
                    this_rdmol = deepcopy(rdmol)

                    # Attach tags to relevant atoms so that a tagged SMILES can be written
                    orig_rdatoms = []
                    for tag_idx, atom_idx in enumerate(orig_atom_indices):
                        rdatom = this_rdmol.GetAtomWithIdx(atom_idx)
                        rdatom.SetAtomMapNum(tag_idx + 1)
                        orig_rdatoms.append(rdatom)    

                    tagged_smiles = Chem.MolToSmiles(this_rdmol)

                    # Make tagged hydrogens into deuteriums so that RemoveHs doesn't get rid of them
                    for rdatom in orig_rdatoms:
                        if rdatom.GetAtomicNum() == 1:
                            rdatom.SetIsotope(2)

                    # Remove hydrogens, since they clutter up the 2D drawing 
                    # (tagged Hs are not removed, since they were converted to deuterium)
                    h_less_rdmol = Chem.RemoveHs(this_rdmol)

                    # Generate 2D coords, since drawing from 3D can look really weird
                    Draw.rdDepictor.Compute2DCoords(h_less_rdmol)

                    # Search over the molecule to find the indices of the tagged atoms
                    # after hydrogen removal
                    h_less_atom_indices = [None for i in orig_atom_indices]
                    for rdatom in h_less_rdmol.GetAtoms():
                        # Convert deuteriums back into hydrogens
                        if rdatom.GetAtomicNum() == 1:
                            rdatom.SetIsotope(1)

                        atom_map_num = rdatom.GetAtomMapNum()
                        if atom_map_num == 0:
                            continue
                        h_less_atom_indices[atom_map_num-1] = rdatom.GetIdx()

                    # Once the new atom indices are found, use them to find the H-less
                    # bond indices
                    h_less_rdbonds = []
                    for i in range(len(h_less_atom_indices)-1):
                        rdbond = h_less_rdmol.GetBondBetweenAtoms( 
                                                h_less_atom_indices[i], 
                                                h_less_atom_indices[i+1])
                        h_less_rdbonds.append(rdbond)
                    h_less_bond_indices = [bd.GetIdx() for bd in h_less_rdbonds]

                    # Create a 2D drawing of the molecule, highlighting the 
                    # parameterization failure
                    highlight_atom_colors = {idx:highlight_color for idx in h_less_atom_indices}
                    highlight_bond_colors = {idx:highlight_color for idx in h_less_bond_indices}
                    image = Draw.MolsToGridImage([h_less_rdmol],
                                                 highlightAtomLists=[h_less_atom_indices],
                                                 highlightBondLists=[h_less_bond_indices],
                                                 molsPerRow=1,
                                                 highlightAtomColors=[highlight_atom_colors],
                                                 highlightBondColors=[highlight_bond_colors],
                                                 subImgSize=(600,600)
                                                )

                    # Structure and append the relevant info to the missing_params dictionary
                    param_description = {'atom_indices': orig_atom_indices,
                                         'image': image, 
                                         'tagged_smiles': tagged_smiles
                                        }
                    missing_params[handler_tagname].append(param_description)

                # Add a "super generic" parameter to the top of this handler's ParameterList, 
                # which will make it always find parameters for each term. This will prevent the same
                # parameterization exception from being raised in the next attempt.
                param_list = forcefield.get_parameter_handler(handler_tagname).parameters
                param_list.insert(0, super_generics[handler_tagname])

        if success != True:
            print(missing_params)
        else:
            print('Parameters succesfully assigned')
        
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

    # changed the default to false
    parser.add_argument('-m', '--merge',
                        dest = "merge",
                        action = 'store_false',
                        default = False,        
                        help = "Use this flag if you do not want the ligand prms to be merged")
    
    args = parser.parse_args()
    run = Run(lig = args.lig,
              FF = args.FF,
              merge = args.merge
             )

    run.charges_and_mapping()
    run.report_missing_parameters()
    run.write_lib_Q()
    run.write_prm_Q()
    run.write_PDB()
