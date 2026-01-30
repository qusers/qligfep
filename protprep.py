#!/usr/bin/env python

import argparse
import os, sys
from time import gmtime, strftime
import numpy as np

import functions as f
import settings as s
import IO

class Run(object):
    """
    Prepare a protein for usage in spherical boundary conditions.
    """
    def __init__(self, prot, radius, center, water, cluster, noclean, forcefield, mutchain, *args, **kwargs):
        self.prot = prot
        self.radius = float(radius)
        self.center = center
        self.water = water
        self.noclean = noclean
        self.PDB = {}
        self.cluster = cluster
        self.original_charges = {}
        self.chains = []
        self.forcefield = forcefield
        self.mutchain = mutchain
        self.log = {'INPUT':        f'protprep.py -p {prot} -r {radius} -c {center} -w={water} -V={noclean}',
                    'CENTER':       None,
                    'DECHARGE':     [],
                    'CHARGE':       [],
                    'TOTAL_CHARGE': 0,
                    'CTERM':        [],
                    'CYX':          [],
                    'QRESN':        {},
                    'PDB2Q':        {},
                    'QRES_LIST':    [],
                    'TIME':         strftime("%Y-%m-%d %H:%M:%S", gmtime()),
                    'DIST':         [],}
    
    # Get the coordinates to use for the simulation sphere center
    def get_center_coordinates(self):
        # extract --center information
        center = self.center.strip('[]').split(':')

        # Open the PDB file and iterate over ATOM and HETATM lines
        with open(self.prot) as infile:
            for line in infile:
                if line.startswith(('ATOM', 'HETATOM')):
                    line = IO.pdb_parse_in(line)

                    id, atn, chain, resi = line[1], line[2], line[5], line[6]
                    xyz   = [float(line[8]), float(line[9]), float(line[10])]
                    
                    # If centered on a residue, find the xyz coordinates for the CB (or HA3 for Glycine)
                    if center[0] == 'RESN':
                        if (resi == int(center[1]) and atn in {'CB', 'HA3'} and chain == self.mutchain):
                            self.center = xyz

                    # If centered on a specific atom, find the xyz coordinates for the atom number
                    elif center[0] == 'ATN':
                        if id == int(center[1]):
                            self.center = xyz
                    
                    # If center coordinates are given directly, store those
                    elif len(center) == 3:
                        self.center = xyz
                        
                    else:
                        print('Could not get center')
                        
        self.log['CENTER'] = '{} {} {}'.format(*self.center)

    # Parse the PDB file to be readable for Q
    def prepwizard_parse(self):

        write = True

        # Open the PDB file and iterate over the ATOM and HETATM lines
        with open(self.prot) as infile, open(self.prot[:-4] + '_noH.pdb', 'w') as outfile:
            for line in infile:
                if line.startswith(('ATOM', 'HETATM')):
                
                    line = IO.pdb_parse_in(line)

                    atn, resn, chain, resi = line[2], line[4], line[5], line[6]

                    # Don't include non-standard terminal residue hydrogens 
                    if (atn in ('H1', 'H2', 'H3')) and (resn not in ('HOH', 'SOL')): continue

                    # rename GROMACS ILE CD to CD1
                    if resn == 'ILE' and atn == 'CD':
                        atn = 'CD1'

                    # Change residue name of GROMACS waters
                    if resn == 'SOL' or resn == 'HOH':
                        
                        # checks if water crystal waters are to be retained
                        if self.water == False: continue

                        resn = 'HOH'
                        if atn == 'OW' or atn == 'O':
                            atn = 'O'
                            
                            # Determines if water is within sphere radius (write=True) or not (write=False)
                            coord1 = self.center
                            coord2 = [float(line[8]), float(line[9]), float(line[10])]
                            write = f.euclidian_overlap(coord1, coord2, self.radius)
                    
                    # Get the charges from the hydrogen connections
                    if resn in IO.charged_res:
                        if atn in IO.charged_res[resn]:
                            if resn not in self.original_charges:
                                self.original_charges[chain] = {}
                            
                            if atn not in IO.charged_res[resn]:
                                self.original_charges[chain][resi] = 'HIP'
                            else:
                                self.original_charges[chain][resi] = IO.charged_res[resn][atn]

                    # Correction of capped termini by maestro
                    if resn == 'NMA' and atn[0] != 'H':
                        outline = IO.pdb_parse_out(line)
                        outline = f.correct_CT(outline)
                        outfile.write(outline  + '\n')
                    elif resn == 'ACE' and atn[0] != 'H':
                        outline = IO.pdb_parse_out(line)
                        outline = f.correct_NT(outline)
                        outfile.write(outline  + '\n')
                    else:
                        outline = IO.pdb_parse_out(line)

                    # Write this atom if all conditions are met
                    if write == True:
                        outfile.write(outline + '\n')
                
                # Simply write out GAP and TER lines
                elif line.startswith(('GAP', 'TER')):
                    outfile.write(line + '\n')
                else:
                    continue
    
    # Read the parsed PDB file to create the self.PDB dict used for applying modifications
    def readpdb(self):
        i = 0
        RES_ref = None
        pdbfile = self.prot[:-4] + '_noH.pdb'
        if self.water == True:
            self.PDB['w'] = {}

        # Open th parsed PDB file and iterate over the ATOM and HETATM lines
        with open(pdbfile) as infile:
            for line in infile:
                if line.startswith(('ATOM', 'HETATM')):
                    line = IO.pdb_parse_in(line)
                    
                    id, resn, chain, resi = line[1], line[4], line[5], line[6]
                    # Store the residue number offset
                    if not RES_ref:
                        RES_ref = resi - 1

                    # Store the chains present in the PDB file
                    if not chain in self.chains:
                        self.chains.append(chain)
                        self.PDB[chain] = {}

                    # Determine the charged state of ionziable residues
                    if chain in self.original_charges:
                        if resi in self.original_charges[chain]:
                            resn = self.original_charges[chain][resi]
                    
                    # Add the atom to the self.PDB dict by chain key, and only keep crystallographic waters if indicated
                    if self.water == True:
                        if resn == 'HOH':
                            self.PDB['w'][id] = line
                        else:
                            self.PDB[chain][id] = line
                    else:
                        if resn != 'HOH':
                            self.PDB[chain][id] = line

                    # Update the logger (writes out to protPREP.log)
                    if resi != RES_ref:
                        RES_ref = resi

                        if resn == 'HOH' or resn == 'SOL':
                            continue
                        else:
                            # check chain ID
                            if chain not in self.log['PDB2Q']:
                                self.log['PDB2Q'][chain] = {}
                            
                            if chain not in self.log['QRESN']:
                                self.log['QRESN'][chain] = {}
                                
                            i += 1
                            self.log['PDB2Q'][chain][i] = resi
                            self.log['QRESN'][chain][resi] = i
        
    def decharge(self):
        charged_res = {'GLU':['GLH', 'CD', -1,  'HD1'], 
                       'ASP':['ASH', 'CG', -1,  'HE1'], 
                       'ARG':['ARN', 'CZ',  1, 'HH12'],  
                       'LYS':['LYN', 'NZ',  1,  'HZ3'],
                       'HIP':['HID', 'CG',  1,  'HE2'],}
                       
        decharge = {}
        
        coord1 = self.center
        shell = float(self.radius) - 3.0    # Shell radius for decharging residues in boundary        

        for chain in [chains for chains in self.chains if chains in self.PDB]:
            for _, at in self.PDB[chain].items():

                atn, resn, resi = at[2], at[4], at[6]
                xyz  = [float(at[8]), float(at[9]), float(at[10])]

                if resn in charged_res:
                    if atn.strip() == charged_res[resn][1]:
                        coord2 = xyz
                        dist = float(self.radius) - np.sqrt(sum((c2 - c1)**2 for c1, c2 in zip(coord1, coord2)))
                        self.log['DIST'].append(f'{resn:<10}{resi:<10}{dist:.2f}')
                        if f.euclidian_overlap(coord1, coord2, shell) == False:
                            if not chain in decharge:
                                decharge[chain] = [resi]
                            else:
                                decharge[chain].append(resi)
                            self.log['DECHARGE'].append(f'{self.log["QRESN"][chain][resi]:<10}{resi:<10}{chain:<10}{resn:<10}')
                            
        # Check if the decharged residue is part of a salt bridge and
        for chain in [chains for chains in self.chains if chains in self.PDB]:
            for _, at in self.PDB[chain].items():
                
                atn, resn, resi = at[2], at[4], at[6]
                xyz  = [float(at[8]), float(at[9]), float(at[10])]

                if chain not in decharge:
                    continue
                if resi in decharge[chain] and resn in charged_res:
                    if atn == charged_res[resn]:
                        coord1 = xyz
                        for chain2 in self.PDB:
                            for _, at2 in self.PDB[chain2].items():

                                atn2, resn2, resi2 = at2[2], at2[4], at2[6]
                                xyz2  = [float(at[8]), float(at[9]), float(at[10])]

                                if (at == at2) or (resi2 in decharge[chain]) or (resi2 in decharge[chain2]):
                                    continue

                                if resn2 in charged_res:
                                    if atn2 == charged_res[resn2][1]:
                                        coord2 = xyz2

                                        if f.euclidian_overlap(coord1, coord2, 4.0) == True:
                                            decharge[chain2].append(resi2)
                                            self.log['DECHARGE'].append(f'{self.log["QRESN"][chain2][resi2]:<10}{resi2:<10}{chain2:<10}{resn2:<10}')

        protons = {}
        # Neutralize designated residues and determine the total protein charge in the sphere
        for chain in [chains for chains in self.chains if chains in self.PDB]:
            for _, at in self.PDB[chain].items():

                atn, resn, resi = at[2], at[4], at[6]

                if (chain in decharge) and (resi in decharge[chain]) and (resn in charged_res):
                    if (charged_res[resn][2] > 0) and (atn == charged_res[resn][3]):
                            if chain not in protons:
                                protons[chain] = []
                            protons[chain].append(_)

                    resn = charged_res[resn][0]
                    self.PDB[chain][_][4] = resn
                   
                if resn in charged_res:
                    if atn == charged_res[resn][1]:
                        self.log['CHARGE'].append(f'{self.log["QRESN"][chain][resi]:<10}{resi:<10}{chain:<10}{resn:<10}')
                        self.log['TOTAL_CHARGE'] += charged_res[resn][2]

        [self.PDB[chain].pop(atom) for chain, atoms in protons.items() for atom in atoms]
                    

    def set_OXT(self):
        decharged = ['LYN','ARN','GLH','ASH']
        c_term = []
        remove = []

        for chain in [chains for chains in self.chains if chains in self.PDB]:
            for _, at in self.PDB[chain].items():

                atn, resn, resi = at[2], at[4], at[6]

                if atn == 'OXT' or atn == 'O1':
                    c_term.append((resi, chain))
                    self.log['CTERM'].append(f'{resi} {resn}')

        for cterm in c_term:
            for chain in [chains for chains in self.chains if chains in self.PDB]:
                for _, at in self.PDB[chain].items():

                    atn, resn, resi = at[2], at[4], at[6]
                
                    if resn in decharged:
                        if atn == 'O1':
                            self.PDB[chain][_][2] = 'O'

                        if atn == 'O2':
                            remove.append([chain, _])
                                
                    else:
                        if resi == cterm[0] and chain == cterm[1]:
                            self.PDB[chain][_][4] = f'C{resn}'

                        if atn == 'O1':
                            self.PDB[chain][_][2] = 'O'

                        if atn == 'O2':
                            self.PDB[chain][_][2] = 'OXT'
    
        # Remove atoms from c or n terminal decharged residues
        for at in remove:
            del self.PDB[at[0]][at[1]]
                
    def get_CYX(self):
        cys = []
        cyx = []
        cys_bond = 2.2
        cysbond_list = []
        i = 0
        
        # Reduce coordinate array
        for chain in [chains for chains in self.chains if chains in self.PDB]:
            for _, at in self.PDB[chain].items():

                atn, resn, resi = at[2], at[4], at[6]
                xyz  = (float(at[8]), float(at[9]), float(at[10]))

                Q_resi = self.log['QRESN'][chain][resi]
                if (resn == 'CYS' or resn == 'CYX') and atn == 'SG':
                    cys.append([resi, xyz, Q_resi])

        # Construct S-S bond matrix
        for i, sg1 in enumerate(cys):
            for sg2 in cys[i+1:]:
                    overlap = f.euclidian_overlap(sg1[1], sg2[1], cys_bond)
                    if overlap:
                        self.log['CYX'].append(f'{sg1[2]} {sg2[2]}')
                        print(f'Overlap between CYS{sg1[0]} and CYS{sg2[0]}')
                        cysbond_list.append(f'{sg1[0]}-{sg2[0]}')

        cysbond_list = list(set(cysbond_list))
        for pair in cysbond_list:
            l1, l2 = map(int, pair.split('-'))
            if l1 < l2:
                cyx.append(l1)
                cyx.append(l2)
        for chain in [chains for chains in self.chains if chains in self.PDB]:
            for _, at in self.PDB[chain].items():

                atn, resn, resi = at[2], at[4], at[6]

                if resi in cyx and resn == 'CYS':
                    self.PDB[chain][_][4] = 'CYX'
                
    def write_tmpPDB(self):
        with open(self.prot[:-4] + '_tmp.pdb', 'w') as outfile:
            for chain in [chains for chains in (self.chains + ['w']) if chains in self.PDB]:
                for _, at in self.PDB[chain].items():
                    outline = IO.pdb_parse_out(at) + '\n'
                    outfile.write(outline)

                if chain != list(self.PDB.keys())[-1] or chain == 'w':
                    outfile.write('TER\n')

        
    def write_qprep(self):
        replacements = {'FF_LIB'    :   s.FF_DIR + '/OPLSAAM.lib',
                        'FF_PRM'    :   s.FF_DIR + '/OPLSAAM.prm',
                        'PROTPDB'   :   self.prot[:-4] + '_tmp.pdb',
                        'CENTER'    :   self.log['CENTER'],
                        'SPHERE'    :   '{:.1f}'.format(self.radius),
                        'SOLVENT'   :   '1 HOH',
                        'FF_LIB'    :   f'{s.FF_DIR}/{self.forcefield}.lib',
                        'FF_PRM'    :   f'{s.FF_DIR}/{self.forcefield}.prm',
                       }

        with open (s.INPUT_DIR + '/qprep_protprep.inp') as infile, open ('qprep.inp', 'w') as outfile:
            for line in infile:
                line = IO.replace(line, replacements)
                outfile.write(line)
                if line[0:8] == '!addbond':
                    for pair in self.log['CYX']:
                        pair = pair.split()
                        cyx_line = f'addbond {pair[0]}:SG {pair[1]}:SG y\n'
                        outfile.write(cyx_line)            
        
    def run_qprep(self):
        hpc_commands = getattr(s, self.cluster)
        qprep = hpc_commands['QPREP']
        os.system(f'{qprep} < qprep.inp > qprep.out')

        error = os.popen("grep 'ERROR' qprep.out").read()
        if error:
            sys.exit('>>>>> ERROR in qprep.out')
        
    def write_pdb_out(self):
        os.system("grep -v HOH top_p.pdb > protein.pdb")
        os.system("grep HOH top_p.pdb > water.pdb")
                
    def write_log(self):
        with open('protPREP.log', 'w') as outfile:
            outfile.write('--------------------------------------------------------------------\n')
            outfile.write('This file contains some output information regarding the preperation\n')
            outfile.write('process of the protein for spherical boundary conditions in Q.\n\n')
            outfile.write('The command line input was:\n')
            outfile.write(f"{self.log['INPUT']}\n\n")
            outfile.write(f"{'Date:':47}{self.log['TIME']:>21}\n")
            outfile.write(f"{'Inputfile:':47}{self.prot:>21}\n")
            outfile.write(f"{'Sphere center:':41}{self.center[0]:>9.3f}{self.center[1]:>9.3f}{self.center[2]:>9.3f}\n")
            outfile.write(f"{'Sphere radius:':47}{self.radius:>21.1f}\n")
            outfile.write(f"{'Total charge in sphere:':47}{self.log['TOTAL_CHARGE']:>21}\n")
            outfile.write('--------------------------------------------------------------------\n')
            outfile.write(f"{'Q_RESN':10}{'PDB_IN':10}{'BOUNDARY_DIST':10}\n")
            for line in self.log['DIST']:
                outfile.write(line + '\n')  
            outfile.write('The following residues have been decharged:\n')
            outfile.write(f"{'Q_RESN':10}{'PDB_IN':10}{'CHAIN':10}{'RESNAME':10}\n")
            for line in self.log['DECHARGE']:
                outfile.write(line + '\n')
            outfile.write('--------------------------------------------------------------------')
            outfile.write('\nThe following charged residues are in the sphere:\n')
            outfile.write(f"{'Q_RESN':10}{'PDB_IN':10}{'CHAIN':10}{'RESNAME':10}\n")
            for line in self.log['CHARGE']:
                outfile.write(line + '\n')
            outfile.write('--------------------------------------------------------------------')  
            outfile.write('\n')
            outfile.write('The following S-S bonds have been found:\n')
            outfile.write(f"{'Q_CYS1':10}{'Q_CYS2':10}\n")
            for line in self.log['CYX']:
                cys1, cys2 = line.split()
                outfile.write(f"{cys1:10}{cys2:10}\n")
            outfile.write('--------------------------------------------------------------------')  
            outfile.write('\n')
            outfile.write('The following is a mapping of residue numbers in Q and the input \n')
            outfile.write(f'pdbfile "{self.prot}":\n')
            outfile.write(f"{'Q_RESN':10}{'PDB_IN':10}{'CHAIN':10}{'RESNAME':10}\n")
            for chain in [chains for chains in self.chains if chains in self.PDB]:
                for _, at in self.PDB[chain].items():

                    atn  = at[2]
                    resn = at[4]
                    resi = at[6]
                    
                    if atn == 'CA':
                        Q_resi = self.log['QRESN'][chain][resi]
                        outfile.write(f"{Q_resi:<10}{resi:<10}{chain:10}{resn:10}\n")
                    
    def cleanup(self):
        if self.noclean == False:
            os.system(f'rm {self.prot[:-4]}_tmp.pdb {self.prot[:-4]}_noH.pdb qprep.inp qprep.out top_p.pdb')
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='protprep',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '===== Prepare a protein file with qprep for further processing; required for QligFEP and QresFEP =====')
    
    parser.add_argument('-p', '--prot',
                        dest = "prot",
                        required = True,
                        help = "PDB file of the protein")

    parser.add_argument('-mc', '--mutchain',
                        dest = "mutchain",
                        required=False,
                        default = False,
                        help = "Use to specficy the PDB chain if intended center is residue number")

    parser.add_argument('-w', '--nowater',
                        dest = "water",
                        required = False,
                        default = True,
                        action = 'store_false',
                        help = "Turn on if crystal waters ought to be removed")

    parser.add_argument('-f', '--forcefield',
                        dest = "forcefield",
                        default = 'OPLSAAM',
                        choices = ['OPLSAAM', 'OPLS2015', 'OPLS2005', 'AMBER14sb', 'CHARMM36'],
                        help = "Use to specficy forcefield")
    
    parser.add_argument('-r', '--sphereradius',
                        dest = "sphereradius",
                        required = True,
                        help = "Radius of intended sphere")    
    
    parser.add_argument('-c', '--spherecenter',
                        dest = "spherecenter",
                        required = True,
                        help = "Center of intended sphere; Either residue number (RESN:$), atomnumber (ATN:$) or explicit coordinates (X:Y:Z)")
    
    parser.add_argument('--noclean',
                        dest = "noclean",
                        default = False,
                        action = 'store_true',
                        help = "Turn on if qprep input files ought not be deleted")
    
    parser.add_argument('-P', '--preplocation',
                        dest = "preplocation",
                        default = s.DEFAULT,
                        help = "Use to specify location of Q executables protprep ought to use")
    
    args = parser.parse_args()
    run = Run(prot       = args.prot,
              radius     = args.sphereradius,
              center     = args.spherecenter,
              water      = args.water,
              noclean    = args.noclean,
              cluster    = args.preplocation,
              forcefield = args.forcefield,
              mutchain   = args.mutchain)
    
    run.get_center_coordinates()        # 00
    run.prepwizard_parse()              # 01
    run.readpdb()                       # 02
    run.decharge()                      # 03
    run.set_OXT()                       # 04
    run.get_CYX()                       # 05
    run.write_tmpPDB()                  # 06
    run.write_qprep()                   # 07
    run.run_qprep()                     # 08
    run.write_pdb_out()                 # 09
    run.write_log()                     # 10
    run.cleanup()                       # 11
