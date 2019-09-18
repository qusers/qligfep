import argparse
from subprocess import check_output
import os
import collections
from time import gmtime, strftime
import numpy as np

import functions as f
import settings as s
import IO

class Run(object):
    """
    Prepare a protein for usage in spherical boundary conditions.
    """
    def __init__(self, prot, sphereradius, spherecenter, include, water, 
                 preplocation, origin, noclean, *args, **kwargs):
        self.prot = prot
        self.radius = float(sphereradius)
        self.center = spherecenter
        self.include = include
        self.water = water
        self.noclean = noclean
        self.PDB = {}
        self.preplocation = preplocation
        self.origin = origin
        self.original_charges = {}
        self.chains = []
            
        self.log = {'INPUT':'protPREP.py -p {} -r {} -c {} -w={} -V={}'.format(prot, 
                                                                               sphereradius, 
                                                                               spherecenter, 
                                                                               water,
                                                                               noclean
                                                                              ),
                    'CENTER':None,
                    'DECHARGE':[],
                    'CHARGE':[],
                    'TOTAL_CHARGE':0,
                    'CTERM':[],
                    'CYX':[],
                    'QRESN':{},
                    'PDB2Q':{},
                    'QRES_LIST':[],
                    'TIME':strftime("%Y-%m-%d %H:%M:%S", gmtime()),
                   }
    
    def get_center_coordinates(self):
        center = self.center
        center = center.strip('[')
        center = center.strip(']')
        center = center.split(':')
        with open(self.prot) as infile:
            for line in infile:
                if line.startswith(self.include):
                    line = IO.pdb_parse_in(line)
                    if center[0] == 'RESN':
                        if line[6] == int(center[1]) \
                        and line[2].strip() == 'CA':
                            self.center = [float(line[8]),
                                           float(line[9]), 
                                           float(line[10])
                                          ]

                    elif center[0] == 'ATN':
                        if line[1] == int(center[1]):
                            self.center = [float(line[8]),
                                           float(line[9]),
                                           float(line[10])
                                          ]
                    
                    elif len(center) == 3:
                        self.center = [float(center[0]),
                                       float(center[1]),
                                       float(center[2])
                                      ]
                        
                    else:
                        print('Could not get center')
                        
        self.log['CENTER'] = '{} {} {}'.format(*self.center)
    
    # The name of this function should be more general when the various inputs are there
    def prepwizard_parse(self):
        if self.origin == 'gromacs':
            with open(self.prot) as infile, \
                 open(self.prot[:-4] + '_noH.pdb', 'w') as outfile:
                for line in infile:
                    if line.startswith(self.include):
                        if line[13] != 'H':
                            line = IO.pdb_parse_in(line)
                            # Change residue name of waters
                            if line[4] == 'SOL':
                                line[4] = 'HOH'
                                line[2] = 'O'
                                if self.water != False:
                                    line_out = IO.pdb_parse_out(line)
                                    
                                else:
                                    continue
                                
                            elif line[4] == 'ILE' and line[2] == 'CD':
                                line[2] = 'CD1'
                                line_out = IO.pdb_parse_out(line)
                                
                            elif line[4] == 'CL-':
                                continue
                                
                            line_out = IO.pdb_parse_out(line)
                            outfile.write(line_out + '\n')

                        else:
                            line = IO.pdb_parse_in(line)
                        
                        # Get the charges from the hydrogen connections
                        # NOTE: this might be more common and thus less lines of code might be
                        # needed, check when implementing MolProbity!!
                        if line[4] in IO.charged_res:
                            if line[2] in IO.charged_res[line[4]]:
                                if line[6] in self.original_charges:
                                    self.original_charges[line[6]] = ['HIP',line[5]]
                                else:
                                    self.original_charges[line[6]] = [IO.charged_res[line[4]][line[2]],
                                                                     line[5]
                                                                     ]
                                    
                            
        elif self.origin == 'maestro':
            with open(self.prot) as infile, \
                 open(self.prot[:-4] + '_noH.pdb', 'w') as outfile:
                for line in infile:
                    if line.startswith(self.include):
                        line = IO.pdb_parse_in(line)
                        if line[2][0] != 'H':
                            outline = IO.pdb_parse_out(line)
                            outfile.write(outline  + '\n')
                        
                        # Get the charges from the hydrogen connections
                        if line[4] in IO.charged_res:
                            if line[2] in IO.charged_res[line[4]]:
                                if line[6] in self.original_charges:
                                    self.original_charges[line[6]] = ['HIP',line[5]]
                                else:
                                    self.original_charges[line[6]] = [IO.charged_res[line[4]][line[2]],
                                                                      line[5]
                                                                     ]

    
    def readpdb(self):
        i = 0
        if self.origin == 'gromacs':
            pdbfile = self.prot[:-4] + '_noH.pdb'
            
        elif self.origin == 'maestro':
            pdbfile = self.prot[:-4] + '_noH.pdb'
            
        with open(pdbfile) as infile:
            for line in infile:
                if line.startswith(self.include):
                    header = IO.pdb_parse_in(line)
                    RES_ref = header[6] - 1
                    break
                    
            for line in infile:
                line = IO.pdb_parse_in(line)
                if not line[5] in self.chains:
                    self.chains.append(line[5])
                # construct chain based container
                self.PDB[line[5]] = {}

        with open(pdbfile) as infile:
            for line in infile:
                if line.startswith(self.include):
                    line = IO.pdb_parse_in(line)
                    RES = line[6]
                    if line[6] in self.original_charges:
                        if self.original_charges[line[6]][1] == line[5]:
                            line[4] = self.original_charges[line[6]][0]
                    
                    if self.water == True:
                        self.PDB[line[5]][line[1]] = line
                        
                    elif self.water == False:
                        if line[4] != 'HOH':
                            self.PDB[line[5]][line[1]] = line

                    if RES != RES_ref:
                        RES_ref = RES
                        i += 1
                        self.log['QRESN'][line[6]] = i
                        self.log['PDB2Q'][i] = [line[6]]
                        if line[4].strip() != 'HOH':
                            self.log['QRES_LIST'].append('{:<10d}{:<10d}{:<10}'.format(i, 
                                                                                  int([line[6]][0]),
                                                                                  line[4]
                                                                                 ))
    def decharge(self):
        charged_res = {'GLU':['GLH', 'CD', -1], 
                       'ASP':['ASH', 'CG', -1], 
                       'ARG':['ARN', 'CZ', 1], 
                       'LYS':['LYN', 'NZ', 1], 
                       'HIP':['HID', 'CG', 1]
                      }
        
        coord1 = self.center
        decharge = {}
        # Distance for decharging residues in boundary
        rest_bound = float(self.radius) - 3.0
        for chain in self.PDB:
            for key in self.PDB[chain]:
                at = self.PDB[chain][key]
                if at[4] in charged_res:
                    if at[2].strip() == charged_res[at[4]][1]:
                        coord2 = [float(at[8]), 
                                  float(at[9]), 
                                  float(at[10])
                                 ]
                        if f.euclidian_overlap(coord1, coord2, rest_bound) == False:
                            if not at[5] in decharge:
                                decharge[at[5]] = [at[6]]
                            else:
                                decharge[at[5]].append(at[6])
                            outline = '{:<10}{:<10}{:<10}'.format(self.log['QRESN'][at[6]],
                                                                                    at[6], 
                                                                                    at[4])
                            self.log['DECHARGE'].append(outline)
        # Check if the decharged residue is part of a salt bridge and
        # neutralize this residue as well
        for chain in self.PDB:
            for key in self.PDB[chain]:
                at = self.PDB[chain][key]
                if chain not in decharge:
                    continue
                if at[6] in decharge[chain] and at[2].strip() == charged_res[at[4]][1]:
                    coord1 = [float(at[8]), 
                              float(at[9]), 
                              float(at[10])
                             ]
                    for chain2 in self.PDB:
                        for key2 in self.PDB[chain2]:
                            at_2 = self.PDB[chain2][key2]
                            if at_2[4] in charged_res:
                                if at_2[2].strip() == charged_res[at_2[4]][1]:
                                    coord2 = [float(at_2[8]), 
                                              float(at_2[9]), 
                                              float(at_2[10])
                                             ]
                                    if at != at_2 and at_2[6] not in decharge:
                                        if f.euclidian_overlap(coord1, coord2, 4.0) == True:
                                            decharge[chain].append(at_2[6])
                                            outline = '{:<10}{:<10}{:<10}'.format(self.log['QRESN'][at_2[6]],
                                                                                                    at_2[6], 
                                                                                                    at_2[4])
                                            self.log['DECHARGE'].append(outline)
                                            
        # Get the charged residues in the sphere and the total charge of these residues in the sphere
        for chain in self.PDB:
            for key in self.PDB[chain]:
                at = self.PDB[chain][key]
                if at[6] in decharge:
                    at[4] = charged_res[at[4]][0]

                else:
                    if at[4] in charged_res:
                        if at[2].strip() == charged_res[at[4]][1]:
                            self.log['CHARGE'].append('{:<10}{:<10}{:<10}'.format(self.log['QRESN'][at[6]],
                                                                                  at[6], 
                                                                                  at[4]))
                            self.log['TOTAL_CHARGE'] += charged_res[at[4]][2]
            
    def set_OXT(self):
        ## NOTE WARNING, ETC: Q AA codes are sometimes 3, sometimes 4, they MUST be updated to pdb
        ## standards ASAP!!!
        CTERM = []
        for chain in self.PDB:
            for key in self.PDB[chain]:
                at = self.PDB[chain][key]
                if at[2].strip() == 'OXT':
                    CTERM.append(at[6])
                    self.log['CTERM'].append('{} {}'.format(at[6], at[4]))

                if self.origin == 'gromacs':
                    if at[2].strip() == 'O1':
                        CTERM.append(at[6])
                        self.log['CTERM'].append('{} {}'.format(at[6], at[4]))
            
            for key in self.PDB[chain]:
                at = self.PDB[chain][key]
                if at[6] in CTERM:
                    self.PDB[chain][key][4] = 'C' + at[4]

                if self.origin == 'gromacs':
                    if at[2] == 'O1':
                        self.PDB[chain][key][2] = 'O'

                    if at[2] == 'O2':
                        self.PDB[chain][key][2] = 'OXT'

                
    def get_CYX(self):
        cys = []
        cyx = []
        cys_bond = 2.2
        cys_mat = []
        i = 0
        k = -1
        
        # Reduce coordinate array
        for chain in self.PDB:
            for key in self.PDB[chain]:
                at = self.PDB[chain][key]
                if at[4] == 'CYS' or at[4] == 'CYX' and at[2].strip() == 'SG':
                    cys.append([at[6], (at[8], at[9], at[10])])

            # Construct S-S bond matrix
            for SG_1 in cys:
                cys_list = [SG_1[0]]
                for SG_2 in cys:
                    cys_list.append(f.euclidian_overlap(SG_1[1], SG_2[1], cys_bond))
                cys_mat.append(cys_list)

            # Fix to better handling
            try:
                total = len(cys_mat[0]) -1

                for line in cys_mat:
                    k += 1
                    for j in range(i, total):
                        if cys_mat[i][j+1] == True and cys_mat[k][0] != cys_mat[j][0]:
                            cyx.append(cys_mat[k][0])
                            cyx.append(cys_mat[j][0])
                            outline = '{:<10}{:<10}{:<10}{:<10}'.format(self.log['QRESN'][cys_mat[k][0]], 
                                                                        self.log['QRESN'][cys_mat[j][0]],
                                                                        cys_mat[k][0],
                                                                        cys_mat[j][0])
                            self.log['CYX'].append(outline)

                    i += 1
                for chain in self.PDB:
                    for key in self.PDB[chain]:
                        at = self.PDB[chain][key]
                        if at[6] in cyx:
                            self.PDB[chain][key][4] = 'CYX'

            except:
                return None
                
    def write_tmpPDB(self):
        with open(self.prot[:-4] + '_tmp.pdb', 'w') as outfile:
            for chain in self.PDB:
                atom_numbers = sorted(list(self.PDB[chain].keys()))
                for atom in atom_numbers:
                    outline = IO.pdb_parse_out(self.PDB[chain][atom]) + '\n'
                    outfile.write(outline)
            outfile.write('GAP\n')
        
    def write_qprep(self):
        replacements = {'FF_LIB'    :   s.FF_DIR + '/OPLS2015.lib',
                        'FF_PRM'    :   s.FF_DIR + '/OPLS2015.prm',
                        'PROTPDB'   :   self.prot[:-4] + '_tmp.pdb',
                        'CENTER'    :   self.log['CENTER'],
                        'SPHERE'    :   '{:.1f}'.format(self.radius),
                        'SOLVENT'   :   '1 HOH'
                       }
        
        with open (s.INPUT_DIR + '/qprep_protprep.inp') as infile, \
            open ('qprep.inp', 'w') as outfile:
            for line in infile:
                line = IO.replace(line, replacements)
                outfile.write(line)
                if line[0:8] == '!addbond':
                    for line in self.log['CYX']:
                        line = line.split()
                        outline = 'addbond {}:SG {}:SG y\n'.format(*line)
                        outfile.write(outline)            
        
    def run_qprep(self):
        qprep = s.Q_DIR[self.preplocation] + 'qprep'
        options = ' < qprep.inp > qprep.out'
        # Somehow Q is very annoying with this < > input style so had to implement
        # another function that just calls os.system instead of using the preferred
        # subprocess module....
        IO.run_command(qprep, options, string = True)
        
    ##### THIS PART COMES AFTER THE TEMP .pdb IS WRITTEN ######
    def write_pdb_out(self):
        waters ={'HOH': ['O', 'H1', 'H2'],
                 'SOL': ['OW1', 'HW1', 'HW2'] 
                }
        waters_tokeep = []
        
        with open('top_p.pdb') as infile:
            for line in infile:
                if line.startswith(self.include):
                    line = IO.pdb_parse_in(line)
                    if line[4].strip() in waters and \
                       line[2].strip() == waters[line[4].strip()][0]:
                        coord1 = self.center
                        coord2 = [float(line[8]), 
                                  float(line[9]), 
                                  float(line[10])
                                 ]
                        if f.euclidian_overlap(coord1, coord2, self.radius) == True:
                            waters_tokeep.append(line[6])
                            
        with open('top_p.pdb') as infile, \
             open('water.pdb', 'w') as watout, \
             open('protein.pdb', 'w') as protout:
                    
            for line in infile:
                if line.startswith(self.include):
                    line = IO.pdb_parse_in(line)
                    if line[6] in waters_tokeep:
                        outline = IO.pdb_parse_out(line) + '\n'
                        watout.write(outline)
                        
                    if line[4] not in waters:
                        outline = IO.pdb_parse_out(line) + '\n'
                        protout.write(outline)
                if line[0:3] == 'TER' or line == 'GAP':
                    protout.write(line)
                
    def write_log(self):
        with open('protPREP.log', 'w') as outfile:
            outfile.write('--------------------------------------------------------------------\n')
            outfile.write('This file contains some output information regarding the preperation\n')
            outfile.write('process of the protein for spherical boundary conditions in Q.\n\n')
            outfile.write('The command line input was:\n')
            outfile.write('{}\n\n'.format(self.log['INPUT']))
            outfile.write('{:47}{:>21}\n'.format('Date:', 
                                                  self.log['TIME']))
            outfile.write('{:47}{:>21}\n'.format('Inputfile:', 
                                                  self.prot))
            outfile.write('{:41}{:>9.3f}{:>9.3f}{:>9.3f}\n'.format('Sphere center:',
                                                                   self.center[0],
                                                                   self.center[1],
                                                                   self.center[2]))
            outfile.write('{:47}{:>21.1f}\n'.format('Sphere radius:',
                                                    self.radius))            
            outfile.write('{:47}{:>21}\n'.format('Total charge in sphere:', 
                                                  self.log['TOTAL_CHARGE']))
            outfile.write('--------------------------------------------------------------------')  
            outfile.write('\nThe following residues have been decharged:\n')
            outfile.write('{:10}{:10}{:10}\n'.format('Q_RESN', 'PDB_IN', 'RESNAME'))
            for line in self.log['DECHARGE']:
                outfile.write(line + '\n')
                
            outfile.write('--------------------------------------------------------------------')
            outfile.write('\nThe following charged residues are in the sphere:\n')
            outfile.write('{:10}{:10}{:10}\n'.format('Q_RESN', 'PDB_IN', 'RESNAME'))
            for line in self.log['CHARGE']:
                outfile.write(line + '\n')
            
            outfile.write('--------------------------------------------------------------------')  
            outfile.write('\n')
            outfile.write('The following S-S bonds have been found:\n')
            outfile.write('{:10}{:10}{:10}{:10}\n'.format('Q_CYS1', 'Q_CYS2', 'PDB_CYS1', 'PDB_CYS2'))
            for line in self.log['CYX']:
                outfile.write(line + '\n')
            
            outfile.write('--------------------------------------------------------------------')  
            outfile.write('\n')
            outfile.write('The following is a mapping of residue numbers in Q and the input \n')
            outfile.write('pdbfile "{}":\n'.format(self.prot))
            outfile.write('{:10}{:10}{:10}\n'.format('Q_RESN', 'PDB_IN', 'RESNAME'))
            for line in self.log['QRES_LIST']:
                outfile.write(line + '\n')
                    
    def cleanup(self):
        if self.noclean == False:
            os.remove(self.prot[:-4] + '_tmp.pdb')
            os.remove(self.prot[:-4] + '_noH.pdb')
            os.remove('qprep.inp')
            os.remove('qprep.out')
            os.remove('top_p.pdb')
            os.remove('complexnotexcluded.pdb')
            os.remove('tmp.top')
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='protPREP',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Generate FEP files for dual topology ligFEP == ')

    
    parser.add_argument('-p', '--prot',
                        dest = "prot",
                        required = True,
                        help = "protein pdb file")
    
    parser.add_argument('-r', '--sphereradius',
                        dest = "sphereradius",
                        required = True,
                        help = "radius of the sphere")    
    
    parser.add_argument('-w', '--nowater',
                        dest = "water",
                        required = False,
                        default = True,
                        action = 'store_false',
                        help = "set no water if crystal waters are NOT to be retained")
    
    parser.add_argument('-c', '--spherecenter',
                        dest = "spherecenter",
                        required = True,
                        help = "center of the sphere, can be residue number (RESN:$)," \
                        "atomnumber (ATN:$) or explicit coordinates (X:Y:Z)")
    
    parser.add_argument('--noclean',
                        dest = "noclean",
                        default = False,
                        action = 'store_true',
                        help = "If turned on, intermediate Q files will not be deleted")
    
    parser.add_argument('-P', '--preplocation',
                        dest = "preplocation",
                        default = 'LOCAL',
                        help = "define this variable if you are setting up your system elsewhere")
    
    parser.add_argument('-O', '--origin',
                        dest = "origin",
                        default = 'maestro',
                        choices = ['maestro', 'gromacs'],
                        help = "Use this flag to specficy with which program the .pdb file was written")
    
    args = parser.parse_args()
    run = Run(prot = args.prot,
              sphereradius = args.sphereradius,
              spherecenter = args.spherecenter,
              water = args.water,
              noclean = args.noclean,
              preplocation = args.preplocation,
              origin = args.origin,
              include = ('ATOM','HETATM')
             )
    
    run.prepwizard_parse()              # 00
    run.readpdb()                       # 01
    run.get_center_coordinates()        # 02
    run.decharge()                      # 03
    run.set_OXT()                       # 04
    run.get_CYX()                       # 05
    run.write_tmpPDB()                  # 06
    run.write_qprep()                   # 07
    run.run_qprep()                     # 08
    run.write_pdb_out()                 # 09
    run.write_log()                     # 10
    run.cleanup()                       # 11
