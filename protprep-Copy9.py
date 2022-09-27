import argparse
from subprocess import check_output
import os
import collections
from time import gmtime, strftime
import numpy as np

import functions as f
import settings as s
import IO
from calendar import c
import re
import shlex
import stat
from xml.etree.ElementInclude import include
from itertools import chain
from collections import defaultdict
from itertools import groupby

distance= 20
coord1 = 36.032,58.593,48.734
to_keep = {}
end_list=[]
final={}
j =1
r3={}
i3=0
i4=0
r4={}
f1=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
f2=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z']

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
        print(center)
        center = center.strip('[')
        center = center.strip(']')
        center = center.split(':')
        print(center)
        with open(self.prot) as infile:
            for line in infile:
                if line.startswith(self.include):
                    line = IO.pdb_parse_in(line)
#                    print(center[0])
                    if center[0] == 'RESN':
                        if line[5] == center[2]:
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
        
        
    def spheremaker(self):
        # Write the new pdb file with the new ordering
        with open(self.prot) as infile, \
             open('try007.pdb', 'w') as outfile2:
            for line2 in infile:
                line2 = IO.pdb_parse_in(line2)
                if type(line2) == list:
                    resnum = line2[6]
                    resname = line2[4]
                    dict_key=()
        #            dict_key=(resname,)#+str(resnum)
        #            resnum=resnum.strip()
        #            dict_key=dict_key+(resnum,)
        #            print(final[('TRP', '59')])
        #            print(final[dict_key])
        #            line2[6] = final[dict_key]
                    outline = IO.pdb_parse_out(line2)
        #            print(outline)
        #            print(line4)

                    outfile2.write(outline + '\n')
        #            print(line4[6])


        with open(self.prot) as infile, \
             open('try52.pdb', 'w') as outfile3:
             for line3 in infile:
                line3 = IO.pdb_parse_in(line3)


                if type(line3) == list:
                    coord2 = tuple(line3[8:11])
                    chain_id = line3[5]
                    resnum = line3[6]
                    combo = line3[5:7]

        #Include only the atoms that are inside the sphere in a new pdb file
                    if f.euclidian_overlap(coord1, coord2, distance) == True:
                        if chain_id not in to_keep:
                            to_keep[chain_id]=[]
                        to_keep[chain_id].append(resnum)

                    if chain_id not in to_keep:
                        to_keep[chain_id]=[]

                    if resnum in to_keep[chain_id]:
                        end_list.append(line3)
        #                print(line3)

        #Exit the loop so as to not to cut the residues that are in the sphere borders.Write the new pdb file.
        with open('try007.pdb') as infile, \
             open('try53.pdb', 'w') as outfile4:
        #     print(to_keep)
        #        print(infile)
             for line4 in infile:
        #        print(line2)
                line4 = IO.pdb_parse_in(line4)
        #        print(line4)
                if type(line4) == list:
                    coord2 = tuple(line4[8:11])
                    chain_id = line4[5]
                    resnum = line4[6]
                    combo = line4[5:7]
                if chain_id in to_keep:
                    if resnum in to_keep[chain_id]:
        #                print(line2)
                        outline2 = IO.pdb_parse_out(line4)

                        outfile4.write(outline2 + '\n')
        #                print(outline2)
        #            print(outline)        

        

    def prepwizard_parse(self):
        if self.origin == 'gromacs':
            print('work')
            with open('try53.pdb') as infile, \
                 open('try67.pdb'[:-4] + '_noH.pdb', 'w') as outfile:
                for line in infile:
                    tmp = line
                    if line.startswith(self.include) == False:
                        continue
                    
                    line = IO.pdb_parse_in(line)
                    if line[13].strip() == 'H':
                        if line[4] == 'SOL'
                            write = True
                        if line[4] == 'POPC':
                            write = True
                            line[4] = 'POP'
                        else:
                            write = False
                        
                    else:
                        write = True

                    # Change residue name of waters
                    if line[4] == 'T3P':
                        line [4]='HOH'
                    if line[4] == 'POPC':
                        line[4] = 'HOH'
                    if line[4] == 'SOL':
                        line[4] = 'HOH'
                        if line[2] == 'OW':
                            line[2] = 'O'
                            coord1 = self.center
                            coord2 = [float(line[8]), 
                                      float(line[9]), 
                                      float(line[10])
                                     ]
                            write = f.euclidian_overlap(coord1, coord2, self.radius + 5)
                            
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
                    
                    # Get the charges from the hydrogen connections
                    # NOTE: this might be more common and thus less lines of code might be
                    # needed, check when implementing MolProbity!!
                    if line[4] in IO.charged_res:
                        if line[2] in IO.charged_res[line[4]]:
                            if line[5] not in self.original_charges:
                                self.original_charges[line[5]] = {}
                            
                            if line[2] not in IO.charged_res[line[4]]:
                                self.original_charges[line[5]][line[6]] = 'HIP'
                                
                            else:
                                self.original_charges[line[5]][line[6]] = \
                                IO.charged_res[line[4]][line[2]]
                
                    if write == True:
                        outfile.write(line_out + '\n')
            
        elif self.origin == 'maestro':
#            print('work')
            with open('try53.pdb') as infile, \
                 open('try67.pdb', 'w') as outfile:
                for line in infile:
                    line = IO.pdb_parse_in(line)
#                    print(line[4])
                    if line[4] == 'T3P':
                        line [4]='HOH'
                    if line[4] == 'POPC':
                        line[4] = 'HOH'
                    if line[4] == 'SOL':
                        line[4] = 'HOH'
                        outline = IO.pdb_parse_out(line)

                        
                    outfile.write(IO.pdb_parse_out(line) + '\n')
            
            
            with open('try67.pdb') as infile, \
                 open('try67.pdb'[:-4] + '_noH.pdb', 'w') as outfile:
                for line in infile:
                    if line.startswith(self.include):
                        line = IO.pdb_parse_in(line)
                        if line[2][0] != 'H':
                            outline = IO.pdb_parse_out(line)
    #                            print(outline)
                            outfile.write(outline  + '\n')

                            # Get the charges from the hydrogen connections
                            if line[4] in IO.charged_res:
                                if line[2] in IO.charged_res[line[4]]:
                                    if line[5] not in self.original_charges:
                                        self.original_charges[line[5]] = {}

                                    if line[2] not in IO.charged_res[line[4]]:
                                        self.original_charges[line[5]][line[6]] = 'HIP'

                                    else:
                                        self.original_charges[line[5]][line[6]] = IO.charged_res[line[4]][line[2]]


    def readpdb(self):
        i = 0
        chain_ref = None
        if self.origin == 'gromacs':
            pdbfile = 'try67.pdb'[:-4] + '_noH.pdb'
            
        elif self.origin == 'maestro':
            pdbfile = 'try67.pdb'[:-4] + '_noH.pdb'
            
        with open(pdbfile) as infile:
            for line in infile:
                if line.startswith('ATOM'):
                    header = IO.pdb_parse_in(line)
#                    print(header)
                    print(header[6])
                    RES_ref = int(header[6]) - 1
                    break
                    
            for line in infile:
                line = IO.pdb_parse_in(line)
                if not line[5] in self.chains:
                    self.chains.append(line[5])
                # construct chain based container
                self.PDB[line[5]] = {}
                
        with open(pdbfile) as infile:
            if self.water == True:
                self.PDB['w'] = {}
            for line in infile:
                if line.startswith('ATOM'):
                    line = IO.pdb_parse_in(line)
                    chain = line[5]
                    RES = line[6]
                    if chain in self.original_charges:
                        if line[6] in self.original_charges[chain]:
                            line[4] = self.original_charges[chain][line[6]]
                    
                    if self.water == True:
                        if line[4] == 'HOH':
                            self.PDB['w'][line[1]] = line
                        else:
                            self.PDB[line[5]][line[1]] = line
                        
                    elif self.water == False:
                        if line[4] != 'HOH':
                            self.PDB[line[5]][line[1]] = line

                    if RES != RES_ref:
                        RES_ref = RES
                        #if line[4].strip() == 'HOH':
                        #    continue
                        if line[4].strip() == 'SOL':
                            continue
                            
                        else:
                            # check chain ID
                            if line[5] not in self.log['PDB2Q']:
                                self.log['PDB2Q'][line[5]] = {}
                            
                            if line[5] not in self.log['QRESN']:
                                self.log['QRESN'][line[5]] = {}
                                
                            i += 1
                            self.log['PDB2Q'][line[5]][i] = line[6]
                            self.log['QRESN'][line[5]][line[6]] = i

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
                            outline = '{:<10}{:<10}{:<10}{:<10}'.format(
                                                                        self.log['QRESN'][chain][at[6]],
                                                                                                 at[6],
                                                                                                 chain,
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
                                    if at == at_2:
                                        continue
                                    if at_2[6] in decharge[chain]:
                                        continue
                                    
                                    if at_2[6] in decharge[chain2]:
                                        continue
                                        
                                    if f.euclidian_overlap(coord1, coord2, 4.0) == True:
                                        decharge[chain2].append(at_2[6])
                                        outline = '{:<10}{:<10}{:<10}{:<10}'.format(
                                                      self.log['QRESN'][chain2][at_2[6]],
                                                                                at_2[6],
                                                                                chain2,
                                                                                at_2[4])
                                        self.log['DECHARGE'].append(outline)
        
        # Get the charged residues in the sphere and the total charge of these residues in the sphere
        for chain in self.PDB:
            for key in self.PDB[chain]:
                at = self.PDB[chain][key]
                if chain not in decharge:
                    continue
        
                if at[6] in decharge[chain]:
                    at[4] = charged_res[at[4]][0]
                    continue

                else: # at[4] in charged_res:
                    if at[4] in charged_res:
                        if at[2].strip() == charged_res[at[4]][1]:
                            self.log['CHARGE'].append('{:<10}{:<10}{:<10}{:<10}'.format(
                                                              self.log['QRESN'][chain][at[6]],
                                                                                        at[6],
                                                                                        chain,
                                                                                        at[4]))
                            self.log['TOTAL_CHARGE'] += charged_res[at[4]][2]                            

    def set_OXT(self):
        ## NOTE WARNING, ETC: Q AA codes are sometimes 3, sometimes 4, they MUST be updated to pdb
        ## standards ASAP!!!
        decharged = ['LYN','ARN','GLH','ASH']
        CTERM = []
        remove = []
        for chain in self.PDB:
            for key in self.PDB[chain]:
                at = self.PDB[chain][key]
                if at[2].strip() == 'OXT':
                    CTERM.append((at[6],at[5]))
                    self.log['CTERM'].append('{} {}'.format(at[6], at[4]))

                if self.origin == 'gromacs':
                    if at[2].strip() == 'O1':
                        CTERM.append((at[6],at[5]))
                        self.log['CTERM'].append('{} {}'.format(at[6], at[4]))
            
            for cterm in CTERM:
                for key in self.PDB[chain]:
                    at = self.PDB[chain][key]
                    
                    #HOTFIX, check if this is ok
                    if at[4] in decharged:
                        if at[2] == 'O1':
                            self.PDB[chain][key][2] = 'O'

                        if at[2] == 'O2':
                            remove.append([chain,key])
                                
                    else:
                        if at[6] == cterm[0] and at[5] == cterm[1]:
                            self.PDB[chain][key][4] = 'C' + at[4]

                        if self.origin == 'gromacs':
                            if at[2] == 'O1':
                                self.PDB[chain][key][2] = 'O'

                            if at[2] == 'O2':
                                self.PDB[chain][key][2] = 'OXT'
        
        # Remove atoms from c or n terminal decharged residues
        for at in remove:
            del self.PDB[at[0]][at[1]]                            

            
            
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
                            outline = '{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}'.format(self.log['QRESN'][chain][cys_mat[k][0]], 
                                                                        self.log['QRESN'][chain][cys_mat[j][0]],
                                                                        cys_mat[k][0],
                                                                                    chain,
                                                                        cys_mat[j][0],
                                                                                   chain)
                            self.log['CYX'].append(outline)

                    i += 1
                    
                
                for chain in self.PDB:
                    for key in self.PDB[chain]:
                        at = self.PDB[chain][key]
                        if at[6] in cyx and at[4] == 'CYS':
                            self.PDB[chain][key][4] = 'CYX'

            except:
                return None            
            

            
    def write_tmpPDB(self):
        with open('try67.pdb'[:-4] + '_tmp.pdb', 'w') as outfile:
            for chain in self.PDB:
                atom_numbers = sorted(list(self.PDB[chain].keys()))
                for atom in atom_numbers:
                    outline = IO.pdb_parse_out(self.PDB[chain][atom]) + '\n'
                    outfile.write(outline)
#                    print(outline)
#                if len(self.PDB) != 1:
#                    outfile.write('GAP\n')            
            
    def write_qprep(self):
        replacements = {'FF_LIB'    :   s.FF_DIR + '/OPLS2015.lib',
                        'FF_PRM'    :   s.FF_DIR + '/OPLS2015.prm',
                        'PROTPDB'   :   'try67.pdb'[:-4] + '_tmp.pdb',
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
#                        print(outline)
            
    def run_qprep(self):
        qprep = s.Q_DIR[self.preplocation] + 'Qprep6'
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
#                        print(outline)

                    if line[4] not in waters:
 #                       print(line)
                        outline = IO.pdb_parse_out(line) + '\n'
                        protout.write(outline)
                if line[0:3] == 'TER' or line == 'GAP':
                    protout.write(line)
#                    print(line)
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
            outfile.write('{:10}{:10}{:10}{:10}\n'.format('Q_RESN', 'PDB_IN', 'CHAIN','RESNAME'))
            for line in self.log['DECHARGE']:
                outfile.write(line + '\n')
                
            outfile.write('--------------------------------------------------------------------')
            outfile.write('\nThe following charged residues are in the sphere:\n')
            outfile.write('{:10}{:10}{:10}{:10}\n'.format('Q_RESN', 'PDB_IN', 'CHAIN','RESNAME'))
            for line in self.log['CHARGE']:
                outfile.write(line + '\n')
            
            outfile.write('--------------------------------------------------------------------')  
            outfile.write('\n')
            outfile.write('The following S-S bonds have been found:\n')
            outfile.write('{:10}{:10}{:10}{:10}{:10}{:10}\n'.format(
                'Q_CYS1', 'Q_CYS2', 'PDB_CYS1', 'CHAIN','PDB_CYS2','CHAIN'))
            for line in self.log['CYX']:
                outfile.write(line + '\n')
            
            outfile.write('--------------------------------------------------------------------')  
            outfile.write('\n')
            outfile.write('The following is a mapping of residue numbers in Q and the input \n')
            outfile.write('pdbfile "{}":\n'.format(self.prot))
            outfile.write('{:10}{:10}{:10}{:10}\n'.format('Q_RESN', 'PDB_IN', 'CHAIN','RESNAME'))
            for chain in self.PDB:
                for key in self.PDB[chain]:
                    if self.PDB[chain][key][2].strip() == 'CB':
                        PDB_resi = self.PDB[chain][key][6]
                        Q_resi = self.log['QRESN'][chain][PDB_resi]
                        resn = self.PDB[chain][key][4]
                        outfile.write('{:<10}{:<10}{:10}{:10}\n'.format(Q_resi,PDB_resi,chain,resn))
    

    def cleanup(self):
        if self.noclean == False:
            os.remove('test_pdb.pdb'[:-4] + '_tmp.pdb')
            os.remove('test_pdb.pdb'[:-4] + '_noH.pdb')
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
                        default = 'gromacs',
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


                
    
    
    

    run.get_center_coordinates()        # 00
    run.spheremaker()
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
#    run.cleanup()                       # 11
    run.spheremaker()                   #12
