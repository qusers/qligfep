#!/usr/bin/env python

import argparse
from pickle import TRUE
import re
import os
import shutil
from ssl import HAS_TLSv1_3
import sys
import glob
import subprocess
from decimal import Decimal
import numpy as np

import functions as f
import settings as s
import IO

class Run(object):
    """
    Setup residue FEPs using either a single or dual topology approach
    """
    def __init__(self, cofactor, mutation, include, forcefield, windows,
                 sampling, system, cluster, temperature, replicates, dual, 
                 preplocation, start, timestep, mutchain, *args, **kwargs):
        
        self.cofactor = cofactor
        # Check whether all required files are there:
        required = ['protein.pdb', 'water.pdb', 'protPREP.log']
        if self.cofactor != None:
            extension = ['.pdb', '.lib', '.prm']
            for filename in self.cofactor:
                for line in extension:
                    required.append(filename + line)
        for filename in required:
            if os.path.exists(filename) == False:
                required = ' '.join(required)
                print("ERROR: {} are required files, exiting now.".format(required))
                sys.exit()
        
        if windows == None:
            self.windows = s.WINDOWS
            
        else:
            self.windows = int(windows)
            
        if sampling == None:
            self.sampling = s.SAMPLING
            
        else:
            self.sampling = sampling
            
        if temperature == None:
            self.temperature = s.TEMPERATURE
            
        else:
            self.temperature = temperature
            
        if replicates == None:
            self.replicates = s.REPLICATES
            
        else:
            self.replicates = replicates
        
        tmp = mutation.split(':')
        if len(tmp) == 2:
            self.mutation = re.split('(\d+)', tmp[1])
        else:
            self.mutation = re.split('(\d+)', tmp[0])

        self.include = include
        self.forcefield = forcefield
        self.preplocation = preplocation
        self.lambdas = []
        self.CYX = []
        self.PDB = {}
        self.PDB2Q = {}
        self.chain = mutchain
        self.systemsize = 0
        self.system = system
        self.cluster = cluster
        self.FEPlist = []
        self.dual = dual
        self.start = start
        self.dualMUT = {'0':[], '1':[]}
        self.nonAA = False
        self.timestep = timestep
        
        self.atoms = {'CA': 0,
                      'HA': 0, 'ha': 0, 'HA2': 0, 'ha2': 0, 'HA3': 0, 'ha3': 0,
                      'CB': 0, 'cb': 0,
                      'HB': 0, 'hb': 0, 'HB1': 0, 'hb1': 0, 'HB2': 0, 'hb2': 0, 'HB3': 0, 'hb3': 0,
                      'CG': 0, 'cg': 0, 'CG1': 0, 'cg1': 0, 'CG2': 0, 'cg2': 0,
                      'OG': 0, 'og': 0, 'OG1': 0, 'og1': 0,
                      'SG': 0, 'sg': 0,
                      'HG': 0, 'hg': 0, 'HG1': 0, 'hg1': 0, 'HG2': 0, 'hg2': 0, 'HG3': 0, 'hg3': 0,
                      'HG11': 0, 'hg11': 0, 'HG12': 0, 'hg12': 0, 'HG13': 0, 'hg13': 0,
                      'HG21': 0, 'hg21': 0, 'HG22': 0, 'hg22': 0, 'HG23': 0, 'hg23': 0,
                      'CD': 0, 'cd': 0, 'CD1': 0, 'cd1': 0, 'CD2': 0, 'cd2': 0,
                      'ND1': 0, 'nd1': 0, 'ND2': 0, 'nd2': 0,
                      'OD1': 0, 'od1': 0, 'OD2': 0, 'od2': 0,
                      'SD': 0, 'sd': 0,
                      'HD': 0, 'hd': 0, 'HD1': 0, 'hd1': 0, 'HD2': 0, 'hd2': 0, 'HD3': 0, 'hd3': 0,
                      'HD11': 0, 'hd11': 0, 'HD12': 0, 'hd12': 0, 'HD13': 0, 'hd13': 0,
                      'HD21': 0, 'hd21': 0, 'HD22': 0, 'hd22': 0, 'HD23': 0, 'hd23': 0,
                      'CE': 0, 'ce': 0, 'CE1': 0, 'ce1': 0, 'CE2': 0, 'ce2': 0, 'CE3': 0, 'ce3': 0,
                      'OE1': 0, 'oe1': 0, 'OE2': 0, 'oe2': 0,
                      'NE': 0, 'ne': 0, 'NE1': 0, 'ne1': 0, 'NE2': 0, 'ne2': 0,
                      'HE': 0, 'he': 0, 'HE1': 0, 'he1': 0, 'HE2': 0, 'he2': 0, 'HE3': 0, 'he3': 0,
                      'HE21': 0, 'he21': 0, 'HE22': 0, 'he22': 0,
                      'CZ': 0, 'cz': 0, 'CZ2': 0, 'cz2': 0, 'CZ3': 0, 'cz3': 0,
                      'NZ': 0, 'nz': 0,
                      'HZ': 0, 'hz': 0, 'HZ1': 0, 'hz1': 0, 'HZ2': 0, 'hz2': 0, 'HZ3': 0, 'hz3': 0,
                      'CH2': 0, 'ch2': 0,
                      'OH': 0, 'oh': 0,
                      'NH1': 0, 'nh1': 0, 'NH2': 0, 'nh2': 0,
                      'HH': 0, 'hh': 0, 'HH2': 0, 'hh2': 0, 
                      'HH11': 0, 'hh11': 0, 'HH12': 0, 'hh12': 0,
                      'HH21': 0, 'hh21': 0, 'HH22': 0, 'hh22': 0}
        
        self.replacements = {'ATOM_START_LIG1':'1',
                             'WATER_RESTRAINT':'',
                             'TEMP_VAR':self.temperature,
                             'RUN_VAR':self.replicates,
                             'RUNFILE':'run' + self.cluster + '.sh'  
                            }
        
        if self.start == '1' or self.start == '0':
            self.replacements['EQ_LAMBDA'] = '1.000 0.000'
            
        if self.start == '0.5':
            self.replacements['EQ_LAMBDA'] = '0.500 0.500'
    
    def checkFEP(self):
        if len(self.mutation[0]) == 1:
            AA_from = IO.AA(self.mutation[0])
        else:
            AA_from = self.mutation[0]
        
        try:
            if len(self.mutation[2]) == 1:
                AA_to = IO.AA(self.mutation[2])
            else:
                AA_to = self.mutation[2]
                
        except:
            print("residue not found in library, assuming non-natural AA")
            AA_to = 'X'
            self.nonAA = True
        mutation = '{}{}{}'.format(*self.mutation)
        FEPdir = s.FF_DIR + '/.FEP/' + self.forcefield +  '/' + AA_from + '_' + AA_to # + '_CI' # Include for charge correction
            
        self.FromGly, self.ToGly = False, False
        if self.mutation[0] == 'G' or self.mutation[0] == 'GLY':
            self.FromGly = True
            self.replacements['EQ_LAMBDA'] = '0.000 1.000'
        if self.mutation[2] == 'G' or self.mutation[2] == 'GLY':
            self.ToGly = True
        
        if self.dual == True:
            print('Generating dual topology inputfiles')
            self.FEPlist = ['FEP1.fep', 'FEP2.fep']
            self.FEPdir = None
        
        elif not os.path.exists(FEPdir):
            print('FATAL: no FEP files found for the {} mutation' \
                      'in {} exiting now.'.format(mutation,
                                                  FEPdir))
            sys.exit()
        
        else:
            for line in sorted(glob.glob(FEPdir + '/FEP*.fep')):
                line = line.split('/')
                self.FEPlist.append(line[-1])
            self.FEPdir = FEPdir
        
    def create_environment(self):
        self.directory = 'FEP_{}{}{}'.format(*self.mutation)
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
            os.makedirs(self.directory + '/inputfiles')
            
        files = {self.directory + '/inputfiles/' + self.forcefield + '.lib': \
                 s.FF_DIR + '/' + self.forcefield + '.lib',
                }
        
        if self.cofactor != None:
            for filename in self.cofactor:
                files[self.directory + '/inputfiles/' + filename + '.lib'] = filename + '.lib'
        
        for key in files:
            shutil.copy(files[key], key)
    
    def read_input(self):
        block = 0
        with open('protPREP.log') as infile:
            for line in infile:
                line = line.split()
                if len(line) > 1:
                    if line[1] == 'center:':
                        self.sphere = line[2:]
                        
                    if line[1] == 'radius:':
                        self.radius = line[2]
                        self.replacements['SPHERE'] = self.radius
                        
                    if line[0] == 'Q_CYS1':
                        block = 1
                        
                    if line[0] == 'pdbfile':
                        block = 2
                        
                    if line[0][0] == '-':
                        block = 0
                        
                    if block == 1:
                        if line[0].isdigit():
                            self.CYX.append([line[0], line[1]])

                    if block == 2 and len(line) ==  3:
                        chain = ' '
                        if chain not in self.PDB2Q:
                            self.PDB2Q[chain] = {}
                        self.PDB2Q[chain][line[1]] = line[0]    
                    if block == 2 and len(line) == 4:
                        chain = line[2]
                        if chain not in self.PDB2Q:
                            self.PDB2Q[chain] = {}
                        self.PDB2Q[chain][line[1]] = line[0]
    def readpdb(self):
        atnr = 0
        if len(self.mutation[2]) == 1:
            if self.nonAA != True:
                self.MUTresn = '{}2{}'.format(self.mutation[0],self.mutation[2])
            else:
                self.MUTresn = '{}2{}'.format(self.mutation[0],'X')
        
        elif len(self.mutation[2]) == 3:
            self.MUTresn = '{}2{}'.format(IO.AA(self.mutation[0]),IO.AA(self.mutation[2]))
        if self.ToGly == True or self.FromGly == True:
            self.backbone = ['C', 'O', 'CA', 'N', 'H']
        else: 
            self.backbone = ['C', 'O', 'CA', 'N', 'H', 'HA']
        # Read in the mutant residue
        if self.dual == True:
            MUT = []
            if len(self.mutation[2]) == 1:
                mut_res = '{}{}.pdb'.format(IO.AA(self.mutation[2]),self.mutation[1])
            else:
                mut_res = '{}{}.pdb'.format(self.mutation[2],self.mutation[1])
            with open(mut_res) as infile:   # Reads the mutant residue pdb as MUTX.pdb
                for line in infile:
                    if line.startswith(self.include):
                        line = IO.pdb_parse_in(line)
                        # The residue name should match the mutant residue number
                        # relative to Q numbering from protPREP.log
                        line[6] = int(self.PDB2Q[self.chain][self.mutation[1]])
                        if line[2] in self.backbone:
                            continue
                        else:
                            MUT.append(line)
        
        pdb_files = ['protein.pdb']
        if self.cofactor != None:
            for line in self.cofactor:
                pdb_files.append(line + '.pdb')
        for pdb_file in pdb_files:
            with open(pdb_file) as infile:
                for line in infile:
                    if line.startswith(self.include):
                        atnr += 1
                        line = IO.pdb_parse_in(line)
                        if pdb_file != 'protein.pdb':
                            line[6] = resoffset + 1
                            line[1] = self.systemsize
                            
                        line[1] = atnr
                        
                        # Change to hybrid residue name
                        if self.dual == True:
                            if str(line[6]) == self.PDB2Q[self.chain][self.mutation[1]]:
                                line[4] = self.MUTresn
                        
                        # Generate the dual topology residue in the PDB dictionary
                        # Construct the PDB dictionary            
                        try:
                            self.PDB[line[6]].append(line)
                        except:
                            self.PDB[line[6]] = [line]                        
                        
                        if self.dual == True:
                            if str(line[6]) == self.PDB2Q[self.chain][self.mutation[1]]:
                                if line[2] == 'O':
                                    for MUTat in MUT:
                                        atnr += 1
                                        if MUTat[2] not in self.backbone:
                                            # Atom name cannot be longer than 4 characters
                                            #if len(MUTat[2]) <= 3: 
                                            #    MUTat[2] = MUTat[2] + 'x'
                                                
                                            #elif len(MUTat[2]) == 4:
                                            #    MUTat[2] = MUTat[2][0:1] + 'x'
                                            MUTat[2] = MUTat[2].lower()
                                            
                                                
                                            MUTat[1] = atnr
                                            MUTat[4] = self.MUTresn

                                            try:
                                                self.PDB[MUTat[6]].append(MUTat)
                                            except:
                                                self.PDB[MUTat[6]] = [MUTat]

                                            self.dualMUT['0'].append(MUTat)
                                            self.systemsize += 1
                            
                                self.dualMUT['1'].append(line)
                        self.systemsize += 1
                
                resoffset = len(self.PDB)
    
    def create_dual_lib(self):
        FFlib = {}
        self.merged_lib = {}
        replacements = {}
        
        if self.dual == False:
            return None        
        
        headers =['[atoms]', 
                  '[bonds]',
                  '[connections]',
                  '[impropers]',
                  '[info]'
                 ]
        if self.forcefield[0:4] == 'OPLS':
            headers.append('[charge_groups]')   
        
        libfiles = [s.FF_DIR + '/' + self.forcefield + '.lib']
        if self.nonAA == True:
            libfiles.append(self.mutation[2] + '.lib')
        for libfile in libfiles:
            with open(libfile) as infile:
                for line in infile:
                    if line[0] == '*':
                        continue

                    if line[0] == '{':
                        RESname = line.split('}')[0][1:]
                        continue

                    try:
                        FFlib[RESname].append(line)

                    except:
                        FFlib[RESname] = [line]
        
        # Read the WT residue
        if len(self.mutation[0]) == 1:
            AA = IO.AA(self.mutation[0])
        else:
            AA = self.mutation[0]        
        for line in FFlib[AA]:
            if 'C223' in line:  #Change atom type of GLY 
                line = line.replace('C223','C224')
            if line.strip() in headers:
                header = line.strip()
                self.merged_lib[header] = []
                
            else:
                self.merged_lib[header].append(line)

        if self.FromGly == True:
            self.merged_lib['[bonds]'].append('       CA     ha\n')
        if self.ToGly == True:
            self.merged_lib['[bonds]'].append('       CA     ha2\n')
            self.merged_lib['[bonds]'].append('       CA     ha3\n')
        else:
            self.merged_lib['[bonds]'].append('       CA     cb\n')
        
        ## Read the mutant residue
        # Construct atom name replacements
        if len(self.mutation[2]) == 1:
            AA = IO.AA(self.mutation[2])
        else:
            AA = self.mutation[2]
        for line in FFlib[AA]:
            line = line.strip()
            line = line.split()
            if line != headers[0]:
                if len(line) > 2:
                    atom = line[1]
                    if atom not in self.backbone:
                        replacements[atom] = atom.lower()
                if line == header[1]:
                    break
                
        for line in FFlib[AA]:
            if line.strip() in headers:
                header = line.strip()
            if line[0] == '!':
                continue
            if len(line.split()) < 2:
                continue
 
            else:
                # connection flag only needs to be there once
                if header == headers[2]:
                    continue
                
                # Remove overlapping backbone definitions
                if line.split()[0] in self.backbone \
                or line.split()[1] in self.backbone:
                    continue
                
                # Merge the library on the reference
                line = IO.replace(line, replacements)
                self.merged_lib[header].append(line)
        
        # Write out the merged library file
        out = self.directory + '/inputfiles/' + self.MUTresn + '.lib'
        res_ID = int(self.PDB2Q[self.chain][self.mutation[1]])
        with open(out, 'w') as outfile:
            self.merged_lib['[charge_groups]'] = []
            cnt = 0
            outfile.write('{' + self.MUTresn + '}\n\n')
            for header in headers:
                outfile.write('\n    ' + header + '\n')
                if not header in self.merged_lib:
                    continue
                for line in self.merged_lib[header]:
                    if header == '[atoms]':
                        cnt += 1
                        line = line.split()
                        outfile.write('    {:2d} {:<10}{:<10}{:<10}\n'.format(cnt,
                                                                        line[1], 
                                                                        line[2], 
                                                                        line[3]))
                    else:
                        outfile.write(line)
                        
            
    def get_zeroforcebonded(self):
        if self.dual == False:
            return None    

        self.zerofk = {'[angles]':[],
                       '[bonds]':[],
                       '[impropers]':[],
                       '[torsions]':[]
                      }
        prmfile = [s.FF_DIR + '/' + self.forcefield + '.prm']
        libfile = self.directory + '/inputfiles/' + self.MUTresn + '.lib'
        headers =['[options]', 
                  '[atom_types]',
                  '[bonds]',
                  '[angles]',
                  '[torsions]',
                  '[impropers]'
                 ]
        
        # Write tmp qprep.inp and prm file
        if self.nonAA == True:
            prmfile.append(self.mutation[2] + '.prm')
        prms = IO.read_prm(prmfile)
        prm_tmp = self.directory + '/inputfiles/tmp.prm'
        qprep_tmp = self.directory + '/inputfiles/tmp.inp'
        pdb_tmp = self.directory + '/inputfiles/tmp.pdb'
        
        # Write out PDB file
        with open(pdb_tmp, 'w') as outfile:
            for item in self.PDB.items():
                for atom in item[1]:
                    outfile.write(IO.pdb_parse_out(atom) + '\n')
        
        with open(prm_tmp, 'w') as outfile:
            for key in headers:
                outfile.write(key + '\n')
                for line in prms[key]:
                    outfile.write(line)
                    
        with open(qprep_tmp, 'w') as outfile:
            outfile.write('rl {}\n'.format(libfile.split('/')[-1]))
            outfile.write('rl {}.lib\n'.format(self.forcefield))
            outfile.write('rprm {}\n'.format(prm_tmp.split('/')[-1]))
            outfile.write('rp {}\n'.format(pdb_tmp.split('/')[-1]))
            outfile.write('boundary 1 {} {} {} {}\n'.format(self.sphere[0],
                                                            self.sphere[1],
                                                            self.sphere[2],
                                                            self.radius))
            outfile.write('maketop tmp.top\n'.format(pdb_tmp))
            outfile.write('q')
            
        # run qprep to get bonded parameters to be set to zero
        os.chdir(self.directory + '/inputfiles/')
        cluster_options = getattr(s, self.preplocation)
        qprep = cluster_options['QPREP']
        options = ' < tmp.inp > tmp.out'
        # Somehow Q is very annoying with this < > input style so had to implement
        # another function that just calls os.system instead of using the preferred
        # subprocess module....
        IO.run_command(qprep, options, string = True)
        
        with open('tmp.out') as infile:
            for line in infile:
                line = line.split()
                if len(line) < 3:
                    continue
                    
                if line[1] != 'Missing':
                    continue
                
                # Create zero fk angles
                if line[2] == 'bond':
                    bond_prm = '{:11}{:11}{: 8.2f}{:>12.7}\n'.format(line[-2],
                                                                     line[-1],
                                                                     0.0,
                                                                     0.0)
                    if bond_prm not in self.zerofk['[bonds]']:
                        self.zerofk['[bonds]'].append(bond_prm)
                
                if line[2] == 'angle':
                    angle_prm = '{:11}{:11}{:11}{: 8.2f}{:>12.7}\n'.format(line[-3],
                                                                           line[-2],
                                                                           line[-1],
                                                                           0.0,
                                                                           110.0)
                    if angle_prm not in self.zerofk['[angles]']:
                        self.zerofk['[angles]'].append(angle_prm)
                    
                # Create zero fk torsions
                if line[2] == 'torsion':
                    torsion_prm='{:11}{:11}{:11}{:11}{:<10.3f}{: d}.000{:>10s}{:>10s}\n'.format(line[-4],
                                                                          line[-3],
                                                                          line[-2],
                                                                          line[-1],
                                                                          0.0,
                                                                          1,
                                                                          '0.000',
                                                                          '1')
                    if torsion_prm not in self.zerofk['[torsions]']:
                        self.zerofk['[torsions]'].append(torsion_prm)

        # ADD REMOVE FUNCTIONS FOR TMP FILES HERE
       # for tmp in glob.glob('tmp*'):
        #    os.remove(tmp)
        os.chdir('../../')
        
    def merge_prm(self):
        headers =['[options]', 
                  '[atom_types]',
                  '[bonds]',
                  '[angles]',
                  '[torsions]',
                  '[impropers]'
                 ]
        
        prmfiles = [s.FF_DIR + '/' + self.forcefield + '.prm']
        if self.cofactor != None:
            for filename in self.cofactor:
                prmfiles.append(filename + '.prm')
                
        if self.nonAA == True:
            prmfiles.append(self.mutation[2] + '.prm')
            
        prms = IO.read_prm(prmfiles)
        prm_merged = self.directory + '/inputfiles/' + self.forcefield + '_merged.prm'
        self.prm_merged = self.forcefield + '_merged.prm'

        with open (prm_merged, 'w') as outfile:
            for key in headers:
                outfile.write(key + '\n')
                for line in prms[key]:
                    outfile.write(line)
                if self.dual == True:
                    if key in self.zerofk:
                        outfile.write('! Zero order {} dual topology\n'.format(key.strip('[]')))
                        for line in self.zerofk[key]:
                            outfile.write(line)
                        outfile.write('\n')
                    
    def write_pdb(self):
        if self.system == 'protein':
            PDBout = self.directory + '/inputfiles/complex.pdb'
            self.PDBout = 'complex.pdb'
            with open(PDBout, 'w') as outfile:
                for key in self.PDB:
                    for line in self.PDB[key]:
                        outline = IO.pdb_parse_out(line) + '\n'
                        outfile.write(outline)
                        
        elif self.system == 'water':
            PDBout = self.directory + '/inputfiles/complex.pdb'
            self.PDBout = 'complex.pdb'
            with open(PDBout, 'w') as outfile:
                tpt = [int(self.PDB2Q[self.chain][self.mutation[1]]) + i for i in range(-1, 2)]
                resi = 1
                i, ter = 0, 0
                for res in tpt:
                    resi += 1
                    for line in self.PDB[res]:
                        i += 1
                        if line[2] == 'CA' and resi == 4:
                            ter = i
                        line[5] = ''
                        line[6] = resi
                        outline = IO.pdb_parse_out(line) + '\n'
                        outfile.write(outline)
            replacements = {'PDB': self.PDBout,
                            'COMPLEX': 'complex',
                            'TER': str(ter),
                            'TRIPEPTIDE': 'tripeptide.pdb'}
            with open(s.INPUT_DIR + '/tripeptide_2.pml') as pml_tmplt, \
            open(self.directory + '/inputfiles/tripeptide.pml', 'w') as pml_trgt:
                for line in pml_tmplt:
                    line = IO.replace(line, replacements)
                    pml_trgt.write(line)
                pml_tmplt.close()
                pml_trgt.close()
            
            cwd = self.directory + '/inputfiles'
            pymol = ['pymol', '-c', 'tripeptide.pml']
            
            result = subprocess.run(pymol, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if result.returncode == 0:
                print("pymol executed successfully")
            else:
                print("pymol error...")
                
            self.systemsize = 0
            with open(self.directory + '/inputfiles/tripeptide.pdb', 'r') as infile, \
                open(PDBout, 'w') as outfile:
                    for line in infile:
                        if line.startswith('ATOM'):
                            line = re.sub('1HB ', ' HB1', line)
                            line = re.sub('2HB ', ' HB2', line)
                            line = re.sub('3HB ', ' HB3', line)
                            line = re.sub('1HH3 ACE', 'HH31 ACE', line)
                            line = re.sub('2HH3 ACE', 'HH32 ACE', line)
                            line = re.sub('3HH3 ACE', 'HH33 ACE', line)
                            line = re.sub(' CH3 NME', ' CA  NMA', line)
                            line = re.sub('1HH3 NME', ' HA1 NMA', line)
                            line = re.sub('2HH3 NME', ' HA2 NMA', line)
                            line = re.sub('3HH3 NME', ' HA3 NMA', line)
                            line = re.sub('NME', 'NMA', line)
                            outfile.write(line)
                            self.systemsize += 1
                    infile.close()
                    outfile.close()
            
        elif self.system == 'vacuum':
            pass
            
                    
    def select_waters(self):
        src = 'water.pdb'
        tgt = self.directory + '/inputfiles/water.pdb'
        cofactor_coordinates = []
        waters_remove = []
        waters = {}
        
        if self.cofactor != None:
            for line in self.cofactor:
                with open(line + '.pdb') as infile:
                    for line in infile:
                        if line.startswith(self.include):
                            line = IO.pdb_parse_in(line)    
                            cofactor_coordinates.append([line[8], line[9], line[10]])
        
        with open (src) as infile:
            with open (tgt, 'w') as outfile:
                outfile.write('{} SPHERE\n'.format(self.radius))
                for line in infile:
                    line = IO.pdb_parse_in(line)    
                    coord_wat = [line[8], line[9], line[10]]
                    try:
                        waters[line[6]].append(line)
                    except:
                        waters[line[6]] = [line]

                    for coord_co in cofactor_coordinates:
                        if f.euclidian_overlap(coord_wat, coord_co, 1.6) == True:
                            waters_remove.append(line[6])
                
                for key in waters:
                    if key not in waters_remove:
                        for water in waters[key]:
                            outfile.write(IO.pdb_parse_out(water) + '\n')
                
    def write_qprep(self):
        replacements = {'PRM':self.prm_merged,
                        'PDB':self.PDBout,
                        'CENTER':'{} {} {}'.format(*self.sphere),
                        'SPHERE':self.radius,
                       }
        if self.system == 'protein':
            replacements['SOLVENT'] = '4 water.pdb'
            
        elif self.system == 'water':
            replacements['SOLVENT'] = '1 HOH'

        elif self.system == 'vacuum':
            replacements['solvate']='!solvate'
        
        src = s.INPUT_DIR + '/qprep_resFEP.inp'
        self.qprep = self.directory + '/inputfiles/qprep.inp'
        libraries = [self.forcefield + '.lib']
        if self.cofactor != None:
            for filename in self.cofactor:
                libraries.append(filename + '.lib')
                
        with open(src) as infile:
            with open(self.qprep, 'w') as outfile:
                if self.dual == True:
                    outfile.write('rl ' + self.MUTresn + '.lib\n')    
                for line in infile:
                    line = IO.replace(line, replacements)
                    if line.split()[0] == '!Added':
                        for libraryfile in libraries:
                            outfile.write('rl ' + libraryfile + '\n')
                        continue
                            
                    if line.split()[0] == '!addbond':
                        for CYX in self.CYX:
                            outline = 'addbond {}:SG {}:SG y\n'.format(*CYX)
                            outfile.write(outline)
                        continue
                    
                    outfile.write(line)
        
    def run_qprep(self):
        os.chdir(self.directory + '/inputfiles/')
        cluster_options = getattr(s, self.preplocation)
        qprep = cluster_options['QPREP']
        options = ' < qprep.inp > qprep.out'
        # Somehow Q is very annoying with this < > input style so had to implement
        # another function that just calls os.system instead of using the preferred
        # subprocess module....
        IO.run_command(qprep, options, string = True)
        os.chdir('../../')
                    
    def get_lambdas(self):
        self.lambdas = IO.get_lambdas(self.windows, self.sampling)

    def write_EQ(self):
        for line in self.PDB[int(self.PDB2Q[self.chain][self.mutation[1]])]:
            if line[2] == 'CA' and self.system == 'water'              \
            or line[2] == 'CA' and self.system == 'vacuum':             \
                self.replacements['WATER_RESTRAINT'] = '19 19 1.0 0 0'
                
        self.replacements['ATOM_END'] = '{}'.format(self.systemsize)
        
        res_ID = int(self.PDB2Q[self.chain][self.mutation[1]])

        top_p = self.directory + '/inputfiles/top_p.pdb'
        with open(top_p, 'r') as file:
                cnt = 0
                for line in file:
                    line = line.split()
                    residue_number = int(line[4])
                    if '2' in line[3]:
#                    if residue_number == res_ID:
                        self.PDB[res_ID][cnt][1] = int(line[1])
                        self.PDB[res_ID][cnt][2] = line[2]
                        self.PDB[res_ID][cnt][8] = line[5]
                        self.PDB[res_ID][cnt][9] = line[6]
                        self.PDB[res_ID][cnt][10] = line[7]
                        cnt += 1
                        if line[2] in self.atoms:
                            self.atoms[line[2]] = line[1] 
        
        wt, mut = IO.restraint_matrix(self.mutation)
        
        wt_ids = sorted([int(self.atoms[x]) for x in wt])
        mut_ids = sorted([int(self.atoms[x]) for x in mut])
        
        self.RES_WT = []
        self.RES_MUT = []
        
        hybrid_res = self.PDB[int(self.PDB2Q[self.chain][self.mutation[1]])]
#        for atom in hybrid_res:
#            if atom[1] in wt_ids:
#                if any(x in atom[2] for x in ['C', 'O', 'N', 'S']):
#                    self.RES_WT.append('{} {} {} {} 10.0 10.0 10.0 0'.format(atom[1], atom[8], atom[9], atom[10]))
#            elif atom[1] in mut_ids:
#                if any(x in atom[2].upper() for x in ['C', 'O', 'N', 'S']):
#                    self.RES_MUT.append('{} {} {} {} 10.0 10.0 10.0 0'.format(atom[1], atom[8], atom[9], atom[10]))                
#            else:
#                continue
        self.replacements['WT_RES'] = ''
        self.replacements['MUT_RES'] = ''
#        if self.FromGly == True:
#            self.replacements['MUT_RES'] = ''
#        else:
#            if mut:
#                self.replacements['MUT_RES'] = '\n'.join(self.RES_MUT)
#            else:
#                self.replacements['MUT_RES'] = ''

        
        self.dist_rest = []
        ha_pairs = IO.heavy_atom_match(wt, mut, hybrid_res)
        for ha_pair in ha_pairs:
            self.dist_rest.append('{} {} 0.0 0.5 10.0 0'.format(self.atoms[ha_pair[0]], self.atoms[ha_pair[1]]))
        self.replacements['DIST'] = '\n'.join(self.dist_rest)
        

        for EQ_file in glob.glob(s.INPUT_DIR + '/eq*.inp'):
            src = EQ_file
            EQ_file = EQ_file.split('/')
            tgt = self.directory + '/inputfiles/' + EQ_file[-1]
            with open(src) as infile:
                with open(tgt, 'w') as outfile:
                    for line in infile:
                        outline = IO.replace(line, self.replacements)
                        outfile.write(outline)
                    
    def write_MD(self):
        # Get restraint for water and vacuum systems
        for line in self.PDB[int(self.PDB2Q[self.chain][self.mutation[1]])]:
            if line[2] == 'CA' and self.system == 'water'              \
            or line[2] == 'CA' and self.system == 'vacuum':             \
                self.replacements['WATER_RESTRAINT'] = '19 19 1.0 0 0'
        
        if self.start == '1':
            for i in range(0, int(self.windows) + 1):
                j = int(self.windows) - i
                lambda1 = self.lambdas[i].replace(".", "")
                lambda2 = self.lambdas[j].replace(".", "")

                if self.lambdas[i] == '1.000':
                    src = s.INPUT_DIR + '/md_1000_0000.inp'
                    if self.FromGly == True:
                        self.replacements['FLOAT_LAMBDA1'] = self.lambdas[j]
                        self.replacements['FLOAT_LAMBDA2'] = self.lambdas[i]
                    else:
                        self.replacements['FLOAT_LAMBDA1'] = self.lambdas[i]
                        self.replacements['FLOAT_LAMBDA2'] = self.lambdas[j]
                    
                else:
                    src = s.INPUT_DIR + '/md_XXXX_XXXX.inp'
                    if self.FromGly == True:
                        self.replacements['FLOAT_LAMBDA1'] = self.lambdas[j]
                        self.replacements['FLOAT_LAMBDA2'] = self.lambdas[i]
                    else:
                        self.replacements['FLOAT_LAMBDA1'] = self.lambdas[i]
                        self.replacements['FLOAT_LAMBDA2'] = self.lambdas[j]
                
                tgt = self.directory + '/inputfiles/md_{}_{}.inp'.format(lambda1,
                                                                         lambda2)
                self.replacements['FILE'] = 'md_{}_{}'.format(lambda1,
                                                              lambda2)
#                self.replacements['MUT_RES'] = 'MUT_RES'
                              
                
                with open(src) as infile:
                    with open(tgt, 'w') as outfile:
                        for line in infile:
                            outline = IO.replace(line, self.replacements)
                            outfile.write(outline)

                ## Store previous file
                self.replacements['FILE_N'] = 'md_{}_{}'.format(lambda1,
                                                                lambda2)

        # This needs to be a function at one point      
        if self.start == '0.5':
            self.mdfiles = ['md_0500_0500']
            half = int((len(self.lambdas)-1)/2)
            lambda1 = self.lambdas[0:half]
            lambda2 = self.lambdas[half+1:]
            
            # Write out the middle file
            self.replacements['FILE_N'] = 'md_0500_0500'
            src = s.INPUT_DIR + '/md_0500_0500.inp'
            tgt = self.directory + '/inputfiles/md_0500_0500.inp'
            with open(src) as infile:
                with open(tgt, 'w') as outfile:
                    for line in infile:
                        outline = IO.replace(line, self.replacements)
                        outfile.write(outline)
            
            # Write out part 1 of the out files
            src = s.INPUT_DIR + '/md_XXXX_XXXX.inp'
            for i in range(0, len(lambda1)):
                l1 = lambda1[(half - i) - 1].replace(".", "")
                l2 = lambda2[i].replace(".", "")
                self.replacements['FLOAT_LAMBDA1'] = lambda1[(half - i) - 1]
                self.replacements['FLOAT_LAMBDA2'] = lambda2[i]
                tgt = self.directory + '/inputfiles/md_{}_{}.inp'.format(l1,
                                                                         l2)
                self.replacements['FILE'] = 'md_{}_{}'.format(l1,
                                                              l2)
                
                with open(src) as infile:
                    with open(tgt, 'w') as outfile:
                        for line in infile:
                            outline = IO.replace(line, self.replacements)
                            outfile.write(outline)
                        
                ## Store previous file
                self.replacements['FILE_N'] = 'md_{}_{}'.format(l1,
                                                                l2)
                
                self.mdfiles.append('md_{}_{}'.format(l1,l2))
            
            # And part 2
            self.replacements['FILE_N'] = 'md_0500_0500'
            for i in range(0, len(lambda1)):
                l1 = lambda2[i].replace(".", "")
                l2 = lambda1[(half - i) - 1].replace(".", "")
                self.replacements['FLOAT_LAMBDA1'] = lambda2[i]
                self.replacements['FLOAT_LAMBDA2'] = lambda1[(half - i) - 1]
                tgt = self.directory + '/inputfiles/md_{}_{}.inp'.format(l1,
                                                                         l2)
                self.replacements['FILE'] = 'md_{}_{}'.format(l1,
                                                              l2)
                
                with open(src) as infile:
                    with open(tgt, 'w') as outfile:
                        for line in infile:
                            outline = IO.replace(line, self.replacements)
                            outfile.write(outline)
                        
                ## Store previous file
                self.replacements['FILE_N'] = 'md_{}_{}'.format(l1,
                                                                l2)
                
                self.mdfiles.append('md_{}_{}'.format(l1,l2))


    def write_runfile(self):
        ntasks = getattr(s, self.cluster)['NTASKS']
        src = s.INPUT_DIR + '/run_benchmark.sh'
        tgt = self.directory + '/inputfiles/run' + self.cluster + '.sh'
        EQ_files = sorted(glob.glob(self.directory + '/inputfiles/eq*.inp'))
        
        if self.start == '1':
            MD_files = reversed(sorted(glob.glob(self.directory + '/inputfiles/md*.inp')))
            
        elif self.start == '0.5':
            MD_files = self.mdfiles[1:]
            half = int((len(MD_files)/2))
            md_1 = MD_files[0:half]
            md_2 = MD_files[half:]
        
        wt, mut = IO.restraint_matrix(self.mutation)
        
        wt_ids = sorted([int(self.atoms[x]) for x in wt])
        mut_ids = sorted([int(self.atoms[x]) for x in mut])
        
        
        if self.FromGly == True:
            self.replacements['RES_WT']  = repr('\n'.join(self.RES_MUT))
            self.replacements['RES_MUT'] = repr('')
        else:
            self.replacements['RES_WT']  = repr('\n'.join(self.RES_WT))
            if mut:
                self.replacements['RES_MUT'] = repr('\n'.join(self.RES_MUT))
            else:
                self.replacements['RES_MUT'] = repr('')
                
        
        replacements = IO.merge_two_dicts(self.replacements, getattr(s, self.cluster))
        replacements['FEPS'] = ' '.join(self.FEPlist)
        replacements['JOBNAME'] = 'QresFEP'
            
        with open(src) as infile:
            with open(tgt, 'w') as outfile:
                for line in infile:
                    if line.strip() == '#SBATCH -A ACCOUNT':
                        try:
                            replacements['ACCOUNT']
                        except:
                            line = ''
                    outline = IO.replace(line, replacements)
                    outfile.write(outline)
                    
                    if line.strip() == '#EQ_FILES':
                        for line in EQ_files:
                            file_base = line.split('/')[-1][:-4]
                            outline = 'time srun $qdyn {}.inp' \
                                       ' > {}.log\n'.format(file_base,
                                                           file_base)
                            outfile.write(outline)
                            
                    if line.strip() == '#RUN_FILES':
                        if self.start == '1':
                            for line in MD_files:
                                file_base = line.split('/')[-1][:-4]
                                outline = 'time srun $qdyn {}.inp'  \
                                          ' > {}.log\n'.format(file_base,
                                                               file_base)
                                outfile.write(outline)
                                
                        elif self.start == '0.5':
                            outline = 'time srun $qdyn {}.inp' \
                                       ' > {}.log\n\n'.format('md_0500_0500',
                                                              'md_0500_0500')
                            outfile.write(outline)
                            for i, md in enumerate(md_1):
                                outline1 = 'time srun $qdyn {}.inp'  \
                                          ' > {}.log &\n'.format(md_1[i],
                                                                 md_1[i])

                                outline2 = 'time srun $qdyn {}.inp'  \
                                          ' > {}.log\n'.format(md_2[i],
                                                               md_2[i])
                                outline3 = 'wait\n'

                                outfile.write(outline1)
                                outfile.write(outline2)
                                outfile.write(outline3)
                                outfile.write('\n')

    def settimestep(self):
        if self.timestep == '1fs':
            self.replacements['NSTEPS1'] = '100000'
            self.replacements['NSTEPS2'] = '10000'
            self.replacements['STEPSIZE'] = '1.0'
            self.replacements['STEPTOGGLE'] = 'off'
            
        if self.timestep == '2fs':
            self.replacements['NSTEPS1'] = '50000'
            self.replacements['NSTEPS2'] = '5000'
            self.replacements['STEPSIZE'] = '2.0'
            self.replacements['STEPTOGGLE'] = 'on' 
                        
    def write_submitfile(self):
        IO.write_submitfile_benchmark(self.directory, self.replacements)
    
    def write_FEPfile(self):
        if self.dual == True:
            run.write_dualFEPfile()
            
        if self.dual == False:
            run.write_singleFEPfile()
        
    def write_singleFEPfile(self):
        RES_list = self.PDB[int(self.PDB2Q[self.chain][self.mutation[1]])]
        RES_dic = {}
        block = 0
        empty = ['[change_charges]', 
                 '[FEP]',
                 '[atom_types]',
                 '[change_atoms]',
                 '[softcore]',
                 '[bond_types]'                
                ]
        atom_map = {}
        for line in RES_list:
            line[2:4] = [''.join(line[2:4])]
            line.insert(3, ' ')
            RES_dic[line[2].strip()] = line

        for src in sorted(glob.glob(self.FEPdir + '/FEP*.fep')):
            tgt = self.directory + '/inputfiles/' + src.split('/')[-1]
            with open (src) as infile:
                with open (tgt, 'w') as outfile:
                    for line in infile:
                        if line.strip() == '[atoms]':
                            outfile.write(line)
                            block = 1
                            
                        if line.strip() == '[change_bonds]':
                            outfile.write(line)
                            block = 2
                            
                        if line.strip() in empty:
                            block = 0
                        
                        if block == 0:
                            outfile.write(line)
                        
                        if block == 1:
                            if line.strip() != '[atoms]':
                                line = line.split()
                                if len(line) > 1 and line[0] != '!charge':
                                    atnr = RES_dic[line[3]][1]
                                    atom_map[line[0]] = atnr

                                    try:
                                        comment = line[4]

                                    except:
                                        comment = ' '
                                    outline = ('{:7s}{:<7d}{:5s}{:5s}{:5s}\n'.format(line[1],
                                                                                  atnr,
                                                                                  '!',
                                                                                  line[3],
                                                                                  comment))
                                
                                    outfile.write(outline)
                                else:
                                    outfile.write('\n')
                                
                        if block == 2:
                            if line.strip != '[change_bonds]':
                                line = line.split()
                                if len(line) > 1 and line[0] != '!charge':
                                    at_1 = atom_map[line[0]]
                                    at_2 = atom_map[line[1]]
                                    outline = '{:<7d}{:<7d}{:5s}{:5s}\n'.format(at_1,
                                                                             at_2,
                                                                             line[2],
                                                                             line[3]) 
                
                                    outfile.write(outline)
                                else:
                                    outfile.write('\n')                           
                            
    def write_dualFEPfile(self):
        # The vdW paramers need to be extracted from the reference .prm file
        prmfiles = [s.FF_DIR + '/' + self.forcefield + '.prm']
        if self.nonAA == True:
            prmfiles.append(self.mutation[2] + '.prm')
        vdw_prms = IO.read_prm(prmfiles)['[atom_types]']
        
        ## Next, we start to construct the FEP file
        fepfiles = self.FEPlist.copy()
        if self.FromGly:
            fepfiles.reverse()
        for i, fepfile in enumerate(self.FEPlist):
            with open(self.directory + '/inputfiles/' + fepfiles[i], 'w') as outfile:
                # First create the file header
                outfile.write('! Dual topology QresFEP: ' + self.MUTresn + '\n\n')
                outfile.write('[FEP]\n')
                outfile.write('states 2\n')
                outfile.write('softcore_use_max_potential on \n\n')

                # Write out [atoms] section
                outfile.write('[atoms]\n')
                cnt = 0
                res_ID = int(self.PDB2Q[self.chain][self.mutation[1]])
                for line in self.PDB[res_ID]:
                    cnt += 1
                    outfile.write('{:4d} {:4d} ! {:<4s}\n'.format(cnt, line[1], line[2]))

                outfile.write('\n')               

                # Write out [change_charges] section
                cnt = 0
                outfile.write('[change_charges]\n')

                for line in self.merged_lib['[atoms]']:
                    cnt += 1
                    line = line.split()                    
                    if line[1][0].islower()  == False:
                        if line[1] in self.backbone:
                            if fepfile == 'FEP1.fep' and line[1] == 'CA' and self.FromGly == True:
                                outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                            line[3], 
                                                                            '0.140'))
                            elif fepfile == 'FEP2.fep' and line[1] == 'CA' and self.FromGly == True:
                                outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                             '0.140', 
                                                                             '0.140'))                                                                
                            elif fepfile == 'FEP2.fep' and line[1] == 'CA' and self.ToGly == True:
                                outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                             '0.140', 
                                                                             '0.080'))                                
                            else:
                                outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                             line[3], 
                                                                             line[3]))                            
                        else:
                            if fepfile == 'FEP1.fep':
                                outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                             line[3], 
                                                                             '0.000'))
                            else:
                                outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                             '0.000', 
                                                                             '0.000'))

                    else:
                        if fepfile == 'FEP1.fep':
                            outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                         '0.000', 
                                                                         '0.000'))
                        else:
                            outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                         '0.000', 
                                                                         line[3]))

                outfile.write('\n')

                # Write out [atom_types] section
                outfile.write('[atom_types]\n')

                for line in vdw_prms:
                    line = line.split()
                    if len(line) >= 7:
                        outfile.write('{:<8s}{:>10s}{:>10s}{:>10s}'\
                                      '{:>10s}{:>10s}{:>10s}{:>10s}\n'.format(line[0],
                                                                            line[1],
                                                                            line[3], 
                                                                            '0.0000', 
                                                                            '0.0000',
                                                                            line[4],
                                                                            line[5],
                                                                            line[6],
                                                                            ))
                outfile.write('\n')

                # Write out [change_atoms] section
                outfile.write('[change_atoms]\n')
                cnt = 0
                for line in self.merged_lib['[atoms]']:
                    cnt += 1
                    line = line.split()
                    if line[1][0].islower()  == False:
                        if line[1] in self.backbone:
                            if fepfile == 'FEP1.fep' and line[1] == 'CA' and self.FromGly == True:
                                outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt,        
                                                                             'C224', 
                                                                             'C224'))
                            elif fepfile == 'FEP2.fep' and line[1] == 'CA' and self.FromGly == True:
                                outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                             'C224', 
                                                                             'C224'))              
                            elif fepfile == 'FEP2.fep' and line[1] == 'CA' and self.ToGly == True:
                                outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                             'C224', 
                                                                             'C224'))         
                            else:
                                outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                             line[2], 
                                                                             line[2]))
                        else:
                            if fepfile == 'FEP1.fep':
                                outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                             line[2], 
                                                                             line[2]))
                            else:
                                outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                             line[2], 
                                                                             'DUM'))

                    else:
                        if fepfile == 'FEP1.fep':    
                            outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                         'DUM', 
                                                                         line[2]))
                        else:
                            outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                         line[2], 
                                                                         line[2]))
                        
                outfile.write('\n')

                # Write out [softcore] section
                outfile.write('[softcore]\n')
                cnt = 0
                for line in self.merged_lib['[atoms]']:
                    cnt += 1
                    line = line.split()
                    if line[1][0].islower()  == False:
                        if line[1] in self.backbone:    
                             outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                          '0', 
                                                                          '0'))                       
                        else:
                            if fepfile == 'FEP1.fep':
                                outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                             '0', 
                                                                             '20'))
                            else:
                                outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                             '20', 
                                                                             '0'))

                    else:
                        if fepfile == 'FEP1.fep':
                            outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                         '0', 
                                                                         '20'))
                        else:
                            outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                         '20', 
                                                                         '0'))
                outfile.write('\n')
                
                # Write excluded pairs section
                outfile.write('[excluded_pairs]\n')
                res_A = []
                res_B = []
                for line in self.PDB[res_ID]:
                    if line[2] in self.backbone:
                        continue
                    if line[2].islower() == False:
                        res_A.append(line[1])
                    else:
                        res_B.append(line[1])
                for atom1 in res_A:
                    for atom2 in res_B:
                        outfile.write('{:6d}{:6d}{:>5s}{:>5s}\n'.format(atom1, atom2, '1', '1'))
                outfile.write('\n')
                
                # Add bond section - if we were to touch bond parameters (in case of changing Gly CA) 
                outfile.write('[bond_types]\n\n')
                outfile.write('[change_bonds]\n\n')
                
                
                # Add angles section
                outfile.write('[angle_types]\n\n')
                outfile.write('[change_angles]\n')
                
                # Set zero angles for wild-type side chain atoms connecting through common backbone CA to mutant side chain atoms for Gly mutations
                if self.FromGly == True or self.ToGly == True:
                    if self.FromGly == True:
                        outfile.write('{:6s}{:6s}{:6s}{:>5d}{:>5d}\n'.format(self.atoms['HA3'], self.atoms['CA'], self.atoms['cb'], 0, 0))
                        outfile.write('{:6s}{:6s}{:6s}{:>5d}{:>5d}\n'.format(self.atoms['HA2'], self.atoms['CA'], self.atoms['cb'], 0, 0))
                        outfile.write('{:6s}{:6s}{:6s}{:>5d}{:>5d}\n'.format(self.atoms['HA3'], self.atoms['CA'], self.atoms['ha'], 0, 0))
                        outfile.write('{:6s}{:6s}{:6s}{:>5d}{:>5d}\n'.format(self.atoms['HA2'], self.atoms['CA'], self.atoms['ha'], 0, 0))

                    elif self.ToGly == True:
                        outfile.write('{:6s}{:6s}{:6s}{:>5d}{:>5d}\n'.format(self.atoms['ha3'], self.atoms['CA'], self.atoms['CB'], 0, 0))
                        outfile.write('{:6s}{:6s}{:6s}{:>5d}{:>5d}\n'.format(self.atoms['ha2'], self.atoms['CA'], self.atoms['CB'], 0, 0))
                        outfile.write('{:6s}{:6s}{:6s}{:>5d}{:>5d}\n'.format(self.atoms['ha3'], self.atoms['CA'], self.atoms['HA'], 0, 0))
                        outfile.write('{:6s}{:6s}{:6s}{:>5d}{:>5d}\n'.format(self.atoms['ha2'], self.atoms['CA'], self.atoms['HA'], 0, 0))
                    outfile.write('\n')

                # Set zero angles for wild-type side chain atoms connecting through common backbone CA to mutant side chain atoms for non-Gly mutations
                if self.FromGly == False and self.ToGly == False:
                    outfile.write('{:6s}{:6s}{:6s}{:>5d}{:>5d}\n\n'.format(self.atoms['CB'], self.atoms['CA'], self.atoms['cb'], 0, 0))
                
                # Add torsions section
                outfile.write('[torsion_types]\n\n')
                outfile.write('[change_torsions]\n')
                
                # Set zero torsions for wild-type side chain atoms connecting through common backbone CA to mutant side chain atoms for Gly mutations
                if self.FromGly == True or self.ToGly == True:
                    if self.FromGly == True:
                        if self.mutation[2] in ['ASP', 'GLU', 'HID', 'HIE', 'ARG', 'LYS', 'PHE', 'LEU', 'MET', 'TRP', 'TYR', 'ASN', 'GLN']:
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['hb2'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA2'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['hb2'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA3'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['hb3'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA2'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['hb3'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA3'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['cg'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA2'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['cg'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA3'], 0, 0))

                        elif self.mutation[2] in ['ILE', 'VAL', 'THR']:
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['hb'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA2'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['hb'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA3'], 0, 0))
                            if self.mutation[2] == 'THR':
                                outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['og1'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA2'], 0, 0))
                                outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['og1'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA3'], 0, 0))
                            else:
                                outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['cg1'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA2'], 0, 0))
                                outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['cg1'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA3'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['cg2'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA2'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['cg2'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA3'], 0, 0))

                        elif self.mutation[2] in ['CYS', 'SER', 'ALA']:
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['hb2'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA2'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['hb2'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA3'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['hb3'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA2'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['hb3'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA3'], 0, 0))
                            if self.mutation[2] == 'CYS':
                                outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['sg'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA2'], 0, 0))
                                outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['sg'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA3'], 0, 0))
                            elif self.mutation[2] == 'SER':
                                outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['og'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA2'], 0, 0))
                                outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['og'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA3'], 0, 0))
                            else:
                                outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['hb1'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA2'], 0, 0))
                                outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['hb1'], self.atoms['cb'], self.atoms['CA'], self.atoms['HA3'], 0, 0))

                    elif self.ToGly == True:
                        if self.mutation[0] in ['ASP', 'GLU', 'HID', 'HIE', 'ARG', 'LYS', 'PHE', 'LEU', 'MET', 'TRP', 'TYR', 'ASN', 'GLN']:
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['HB2'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha2'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['HB2'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha3'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['HB3'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha2'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['HB3'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha3'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['CG'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha2'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['CG'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha3'], 0, 0))

                        elif self.mutation[0] in ['ILE', 'VAL', 'THR']:
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['HB'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha2'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['HB'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha3'], 0, 0))
                            if self.mutation[0] == 'THR':
                                outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['OG1'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha2'], 0, 0))
                                outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['OG1'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha3'], 0, 0))
                            else:
                                outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['CG1'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha2'], 0, 0))
                                outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['CG1'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha3'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['CG2'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha2'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['CG2'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha3'], 0, 0))

                        elif self.mutation[0] in ['CYS', 'SER', 'ALA']:
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['HB2'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha2'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['HB2'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha3'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['HB3'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha2'], 0, 0))
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['HB3'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha3'], 0, 0))
                            if self.mutation[0] == 'CYS':
                                outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['SG'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha2'], 0, 0))
                                outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['SG'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha3'], 0, 0))
                            elif self.mutation[0] == 'SER':
                                outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['OG'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha2'], 0, 0))
                                outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['OG'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha3'], 0, 0))
                            else:
                                outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['HB1'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha2'], 0, 0))
                                outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['HB1'], self.atoms['CB'], self.atoms['CA'], self.atoms['ha3'], 0, 0))
                    outfile.write('\n')
                    
                # Set zero torsions for wild-type side chain atoms connecting through common backbone CA to mutant side chain atoms for non-Gly mutations
                if self.FromGly == False and self.ToGly == False:
                    
                    if self.mutation[0] in ['ASP', 'GLU', 'HID', 'HIE', 'ARG', 'LYS', 'PHE', 'LEU', 'MET', 'TRP', 'TYR', 'ASN', 'GLN']:
                        outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['HB2'], self.atoms['CB'], self.atoms['CA'], self.atoms['cb'], 0, 0))
                        outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['HB3'], self.atoms['CB'], self.atoms['CA'], self.atoms['cb'], 0, 0))
                        outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['CG'], self.atoms['CB'], self.atoms['CA'], self.atoms['cb'], 0, 0))

                    elif self.mutation[0] in ['ILE', 'VAL', 'THR']:
                        outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['HB'], self.atoms['CB'], self.atoms['CA'], self.atoms['cb'], 0, 0))
                        if self.mutation[0] == 'THR':
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['OG1'], self.atoms['CB'], self.atoms['CA'], self.atoms['cb'], 0, 0))
                        else:
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['CG1'], self.atoms['CB'], self.atoms['CA'], self.atoms['cb'], 0, 0))
                        outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['CG2'], self.atoms['CB'], self.atoms['CA'], self.atoms['cb'], 0, 0))

                    elif self.mutation[0] in ['CYS', 'SER', 'ALA']:
                        outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['HB2'], self.atoms['CB'], self.atoms['CA'], self.atoms['cb'], 0, 0))
                        outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['HB3'], self.atoms['CB'], self.atoms['CA'], self.atoms['cb'], 0, 0))
                        if self.mutation[0] == 'CYS':
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['SG'], self.atoms['CB'], self.atoms['CA'], self.atoms['cb'], 0, 0))
                        elif self.mutation[0] == 'SER':
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['OG'], self.atoms['CB'], self.atoms['CA'], self.atoms['cb'], 0, 0))
                        else:
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['HB1'], self.atoms['CB'], self.atoms['CA'], self.atoms['cb'], 0, 0))

                    # Set zero torsions for mutant side chain atoms connecting through common backbone CA to wild-type side chain atoms
                    if self.mutation[2] in ['ASP', 'GLU', 'HID', 'HIE', 'ARG', 'LYS', 'PHE', 'LEU', 'MET', 'TRP', 'TYR', 'ASN', 'GLN']:
                        outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['hb2'], self.atoms['cb'], self.atoms['CA'], self.atoms['CB'], 0, 0))
                        outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['hb3'], self.atoms['cb'], self.atoms['CA'], self.atoms['CB'], 0, 0))
                        outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['cg'], self.atoms['cb'], self.atoms['CA'], self.atoms['CB'], 0, 0))

                    elif self.mutation[2] in ['ILE', 'VAL', 'THR']:
                        outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['hb'], self.atoms['cb'], self.atoms['CA'], self.atoms['CB'], 0, 0))
                        if self.mutation[2] == 'THR':
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['og1'], self.atoms['cb'], self.atoms['CA'], self.atoms['CB'], 0, 0))
                        else:
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['cg1'], self.atoms['cb'], self.atoms['CA'], self.atoms['CB'], 0, 0))
                        outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['cg2'], self.atoms['cb'], self.atoms['CA'], self.atoms['CB'], 0, 0))

                    elif self.mutation[2] in ['CYS', 'SER', 'ALA']:
                        outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['hb2'], self.atoms['cb'], self.atoms['CA'], self.atoms['CB'], 0, 0))
                        outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['hb3'], self.atoms['cb'], self.atoms['CA'], self.atoms['CB'], 0, 0))
                        if self.mutation[2] == 'CYS':
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['sg'], self.atoms['cb'], self.atoms['CA'], self.atoms['CB'], 0, 0))
                        elif self.mutation[2] == 'SER':
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['og'], self.atoms['cb'], self.atoms['CA'], self.atoms['CB'], 0, 0))
                        else:
                            outfile.write('{:6s}{:6s}{:6s}{:6s}{:4d}{:4d}\n'.format(self.atoms['hb1'], self.atoms['cb'], self.atoms['CA'], self.atoms['CB'], 0, 0))
                    outfile.write('\n')
            
                    
                    
    def write_qfep(self):
        qfep_in = s.ROOT_DIR + '/INPUTS/qfep.inp'
        qfep_out = self.directory + '/inputfiles/qfep.inp'
        i = 0
        total_l = len(self.lambdas)

        # TO DO: multiple files will need to be written out for temperature range
        kT = f.kT(float(self.temperature))
        replacements = {}
        replacements['kT']=kT
        replacements['WINDOWS']=str(self.windows)
        replacements['TOTAL_L']=str(total_l)
        
        with open(qfep_in) as infile:
            with open(qfep_out, 'w') as outfile:
                for line in infile:
                    line = IO.replace(line, replacements)
                    outfile.write(line)

                if line == '!ENERGY_FILES\n':
                    for i in range(0, total_l):
                        j = -(i + 1)
                        lambda1 = self.lambdas[i]
                        lambda2 = self.lambdas[j]
                        filename = 'md_' +                          \
                                    lambda1.replace('.', '') +      \
                                    '_' +                           \
                                    lambda2.replace('.', '') +      \
                                    '.en\n'

                        outfile.write(filename)
                            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='QresPREP',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Generate FEP files for resFEP, takes output from protPREP.py == ')

    
    parser.add_argument('-c', '--cofactor',
                        nargs='*',
                        dest = "cofactor",
                        required = False,
                        help = "PDB files of cofactors (e.g. ligand) with *.prm and *.lib files")
    
    parser.add_argument('-m', '--mutation',
                        dest = "mutation",
                        required = True,
                        help = "The desired mutation given as $WT$RESN$MUT,e.g. Y197A, resn is taken from original .pdb")
    
    parser.add_argument('-f', '--forcefield',
                        dest = "forcefield",
                        required = True,
                        choices = ['OPLS2015', 'OPLS2005', 'SIDECHAIN', 'AMBER14sb','CHARMM36'],
                        help = "Forcefield to use.")
    
    parser.add_argument('-s', '--sampling',
                        dest = "sampling",
                        required = False,
                        choices = ['linear',
                                   'sigmoidal', 
                                   'exponential', 
                                   'reverse_exponential'],
                        help = "Sampling type to be used, default taken from settings.py.")
    
    parser.add_argument('-w', '--windows',
                        dest = "windows",
                        required = False,
                        help = "Amount of windows to be used per FEP window, default taken from settings.py")    
    
    parser.add_argument('-S', '--system',
                        dest = "system",
                        required = True,
                        choices = ['protein', 'water', 'vacuum'],
                        help = "System type, can be protein, water or vacuum")
    
    parser.add_argument('-T', '--temperature',
                        dest = "temperature",
                        required = False,
                        help = "List of temperatures, default taken from settings.py")
    
    parser.add_argument('-r', '--replicates',
                        dest = "replicates",
                        required = False,
                        help = "Amount of replicates, default taken from settings.py")   
    
    parser.add_argument('-C', '--cluster',
                        dest = "cluster",
                        required = True,
                        help = "Cluster to use.")
    
    parser.add_argument('-P', '--preplocation',
                    dest = "preplocation",
                    default = 'CSB',
                    help = "define this variable if you are setting up your system elsewhere")
    
    parser.add_argument('-d', '--dual',
                        dest = "dual",
                        default = False,
                        action = 'store_true',
                        help = "Turn on if you want QresFEP to generate a dual topology hybrid")
    
    parser.add_argument('-l', '--start',
                        dest = "start",
                        default = '1',
                        choices = ['1', '0.5', '0'],
                        help = "Starting FEP in the middle or endpoint, 0.5 recommended for dual"
                       )
    
    parser.add_argument('-ts', '--timestep',
                        dest = "timestep",
                        choices = ['1fs','2fs'],
                        default = "2fs",
                        help = "Simulation timestep, default 2fs"
                       )

    parser.add_argument('-mc', '--mutchain',
                        dest = "mutchain",
                        default = ' ',
                        help = "Chain of the mutation"
                       )

    args = parser.parse_args()
    run = Run(mutation = args.mutation,
              cofactor = args.cofactor,
              forcefield = args.forcefield,
              sampling = args.sampling,
              windows = args.windows,
              system = args.system,
              cluster = args.cluster,
              temperature = args.temperature,
              replicates = args.replicates,
              preplocation = args.preplocation,
              dual = args.dual,
              start = args.start,
              timestep = args.timestep,              
              mutchain = args.mutchain,
              include = ('ATOM','HETATM'))
    
    run.checkFEP()                      # 00
    run.create_environment()            # 01
    run.read_input()                    # 02    
    run.readpdb()                       # 03
    run.create_dual_lib()               # 04
    run.get_zeroforcebonded()           # 05
    run.merge_prm()                     # 06
    run.write_pdb()                     # 07
    run.select_waters()                 # 08
    run.write_qprep()                   # 09
    run.run_qprep()                     # 10
    run.get_lambdas()                   # 11
    run.settimestep()                   # 12
    run.write_EQ()                      # 13
    run.write_MD()                      # 14
    run.settimestep()                   # 15
    run.write_submitfile()              # 16
    run.write_runfile()                 # 17
    run.write_FEPfile()                 # 18
    run.write_qfep()                    # 19

