#!/usr/bin/env python

import argparse
from pickle import TRUE
import re
import os
import stat
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
    def __init__(self, cofactor, include, forcefield,
                 cluster, temperature, replicates, production,
                 preplocation, timestep, *args, **kwargs):
        
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
            
        if temperature == None:
            self.temperature = s.TEMPERATURE
        else:
            self.temperature = temperature
            
        if replicates == None:
            self.replicates = s.REPLICATES
        else:
            self.replicates = replicates
        
        if production == None:
            self.production = 5
        else:
            self.production = production

        self.include = include
        self.forcefield = forcefield
        self.preplocation = preplocation
        self.CYX = []
        self.PDB = {}
        self.PDB2Q = {}
        self.systemsize = 0
        self.cluster = cluster
        self.timestep = timestep
        
        self.replacements = {'ATOM_START_LIG1':'1',
                             'WATER_RESTRAINT':'',
                             'TEMP_VAR':self.temperature,
                             'RUN_VAR':self.replicates,
                             'RUNFILE':'run' + self.cluster + '.sh',
                             'EQ_LAMBDA': '',
                             'FLOAT_LAMBDA1': '',
                             'FLOAT_LAMBDA2': '',
                             'DIST': ''
        }

        
    def create_environment(self):
        self.directory = 'sampling'
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
            os.makedirs(self.directory + '/inputfiles')
            
        files = {self.directory + '/inputfiles/' + self.forcefield + '.lib': \
                 s.FF_DIR + '/' + self.forcefield + '.lib'}
        
        if self.cofactor != None:
            for filename in self.cofactor:
                files[self.directory + '/inputfiles/' + filename + '.lib'] = filename + '.lib'
        
        for key in files:
            shutil.copy(files[key], key)
    
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
            
        prms = IO.read_prm(prmfiles)
        prm_merged = self.directory + '/inputfiles/' + self.forcefield + '_merged.prm'
        self.prm_merged = self.forcefield + '_merged.prm'

        with open (prm_merged, 'w') as outfile:
            for key in headers:
                outfile.write(key + '\n')
                for line in prms[key]:
                    outfile.write(line)
    
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

                        try:
                            self.PDB[line[6]].append(line)
                        except:
                            self.PDB[line[6]] = [line]                        
                        self.systemsize += 1
                
                resoffset = len(self.PDB)
                        
                    
    def write_pdb(self):
        PDBout = self.directory + '/inputfiles/complex.pdb'
        self.PDBout = 'complex.pdb'
        with open(PDBout, 'w') as outfile:
            for key in self.PDB:
                for line in self.PDB[key]:
                    outline = IO.pdb_parse_out(line) + '\n'
                    outfile.write(outline)
            
                    
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
        replacements['SOLVENT'] = '4 water.pdb'

        
        src = s.INPUT_DIR + '/qprep_resFEP.inp'
        self.qprep = self.directory + '/inputfiles/qprep.inp'
        libraries = [self.forcefield + '.lib']
        if self.cofactor != None:
            for filename in self.cofactor:
                libraries.append(filename + '.lib')
                
        with open(src) as infile:
            with open(self.qprep, 'w') as outfile:    
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
        IO.run_command(qprep, options, string = True)
        os.chdir('../../')

    def write_EQ(self):
        self.replacements['ATOM_END'] = '{}'.format(self.systemsize)
        for EQ_file in glob.glob(s.INPUT_DIR + '/eq*.inp'):
            src = EQ_file
            EQ_file = EQ_file.split('/')
            tgt = self.directory + '/inputfiles/' + EQ_file[-1]
            with open(src) as infile:
                with open(tgt, 'w') as outfile:
                    for line in infile:
                        if 'fep' in line or 'lambdas' in line:
                            continue
                        outline = IO.replace(line, self.replacements)
                        outfile.write(outline)
                    
    def write_MD(self):
        
        for i in range(1, int(self.production) + 1):
            src = s.INPUT_DIR + '/md_XXXX_XXXX.inp'
            tgt = self.directory + '/inputfiles/md_{}.inp'.format(i)
            self.replacements['FILE'] = 'md_{}'.format(i)
            
            with open(src) as infile:
                with open(tgt, 'w') as outfile:
                    for line in infile:
                        if 'fep' in line or 'lambdas' in line:
                            continue
                        outline = IO.replace(line, self.replacements)
                        outfile.write(outline)

            ## Store previous file
            self.replacements['FILE_N'] = 'md_{}'.format(i)


    def write_runfile(self):
        src = s.INPUT_DIR + '/run_sample.sh'
        tgt = self.directory + '/inputfiles/run' + self.cluster + '.sh'
        EQ_files = sorted(glob.glob(self.directory + '/inputfiles/eq*.inp'))
        MD_files = sorted(glob.glob(self.directory + '/inputfiles/md*.inp'))
                
        replacements = IO.merge_two_dicts(self.replacements, getattr(s, self.cluster))
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
                        for line in MD_files:
                            file_base = line.split('/')[-1][:-4]
                            outline = 'time srun $qdyn {}.inp'  \
                                        ' > {}.log\n'.format(file_base,
                                                            file_base)
                            outfile.write(outline)

    def settimestep(self):
        if self.timestep == '1fs':
            self.replacements['NSTEPS1'] = '1000000'
            self.replacements['NSTEPS2'] = '2000000'
            self.replacements['STEPSIZE'] = '1.0'
            self.replacements['STEPTOGGLE'] = 'off'
            
        if self.timestep == '2fs':
            self.replacements['NSTEPS1'] = '500000'
            self.replacements['NSTEPS2'] = '1000000'
            self.replacements['STEPSIZE'] = '2.0'
            self.replacements['STEPTOGGLE'] = 'on' 
                        
    def write_submitfile(self):
        submit_in = s.ROOT_DIR + '/INPUTS/FEP_submit_benchmark.sh'
        submit_out = self.directory + ('/FEP_submit.sh')
        with open(submit_in) as infile, open (submit_out, 'w') as outfile:
            for line in infile:
                if 'restart' in line:
                    continue
                line = IO.replace(line, self.replacements)
                outfile.write(line)

        try:
            st = os.stat(submit_out)
            os.chmod(submit_out, st.st_mode | stat.S_IEXEC)

        except:
            print("WARNING: Could not change permission for " + submit_out)
                            
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
    
    parser.add_argument('-f', '--forcefield',
                        dest = "forcefield",
                        required = True,
                        choices = ['OPLS2015', 'OPLS2005', 'SIDECHAIN', 'AMBER14sb','CHARMM36'],
                        help = "Forcefield to use.")
    
    parser.add_argument('-p', '--production',
                        dest = "production",
                        required = False,
                        help = "Sampling length of production phase x2 nanoseconds. Default is 5 (x2 ns = 10 ns).")
    
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
    
    parser.add_argument('-ts', '--timestep',
                        dest = "timestep",
                        choices = ['1fs','2fs'],
                        default = "2fs",
                        help = "Simulation timestep, default 2fs"
                       )

    args = parser.parse_args()
    run = Run(cofactor = args.cofactor,
              forcefield = args.forcefield,
              production = args.production,
              cluster = args.cluster,
              temperature = args.temperature,
              replicates = args.replicates,
              preplocation = args.preplocation,
              timestep = args.timestep,              
              include = ('ATOM','HETATM'))
    
    run.create_environment()            # 01
    run.read_input()                    # 02    
    run.readpdb()                       # 03
    run.merge_prm()                     # 06
    run.write_pdb()                     # 07
    run.select_waters()                 # 08
    run.write_qprep()                   # 09
    run.run_qprep()                     # 10
    run.settimestep()                   # 12
    run.write_EQ()                      # 13
    run.write_MD()                      # 14
    run.settimestep()                   # 15
    run.write_submitfile()              # 16
    run.write_runfile()                 # 17
