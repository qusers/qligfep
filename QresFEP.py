import argparse
import re
import os
import shutil
import sys
import glob

import functions as f
import settings as s
import IO

class Run(object):
    """
    Setup residue FEPs using either a single or dual topology approach
    """
    def __init__(self, cofactor, mutation, include, forcefield, windows,
                 sampling, system, cluster, temperature, replicates, dual, 
                 preplocation, start, *args, **kwargs):
        
        self.cofactor = [cofactor]
        # Check whether all required files are there:
        required = ['protein.pdb', 'water.pdb', 'protPREP.log']
        if self.cofactor[0] != None:
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
        
        self.mutation = re.split('(\d+)', mutation)
        self.include = include
        self.forcefield = forcefield
        self.preplocation = preplocation
        self.lambdas = []
        self.CYX = []
        self.PDB = {}
        self.PDB2Q = {}
        self.systemsize = 0
        self.system = system
        self.cluster = cluster
        self.FEPlist = []
        self.dual = dual
        self.start = start
        self.dualMUT = {'0':[], '1':[]}
        
        self.replacements = {'ATOM_START_LIG1':'1',
                             'WATER_RESTRAINT':'',
                             'TEMP_VAR':self.temperature,
                             'RUN_VAR':self.replicates,
                             'RUNFILE':'run' + self.cluster + '.sh'  
                            }
        
        if self.start == '1':
            self.replacements['EQ_LAMBDA'] = '1.000 0.000'
            
        if self.start == '0.5':
            self.replacements['EQ_LAMBDA'] = '0.500 0.500'
    
    def checkFEP(self):
        if len(self.mutation[0]) == 1:
            AA_from = IO.AA(self.mutation[0])
            
        else:
            AA_from = self.mutation[0]
            
        AA_to = IO.AA(self.mutation[2])
        mutation = '{}{}{}'.format(*self.mutation)
        FEPdir = s.FF_DIR + '/.FEP/' + self.forcefield +  '/' + AA_from + '_' + AA_to
            
        if not os.path.exists(FEPdir):
            if self.dual == True:
                'Generating dual topology inputfiles'
                self.FEPlist = ['FEP1.fep']
                self.FEPdir = None
                
            else:    
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
        
        if self.cofactor[0] != None:
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
                        
                    if block == 2:
                        self.PDB2Q[line[1]] = line[0]

    def readpdb(self):
        atnr = 0
        self.MUTresn = '{}2{}'.format(self.mutation[0],self.mutation[2])
        self.backbone = ['C', 'O', 'CA', 'N', 'H', 'HA']
        # Read in the mutant residue
        if self.dual == True:
            MUT = []
            with open(self.mutation[2] + '.pdb') as infile:
                for line in infile:
                    line = IO.pdb_parse_in(line)
                    if line[2] in self.backbone:
                        continue
                        
                    else:
                        MUT.append(line)
        
        pdb_files = ['protein.pdb']
        if self.cofactor[0] != None:
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
                            if str(line[6]) == self.PDB2Q[self.mutation[1]]:
                                line[4] = self.MUTresn                        
                        
                        # Generate the dual topology residue in the PDB dictionary
                        # Construct the PDB dictionary            
                        try:
                            self.PDB[line[6]].append(line)

                        except:
                            self.PDB[line[6]] = [line]                        
                        
                        if self.dual == True:
                            if str(line[6]) == self.PDB2Q[self.mutation[1]]:
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
                            
                                self.dualMUT['1'].append(line)
                        
                        self.systemsize += 1
                
                resoffset = len(self.PDB)
    
    def create_dual_lib(self):
        # Perhaps we should create json .lib/.prm files at one point
        FFlib = {}
        self.merged_lib = {}
        replacements = {}
        
        if self.dual == False:
            return None        
        
        headers =['[atoms]', 
                  '[bonds]',
                  '[connections]',
                  '[impropers]',
                  '[charge_groups]',
                 ]   
        
        with open(s.FF_DIR + '/' + self.forcefield + '.lib') as infile:
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
        for line in FFlib[IO.AA(self.mutation[0])]:
            if line.strip() in headers:
                header = line.strip()
                self.merged_lib[header] = []
                
            else:
                self.merged_lib[header].append(line)
        
        # Add CA CBx bond
        self.merged_lib['[bonds]'].append('       CA     cb\n')
        
        ## Read the mutant residue
        # Construct atom name replacements
        for line in FFlib[IO.AA(self.mutation[2])]:
            line = line.strip()
            line = line.split()
            if line != headers[0]:
                if len(line) > 2:
                    atom = line[1]
                    replacements[atom] = atom.lower()
                
                if line == header[1]:
                    break
                
        for line in FFlib[IO.AA(self.mutation[2])]:
            if line.strip() in headers:
                header = line.strip()
                
            else:
                # connection flag only needs to be there once
                if header == headers[2]:
                    continue
                
                # Remove overlapping backbone definitions
                if 'O ' in line     \
                or 'C ' in line     \
                or 'CA ' in line     \
                or 'N ' in line     \
                or 'H ' in line     \
                or 'HA ' in line :
                    continue
                
                # Merge the library on the reference
                line = IO.replace(line, replacements)
                self.merged_lib[header].append(line)
        
        # Write out the merged library file
        out = self.directory + '/inputfiles/' + self.MUTresn + '.lib'
        with open(out, 'w') as outfile:
            cnt = 0
            outfile.write('{' + self.MUTresn + '}\n\n')
            for header in headers:
                outfile.write('\n    ' + header + '\n')
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

            
    def merge_prm(self):
        headers =['[options]', 
                  '[atom_types]',
                  '[bonds]',
                  '[angles]',
                  '[torsions]',
                  '[impropers]'
                 ]
        
        prmfiles = [s.FF_DIR + '/' + self.forcefield + '.prm']
        if self.cofactor[0] != None:
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
                    
    def write_pdb(self):
        PDBout = self.directory + '/inputfiles/complex.pdb'
        self.PDBout = 'complex.pdb'
        with open(PDBout, 'w') as outfile:
            for key in self.PDB:
                #print key
                for line in self.PDB[key]:
                    #print line
                    outline = IO.pdb_parse_out(line) + '\n'
                    outfile.write(outline)
                    
    def select_waters(self):
        src = 'water.pdb'
        tgt = self.directory + '/inputfiles/water.pdb'
        cofactor_coordinates = []
        waters_remove = []
        waters = {}
        
        if self.cofactor[0] != None:
            for line in self.cofactor:
                with open(line + '.pdb') as infile:
                    for line in infile:
                        if line.startswith(self.include):
                            line = IO.pdb_parse_in(line)    
                            cofactor_coordinates.append([line[8], line[9], line[10]])
        
        with open (src) as infile, open (tgt, 'w') as outfile:
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
        if self.cofactor[0] != None:
            for filename in self.cofactor:
                libraries.append(filename + '.lib')
                
        with open(src) as infile, open(self.qprep, 'w') as outfile:
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
        qprep = s.Q_DIR[self.preplocation] + 'qprep'
        options = ' < qprep.inp > qprep.out'
        # Somehow Q is very annoying with this < > input style so had to implement
        # another function that just calls os.system instead of using the preferred
        # subprocess module....
        IO.run_command(qprep, options, string = True)
        os.chdir('../../')
                    
    def get_lambdas(self):
        self.lambdas = IO.get_lambdas(self.windows, self.sampling)

    def write_EQ(self):
        for line in self.PDB[int(self.mutation[1])]:
            if line[2] == 'CA' and              \
            self.system == 'water'              \
            or self.system == 'vacuum':
                self.replacements['WATER_RESTRAINT'] = '{} {} 1.0 0 0'.format(line[1], 
                                                                              line[1])
        self.replacements['ATOM_END'] = '{}'.format(self.systemsize)
        
        for EQ_file in glob.glob(s.INPUT_DIR + '/eq*.inp'):
            src = EQ_file
            EQ_file = EQ_file.split('/')
            tgt = self.directory + '/inputfiles/' + EQ_file[-1]
            with open(src) as infile, open(tgt, 'w') as outfile:
                for line in infile:
                    outline = IO.replace(line, self.replacements)
                    outfile.write(outline)
                    
    def write_MD(self):
        if self.start == '1.0':
            for line in self.PDB[int(self.mutation[1])]:
                if line[2] == 'CA' and self.system == 'water' or self.system == 'vacuum':
                    self.replacements['WATER_RESTRAINT'] = '{} {} 1.0 0 0'.format(line[1], 
                                                                                  line[1])

            for i in range(0, int(self.windows) + 1):
                j = int(self.windows) - i
                lambda1 = self.lambdas[i].replace(".", "")
                lambda2 = self.lambdas[j].replace(".", "")

                if self.lambdas[i] == '1.000':
                    src = s.INPUT_DIR + '/md_1000_0000.inp'
                    
                else:
                    src = s.INPUT_DIR + '/md_XXXX_XXXX.inp'
                    self.replacements['FLOAT_LAMBDA1'] = self.lambdas[i]
                    self.replacements['FLOAT_LAMBDA2'] = self.lambdas[j]

                tgt = self.directory + '/inputfiles/md_{}_{}.inp'.format(lambda1,
                                                                         lambda2)
                self.replacements['FILE'] = 'md_{}_{}'.format(lambda1,
                                                              lambda2)

                with open(src) as infile, open(tgt, 'w') as outfile:
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
            with open(src) as infile, open(tgt, 'w') as outfile:
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
                
                with open(src) as infile, open(tgt, 'w') as outfile:
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
                
                with open(src) as infile, open(tgt, 'w') as outfile:
                    for line in infile:
                        outline = IO.replace(line, self.replacements)
                        outfile.write(outline)
                        
                ## Store previous file
                self.replacements['FILE_N'] = 'md_{}_{}'.format(l1,
                                                                l2)
                
                self.mdfiles.append('md_{}_{}'.format(l1,l2))


    def write_runfile(self):
        ntasks = getattr(s, self.cluster)['NTASKS']
        src = s.INPUT_DIR + '/run.sh'
        tgt = self.directory + '/inputfiles/run' + self.cluster + '.sh'
        EQ_files = sorted(glob.glob(self.directory + '/inputfiles/eq*.inp'))
        
        if self.start == '1':
            MD_files = reversed(sorted(glob.glob(self.directory + '/inputfiles/md*.inp')))
            
        elif self.start == '0.5':
            MD_files = self.mdfiles

        replacements = IO.merge_two_dicts(self.replacements, getattr(s, self.cluster))
        replacements['FEPS'] = ' '.join(self.FEPlist)
            
        with open(src) as infile, open(tgt, 'w') as outfile:
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
                        outline = 'time mpirun -np {} $qdyn {}.inp' \
                                   ' > {}.log\n'.format(ntasks,
                                                       file_base,
                                                       file_base)
                        outfile.write(outline)
                        
                if line.strip() == '#RUN_FILES':
                    for line in MD_files:
                        if self.start == '1':
                            file_base = line.split('/')[-1][:-4]
                            
                        elif self.start == '0.5':
                            file_base=line
                            
                        outline = 'time mpirun -np {} $qdyn {}.inp'  \
                                  ' > {}.log\n'.format(ntasks,
                                                      file_base,
                                                      file_base)
                        outfile.write(outline)
                        
    def write_submitfile(self):
        IO.write_submitfile(self.directory, self.replacements)
    
    def write_FEPfile(self):
        if self.dual == True:
            run.write_dualFEPfile()
            
        if self.dual == False:
            run.write_singleFEPfile()
        
    def write_singleFEPfile(self):
        RES_list = self.PDB[int(self.PDB2Q[self.mutation[1]])]
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
            RES_dic[line[2].strip()] = line

        
        for src in sorted(glob.glob(self.FEPdir + '/FEP*.fep')):
            tgt = self.directory + '/inputfiles/' + src.split('/')[-1]
            with open (src) as infile, open (tgt, 'w') as outfile:
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
        parameterfile = s.FF_DIR + '/' + self.forcefield + '.prm'
        vdw_prms = IO.read_prm([parameterfile])['[atom_types]']
        
        ## Next, we start to construct the FEP file
        with open(self.directory + '/inputfiles/FEP1.fep', 'w') as outfile:
            # First create the file header
            outfile.write('! Dual topology QresFEP: ' + self.MUTresn + '\n\n')
            outfile.write('[FEP]\n')
            outfile.write('states 2\n')
            outfile.write('softcore_use_max_potential on \n\n')
            
            # Write out [atoms] section
            outfile.write('[atoms]\n')
            cnt = 0
            for line in self.PDB[int(self.mutation[1])]:
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
                        outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                     line[3], 
                                                                     line[3]))
                        
                    else:
                        outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                     line[3], 
                                                                     '0.0000'))
                    
                else:
                    outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                 '0.0000', 
                                                                 line[3]))
                    
            outfile.write('\n')
            
            # Write out [atom_types] section
            outfile.write('[atom_types]\n')
            for line in vdw_prms:
                line = line.split()
                if len(line) == 8:
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
            # Write out DUM entry        
            outfile.write('{:<8s}{:>10s}{:>10s}{:>10s}'\
                          '{:>10s}{:>10s}{:>10s}{:>10s}\n'.format('DUM',
                                                                  '0.0000',
                                                                  '0.0000',
                                                                  '0.0000',
                                                                  '0.0000',
                                                                  '0.0000',
                                                                  '0.0000',
                                                                  '1.0080',
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
                        outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                     line[2], 
                                                                     line[2]))
                    else:
                        outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                     line[2], 
                                                                     'DUM'))
                        
                else:
                    outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                 'DUM', 
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
                        outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                     '0', 
                                                                     '20'))
                    
                else:
                    outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                 '20', 
                                                                 '0'))
            
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
        
        with open(qfep_in) as infile, open(qfep_out, 'w') as outfile:
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
                        choices = ['OPLS2015', 'OPLS2005', 'SIDECHAIN'],
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
                    default = 'LOCAL',
                    help = "define this variable if you are setting up your system elsewhere")
    
    parser.add_argument('-d', '--dual',
                        dest = "dual",
                        default = False,
                        action = 'store_true',
                        help = "Turn on if you want QresFEP to generate a dual topology hybrid")
    
    parser.add_argument('-l', '--start',
                        dest = "start",
                        default = '1',
                        choices = ['1', '0.5'],
                        help = "Starting FEP in the middle or endpoint, 0.5 recommended for dual"
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
              include = ('ATOM','HETATM')
             )
    
    run.checkFEP()                      # 00
    run.create_environment()            # 01
    run.read_input()                    # 02    
    run.readpdb()                       # 03
    run.create_dual_lib()               # 04  
    run.merge_prm()                     # 05
    run.write_pdb()                     # 06
    run.select_waters()                 # 07
    run.write_qprep()                   # 08
    run.run_qprep()                     # 09
    run.get_lambdas()                   # 10
    run.write_EQ()                      # 11
    run.write_MD()                      # 12
    run.write_submitfile()              # 13
    run.write_runfile()                 # 14
    run.write_FEPfile()                 # 15
    run.write_qfep()                    # 16
