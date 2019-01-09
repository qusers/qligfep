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
    Setup residue FEPs using either a single or dual topology approach (latter under
    construction).
    """
    def __init__(self, cofactor, mutation, include, forcefield, windows,
                 sampling, system, cluster, temperature, replicates, 
                 preplocation, *args, **kwargs):
        
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
                print "ERROR: {} are required files, exiting now.".format(required)
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
        self.systemsize = 1
        self.system = system
        self.cluster = cluster
        self.FEPlist = []
        
        self.replacements = {'EQ_LAMBDA': '1.000 0.000',
                     'ATOM_START_LIG1':'1',
                     'WATER_RESTRAINT':'',
                     'TEMP_VAR':self.temperature,
                     'RUN_VAR':self.replicates,
                     'RUNFILE':'run' + self.cluster + '.sh'  
                    }
    
    def checkFEP(self):
        if len(self.mutation[0]) == 1:
            AA_from = IO.AA(self.mutation[0])
            
        else:
            AA_from = self.mutation[0]
            
        AA_to = IO.AA(self.mutation[2])
        mutation = '{}{}{}'.format(*self.mutation)
        FEPdir = s.FF_DIR + '/.FEP/' + AA_from + '_' + AA_to
        if self.forcefield == 'OPLS2005':
            FEPdir = FEPdir + '-aa'
            
        if not os.path.exists(FEPdir):
            print 'FATAL: no FEP files found for the {} mutation' \
                  'in {} exiting now.'.format(mutation,
                                              FEPdir)
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
        pdb_files = ['protein.pdb']
        if self.cofactor[0] != None:
            for line in self.cofactor:
                pdb_files.append(line + '.pdb')
        
        for pdb_file in pdb_files:
            with open(pdb_file) as infile:
                for line in infile:
                    if line.startswith(self.include):
                        line = IO.pdb_parse_in(line)
                        if pdb_file != 'protein.pdb':
                            line[6] = resoffset + 1
                            line[1] = self.systemsize
                        try:
                            self.PDB[line[6]].append(line)

                        except:
                            self.PDB[line[6]] = [line]
                        
                        self.systemsize += 1
                
                resoffset = len(self.PDB)

                       
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
                for line in self.PDB[key]:
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
                        line = IO.pdb_parse_in(line)    
                        cofactor_coordinates.append([line[8], line[9], line[10]])
        
        with open (src) as infile, open (tgt, 'w') as outfile:
            for line in infile:
                if line.split()[-1] == 'SPHERE':
                    outfile.write(line)
                    continue
                    
                else:
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
                        'SOLVENT':'4 water.pdb'
                       }
        src = s.INPUT_DIR + '/qprep_resFEP.inp'
        self.qprep = self.directory + '/inputfiles/qprep.inp'
        libraries = [self.forcefield + '.lib']
        if self.cofactor[0] != None:
            for filename in self.cofactor:
                libraries.append(filename + '.lib')
                
        with open(src) as infile, open(self.qprep, 'w') as outfile:
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
        for line in self.PDB[1]:
            if line[2] == 'CA' and self.system == 'water' or self.system == 'vacuum':
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
            
    def write_runfile(self):
        ntasks = getattr(s, self.cluster)['NTASKS']
        src = s.INPUT_DIR + '/run.sh'
        tgt = self.directory + '/inputfiles/run' + self.cluster + '.sh'
        EQ_files = sorted(glob.glob(self.directory + '/inputfiles/eq*.inp'))
        MD_files = reversed(sorted(glob.glob(self.directory + '/inputfiles/md*.inp')))
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
                        file_base = line.split('/')[-1][:-4]
                        outline = 'time mpirun -np {} $qdyn {}.inp'  \
                                  '> {}.log\n'.format(ntasks,
                                                      file_base,
                                                      file_base)
                        outfile.write(outline)
                        
    def write_submitfile(self):
        IO.write_submitfile(self.directory, self.replacements)
                        
    def write_singletop_FEPs(self):
        # This piece needs to be cleaner at some point
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

    def write_qfep(self, inputdir, windows, lambdas):
        qfep_in = s.ROOT_DIR + '/INPUTS/qfep.inp'
        qfep_out = writedir + '/inputfiles/qfep.inp'
        i = 0
        total_l = len(lambdas)

        # TO DO: multiple files will need to be written out for temperature range
        kT = f.kT(float(self.temperature))
        replacements = {}
        replacements['kT']=kT
        replacements['WINDOWS']=windows
        replacements['TOTAL_L']=str(total_l)
        with open(qfep_in) as infile, open(qfep_out, 'w') as outfile:
            for line in infile:
                line = run.replace(line, replacements)
                outfile.write(line)

            if line == '!ENERGY_FILES\n':
                for i in range(0, total_l):
                    j = -(i + 1)
                    lambda1 = lambdas[i]
                    lambda2 = lambdas[j]
                    filename = 'md_' +                          \
                                lambda1.replace('.', '') +      \
                                '_' +                           \
                                lambda2.replace('.', '') +      \
                                '.en\n'

                    outfile.write(filename)
                            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='protPREP',
        version='1.0',
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
                        choices = ['OPLS2015', 'OPLS2005'],
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
              include = ('ATOM','HETATM')
             )
    
    run.checkFEP()
    run.create_environment()
    run.readpdb()
    run.read_input()
    run.merge_prm()
    run.write_pdb()
    run.select_waters()
    run.write_qprep()
    run.run_qprep()
    run.get_lambdas()
    run.write_EQ()
    run.write_MD()
    run.write_submitfile()
    run.write_runfile()
    run.write_singletop_FEPs()
