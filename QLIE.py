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
    def __init__(self, ligand, cofactor, forcefield, include, system,
                 preplocation, cluster, temperature, replicates,
                 *args, **kwargs):
        # Argparse arguments
        self.ligand         = ligand
        self.cofactor       = [ligand]
        self.forcefield     = forcefield
        self.system         = system
        self.preplocation   = preplocation
        self.cluster        = cluster
        if temperature == None:
            self.temperature = s.TEMPERATURE
            
        else:
            self.temperature = temperature
            
        if replicates == None:
            self.replicates = s.REPLICATES
            
        else:
            self.replicates = replicates
            
        # Constructs for local variables
        self.include        = include
        self.replacements   = {'TEMP_VAR':self.temperature,
                               'RUN_VAR':self.replicates,
                               'RUNFILE':'run' + self.cluster + '.sh'  
                              }
        self.CYX            = []
        self.PDB2Q          = {}
        self.PDB            = {}
        self.systemsize     = 0

        
        # Add cofactors to list for further pdb/prm parsing
        if cofactor != None:
            self.cofactor.append(cofactor)
    
    def create_environment(self):
        self.directory = 'LIE_{}'.format(self.ligand)
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

    def read_pdb(self):
        atnr = 0
        if self.system == 'protein':
            pdb_files = ['protein.pdb']
        else:
            pdb_files = []
        
        if self.cofactor[0] != None:
            for line in self.cofactor:
                pdb_files.append(line + '.pdb')
        
        for pdb_file in pdb_files:
            with open(pdb_file) as infile:
                for line in infile:
                    if line.startswith(self.include):
                        atnr += 1
                        line = IO.pdb_parse_in(line)
                        if pdb_file != 'protein.pdb' and self.system == 'protein':
                                line[6] = resoffset + 1
                                line[1] = self.systemsize

                        line[1] = atnr
                        
                        # Construct the PDB dictionary            
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
    def write_EQ(self):
        # If water or vacuum system, apply restrain on ligand
        if self.system != 'protein':
            at1 = self.PDB[len(self.PDB)][0][1]
            at2 = self.PDB[len(self.PDB)][-1][1]
            self.replacements['WATER_RESTRAINT'] = '{} {} 1.0 0 0'.format(at1, 
                                                                          at2)
        
        else:
            self.replacements['WATER_RESTRAINT'] = ''
            
        self.replacements['ATOM_END'] = '{}'.format(self.systemsize)
        self.replacements['ATOM_START_LIG1'] = '{}'.format(1)
        self.replacements['EQ_LAMBDA'] = '1.000 0.000'
        
        for EQ_file in glob.glob(s.INPUT_DIR + '/eq*.inp'):
            src = EQ_file
            EQ_file = EQ_file.split('/')
            tgt = self.directory + '/inputfiles/' + EQ_file[-1]
            with open(src) as infile, open(tgt, 'w') as outfile:
                for line in infile:
                    outline = IO.replace(line, self.replacements)
                    outfile.write(outline)
                    
    def write_MD(self):
        self.replacements['FILE_N'] = 'eq5'
        src = s.INPUT_DIR + '/md_LIE_XX.inp'
        for i in range(1, 11):
            tgt = self.directory + '/inputfiles/md_LIE_{:02d}.inp'.format(i)
            self.replacements['FILE'] = 'md_LIE_{:02d}'.format(i)
            
            with open(src) as infile, open(tgt, 'w') as outfile:
                for line in infile:
                    outline = IO.replace(line, self.replacements)
                    outfile.write(outline)                
        
            self.replacements['FILE_N'] = 'md_LIE_{:02d}'.format(i)
        
    def write_runfile(self):
        ntasks = getattr(s, self.cluster)['NTASKS']
        src = s.INPUT_DIR + '/run.sh'
        tgt = self.directory + '/inputfiles/run' + self.cluster + '.sh'
        EQ_files = sorted(glob.glob(self.directory + '/inputfiles/eq*.inp'))

        MD_files = sorted(glob.glob(self.directory + '/inputfiles/md*.inp'))

        replacements = IO.merge_two_dicts(self.replacements, getattr(s, self.cluster))
        replacements['FEPS'] = 'FEP1.fep'
            
        with open(src) as infile, open(tgt, 'w') as outfile:
            for line in infile:
                if line.strip() == '#SBATCH -A ACCOUNT':
                    try:
                        replacements['ACCOUNT']
                        
                    except:
                        line = ''
                        
                outline = IO.replace(line, replacements)
                if line[0:7] != 'timeout':
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
                                  ' > {}.log\n'.format(ntasks,
                                                      file_base,
                                                      file_base)
                            
                        outfile.write(outline)     

    def write_submitfile(self):
        IO.write_submitfile(self.directory, self.replacements)                        
                        
    def write_FEPfile(self):
        vdw_prms = IO.read_prm([self.ligand + '.prm'])['[atom_types]']
        with open(self.directory + '/inputfiles/FEP1.fep', 'w') as outfile:
            # First create the file header
            outfile.write('! LIE ' + self.ligand + '\n\n')
            outfile.write('[FEP]\n')
            outfile.write('states 2\n')
            # Write out [atoms] section
            outfile.write('[atoms]\n')
            cnt = 0
            for line in self.PDB[len(self.PDB)]:
                cnt += 1
                outfile.write('{:4d} {:4d} ! {:<4s}\n'.format(cnt, line[1], line[2]))
            outfile.write('\n')

            outfile.write('[atom_types]\n')
            for line in vdw_prms:
                line = line.split()
                if len(line) == 7:
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

            outfile.write('[change_atoms]\n')
            cnt = 0
            for line in vdw_prms:
                line = line.split()
                if len(line) == 7:
                    cnt += 1
                    outfile.write('{:4d}{:>10s}{:>10s}\n'.format(cnt, 
                                                                 line[0], 
                                                                 line[0]))
            
            outfile.write('\n')
            
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='QLIE',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = ' == Generate LIE inputfiles for simulation and analysis == ')

    parser.add_argument('-l', '--ligand',
                        dest = "ligand",
                        required = True,
                        help = "Name of the ligand to perform LIE calculations on")
    
    parser.add_argument('-c', '--cofactor',
                        dest = "cofactor",
                        required = False,
                        help = "PDB files of cofactors other than the ligand with *.prm and *.lib files")
    
    parser.add_argument('-f', '--forcefield',
                        dest = "forcefield",
                        required = True,
                        choices = ['OPLS2015', 'OPLS2005'],
                        help = "Forcefield to use.")
    
    parser.add_argument('-S', '--system',
                        dest = "system",
                        required = True,
                        choices = ['protein', 'water', 'vacuum'],
                        help = "System type, can be protein, water or vacuum")
    
    parser.add_argument('-P', '--preplocation',
                    dest = "preplocation",
                    default = 'LOCAL',
                    help = "define this variable if you are setting up your system elsewhere")
    
    parser.add_argument('-C', '--cluster',
                        dest = "cluster",
                        required = True,
                        help = "Cluster to use.")    
    
    parser.add_argument('-T', '--temperature',
                        dest = "temperature",
                        required = False,
                        help = "List of temperatures, default taken from settings.py")
    
    parser.add_argument('-r', '--replicates',
                        dest = "replicates",
                        required = False,
                        help = "Amount of replicates, default taken from settings.py")     
    
    args = parser.parse_args()
    run = Run(ligand        = args.ligand,
              cofactor      = args.cofactor,
              forcefield    = args.forcefield,
              system        = args.system,
              preplocation  = args.preplocation,
              cluster       = args.cluster,
              temperature   = args.temperature,
              replicates    = args.replicates,
              include       = ('ATOM', 'HETATM')
             )
    
    run.create_environment()        # 00
    run.read_input()                # 01
    run.read_pdb()                  # 02
    run.merge_prm()                 # 03
    run.write_pdb()                 # 04
    run.select_waters()             # 05
    run.write_qprep()               # 06
    run.run_qprep()                 # 07
    run.write_EQ()                  # 08
    run.write_MD()                  # 09
    run.write_runfile()             # 10
    run.write_submitfile()          # 11
    run.write_FEPfile()             # 12
