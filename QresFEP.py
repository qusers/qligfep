#!/usr/bin/env python

import argparse
import re
import os
import shutil
import sys
import glob
import subprocess
import itertools

import functions as f
import settings as s
import IO

from scripts import counter_ions

class Run(object):
    """
    Setup residue FEPs using either a single or dual topology approach
    """
    def __init__(self, mutation, mutchain, system, tripeptide, dual, cofactors,
                 forcefield, windows, sampling, start, timestep, temperature,
                 replicates, cluster, preplocation, *args, **kwargs):
        
        self.mutation = re.split('(\d+)', mutation)
        self.chain = mutchain
        self.system = system
        self.tpt = tripeptide
        self.dual = dual
        self.cofactors = cofactors
        self.forcefield = forcefield
        self.sampling = sampling
        self.start = start
        self.windows = int(windows)
        self.timestep = timestep
        self.temperature = temperature
        self.replicates = replicates
        self.cluster = cluster
        self.preplocation = preplocation

        self.lambdas = []
        self.CYX = []
        self.PDB = {}
        self.PDB2Q = {}
        self.systemsize = 0
        self.FEPlist = []
        self.nonAA = False

        # Check whether all required files are there:
        required = ['protein.pdb', 'water.pdb', 'protPREP.log']
        if self.cofactors:
            extensions = ['.pdb', '.lib', '.prm']
            for cofactor in self.cofactors:
                for extension in extensions:
                    required.append(cofactor + extension)
        for filename in required:
            if os.path.exists(filename) == False:
                required = ' '.join(required)
                print(f"ERROR: {required} are required files, exiting now.")
                sys.exit()

        # Create a dict for all possible side-chain atom names in the mutation
        self.atoms = {atom: 0 for atom in IO.atoms} # Add all wild-type atoms
        self.atoms.update({atom.lower(): 0 for atom in IO.atoms}) # mutant atoms
        
        self.replacements = {
            'ATOM_START_LIG1': '1',
            'WATER_RESTRAINT': '',
            'TEMP_VAR': self.temperature,
            'RUN_VAR': self.replicates,
            'RUNFILE': f'run{self.cluster}.sh'}
        
        if self.start == '1' or self.start == '0':
            self.replacements['EQ_LAMBDA'] = '1.000 0.000'
            
        if self.start == '0.5':
            if self.dual:
                self.replacements['EQ_LAMBDA'] = '0.000 1.000'
            else:
                self.replacements['EQ_LAMBDA'] = '0.500 0.500'
    
    def checkFEP(self):
        # Make sure the given mutation will have 3 letter-code amino acid names
        if len(self.mutation[0]) == 1:
            self.mutation[0] = IO.AA(self.mutation[0])
        try:
            if len(self.mutation[2]) == 1:
                self.mutation[2] = IO.AA(self.mutation[2])
        except: # Non-natural amino acid
            print("residue not found in library, assuming non-natural AA")
            self.mutation[2] = 'X'
            self.nonAA = True
        
        mutation = '{}{}{}'.format(*self.mutation)
        FEPdir = f"{s.FF_DIR}/.FEP/{self.forcefield}/{self.mutation[0]}_{self.mutation[2]}"

        # Check wether the mutation involves Glycine    
        self.FromGly, self.ToGly = False, False
        if self.mutation[0] == 'G' or self.mutation[0] == 'GLY':
            self.FromGly = True
            self.replacements['EQ_LAMBDA'] = '0.000 1.000'
        if self.mutation[2] == 'G' or self.mutation[2] == 'GLY':
            self.ToGly = True
        
        # FEP input files will be generated for dual topology mutations
        if self.dual == True:
            print('Generating dual topology inputfiles')
            self.FEPlist = ['FEP1.fep', 'FEP2.fep']
            self.FEPdir = None
        
        elif not os.path.exists(FEPdir):
            print(f'FATAL: no FEP files found for the {mutation} mutation in {FEPdir} exiting now.')
            sys.exit()

        else: # FEP input files will be taken from the template location (qligfep/.FEP/)
            print('Generating single topology inputfiles')
            for fep in sorted(glob.glob(f"{FEPdir}/FEP*.fep")):
                self.FEPlist.append(fep.split('/')[-1])
            self.FEPdir = FEPdir
    
    # Create FEP directory and copy template files
    def create_environment(self):
        self.directory = 'FEP_{}{}{}'.format(*self.mutation)
        print(self.directory)
        # Create FEP dir
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
            os.makedirs(self.directory + '/inputfiles')
            
        # Copy force field and cofactor library files to FEP dir
        shutil.copy(f"{s.FF_DIR}/{self.forcefield}.lib", f"{self.directory}/inputfiles/{self.forcefield}.lib")
        if self.cofactors:
            for cofactor in self.cofactors:
                shutil.copy(f"{cofactor}.lib", f"{self.directory}/inputfiles/{cofactor}.lib")
    
    # Extract protprep.py information from protPREP.log
    def read_input(self):
        block = 0
        with open('protPREP.log') as infile:
            for line in infile:
                line = line.split()
                if len(line) > 1:
                    if line[1] == 'center:':
                        self.sphere = [float(coord) for coord in line[2:]]
                        
                    if line[1] == 'radius:':
                        self.radius = line[2]
                        self.replacements['SPHERE'] = self.radius
                    
                    if line[1] == 'charge':
                        self.charge = int(line[4])
                        
                    if line[0] == 'Q_CYS1':
                        block = 1
                        
                    if line[0] == 'pdbfile':
                        block = 2
                        
                    if line[0][0] == '-':
                        block = 0
                        
                    if block == 1:
                        if line[0].isdigit():
                            self.CYX.append([line[0], line[1]])

                    if block == 2 and (len(line) == 3 or len(line) == 4):
                        chain = ' ' if len(line) == 3 else line[2]
                        self.PDB2Q.setdefault(chain, {})[line[1]] = line[0]

    # Read protprep.py generated protein.pdb
    def readpdb(self):
        # Determine the hybrid residue name
        if self.nonAA != True:
            self.MUTresn = f"{IO.AA(self.mutation[0])}2{IO.AA(self.mutation[2])}"
        else:
            self.MUTresn = f"{self.mutation[0]}2X"

        # Determine backbone atoms
        if self.ToGly == True or self.FromGly == True:
            self.backbone = ['C', 'O', 'CA', 'N', 'H']
        else: 
            self.backbone = ['C', 'O', 'CA', 'N', 'H', 'HA']

        # Read in the mutant residue
        if self.dual == True:
            mutation = []
            mut_pdb = f"{self.mutation[2]}{self.mutation[1]}.pdb"
            with open(mut_pdb) as infile:   # Reads the mutant residue pdb as MUTX.pdb
                for at in infile:
                    if at.startswith('ATOM'):
                        at = IO.pdb_parse_in(at)

                        atn = at[2]
                        resi = at[6] = int(self.PDB2Q[self.chain][self.mutation[1]])
                        
                        if atn in self.backbone:
                            continue
                        else:
                            mutation.append(at)
        
        pdb_files = ['protein.pdb']
        if self.cofactors:
            for cofactor in self.cofactors:
                pdb_files.append(f'{cofactor}.pdb')

        id = itertools.count(1)
        for pdb_file in pdb_files:
            with open(pdb_file) as infile:
                for at in infile:
                    if at.startswith(('ATOM', 'HETATM')):
                        at = IO.pdb_parse_in(at)

                        if pdb_file != 'protein.pdb':
                            at[6] = len(self.PDB) + 1
                            
                        at[1] = next(id)
                        atn  = at[2]
                        resi = str(at[6])
                        
                        # Change to hybrid residue name
                        if self.dual == True:
                            if resi == self.PDB2Q[self.chain][self.mutation[1]]:
                                at[4] = self.MUTresn
                        
                        # Construct the PDB dictionary            
                        self.PDB.setdefault(int(resi), []).append(at)
                        
                        if self.dual == True:
                            if resi == self.PDB2Q[self.chain][self.mutation[1]]:
                                if atn == 'O':
                                    # Add mutant side chain atoms to hybrid residue
                                    for at_mut in mutation:
                                        if at_mut[2] not in self.backbone:

                                            at_mut[1] = next(id)
                                            at_mut[2] = at_mut[2].lower()
                                            at_mut[4] = self.MUTresn
                                            
                                            # Generate the dual topology residue in the PDB dictionary
                                            self.PDB.setdefault(at_mut[6], []).append(at_mut)

                                            self.systemsize += 1
                            
                        self.systemsize += 1
                    else: # Just adds TER lines
                        self.PDB[int(at[23:26])].append(at)
    
    # Function to generate a library file for the hybrid residue if using dual topology
    def create_dual_lib(self):
        if self.dual == False:
            return None

        self.merged_lib = {} # Dictionairy for storing hybrid residue .lib file
        FFlib = {} # Dictionairy for storing residues from FF.lib file
        replacements = {}   
        
        # Headers for different sections in the .lib file
        headers =['[atoms]', '[bonds]', '[connections]', '[impropers]', '[info]']
        if self.forcefield.startswith('OPLS'):
            headers.append('[charge_groups]')
        
        # The FF.lib file (and possible the non-natural amino acid .lib file)
        libfiles = [f"{s.FF_DIR}/{self.forcefield}.lib"]
        if self.nonAA == True:
            libfiles.append(f"{self.mutation[2]}.lib")
            
        for libfile in libfiles:
            with open(libfile) as infile:
                for line in infile:
                    
                    if line.startswith('*'): # Skip comments
                        continue

                    if line.startswith('{'): # Find residues
                        RESname = line[1:line.index('}')]
                        continue
                    
                    # Only store residues involved in the mutation
                    if RESname not in [self.mutation[0], self.mutation[2]]:
                        continue
                    else:
                        FFlib.setdefault(RESname, []).append(line)
        
        header = None
        # Read the WT residue      
        for line in FFlib[self.mutation[0]]:
            if 'C223' in line:  #Change atom type of GLY 
                line = line.replace('C223', 'C224')
            # Add lines of each section for the different headers in th .lib file
            if any(line.lstrip().startswith(header) for header in headers):
                header = line.strip().split()[0]
                self.merged_lib[header] = []
            else:
                self.merged_lib[header].append(line)

        # Add a bond line to connect the main chain CA to the mutant side chain
        if self.FromGly == True:
            self.merged_lib['[bonds]'].append('\tCA    ha\n')
        if self.ToGly == True:
            self.merged_lib['[bonds]'].append('\tCA    ha2\n')
            self.merged_lib['[bonds]'].append('\tCA    ha3\n')
        else:
            self.merged_lib['[bonds]'].append('\tCA    cb\n')
        
        # Read the mutant residue
        for line in FFlib[self.mutation[2]]:
            if any(line.lstrip().startswith(header) for header in headers):
                header = line.strip().split()[0]
                continue
            if header == '[atoms]':
                atn = line.strip().split()[1]
                if atn not in self.backbone:
                    replacements[atn] = atn.lower() # Set non-backbone mutant atoms to lowercase

            if header == '[connections]': # connection flag only needs to be there once
                continue
            # Remove overlapping backbone definitions
            if line.strip().split()[0] in self.backbone or line.strip().split()[1] in self.backbone:
                continue

            # Merge the library on the reference
            line = IO.replace(line, replacements)
            self.merged_lib[header].append(line)
        
        # Write out the merged library file
        hyb_lib = self.directory + '/inputfiles/' + self.MUTresn + '.lib'
        with open(hyb_lib, 'w') as lib:
            self.merged_lib['[charge_groups]'] = [] # No charge groups for Q atoms required

            lib.write(f'{{{self.MUTresn}}}\n')
            for header in headers:
                lib.write(f'\n    {header}\n')
                if not header in self.merged_lib:
                    continue
                for id, line in enumerate(self.merged_lib[header], 1):
                    if header == '[atoms]':
                        line = line.split()
                        lib.write(f"    {id:2d} {line[1]:<10}{line[2]:<10}{line[3]:<10}\n")
                    else:
                        lib.write(line)
                        

    # Function to obtain the missing bonded parameter terms for the hybrid residue
    def get_zeroforcebonded(self):
        if self.dual == False:
            return None    

        # Dictionary for .prm file terms
        self.zerofk = {'[angles]':[], '[bonds]':[], '[impropers]':[], '[torsions]':[]}
        prmfiles = [s.FF_DIR + '/' + self.forcefield + '.prm']
        if self.nonAA == True:
            prmfiles.append(f"{self.mutation[2]}.prm")
        hyb_lib = self.directory + '/inputfiles/' + self.MUTresn + '.lib'
        headers =['[options]', '[atom_types]', '[bonds]', '[angles]', '[torsions]', '[impropers]']
    
        # Write tmp qprep.inp and prm file
        prms = IO.read_prm(prmfiles)
        prm_tmp = f"{self.directory}/inputfiles/tmp.prm"
        qprep_tmp = f"{self.directory}/inputfiles/tmp.inp"
        pdb_tmp = f"{self.directory}/inputfiles/tmp.pdb"
        
        with open(pdb_tmp, 'w') as pdb: # Write out tmp.pdb file
            for resi, atoms in self.PDB.items():
                for atom in atoms:
                    try:
                        pdb.write(f"{IO.pdb_parse_out(atom)}\n")
                    except:
                        pdb.write(atom)
        
        with open(prm_tmp, 'w') as prm: # Write out tmp.prm file
            for header in headers:
                prm.write(f"{header}\n")
                for line in prms[header]:
                    prm.write(line)
                    
        with open(qprep_tmp, 'w') as qprep:
            qprep.write(f"rl {hyb_lib.split('/')[-1]}\n")
            qprep.write(f"rl {self.forcefield}.lib\n")
            qprep.write(f"rprm {prm_tmp.split('/')[-1]}\n")
            qprep.write(f"rp {pdb_tmp.split('/')[-1]}\n")
            qprep.write(f"boundary 1 {self.sphere[0]} {self.sphere[1]} {self.sphere[2]} {self.radius}\n")
            qprep.write(f"maketop tmp.top\n")
            qprep.write("q")
            
        # Run qprep to find missing parameters to be set to zero
        os.chdir(self.directory + '/inputfiles/')
        q_commands = getattr(s, self.preplocation)
        qprep = q_commands['QPREP']
        os.system(f'{qprep} < tmp.inp > tmp.out')
        
        # Read qprep tmp.out to retrieve the missing parameters
        with open('tmp.out') as tmp_out:
            for line in tmp_out:
                line = line.split()
                if len(line) < 3:
                    continue
                    
                if line[1] != 'Missing':
                    continue
                
                # Create zero fk bonds
                if line[2] == 'bond':
                    bond_prm = f"{line[-2]:11}{line[-1]:11}{0.0: 8.2f}{0.0:>12.7}\n"
                    if bond_prm not in self.zerofk['[bonds]']:
                        self.zerofk['[bonds]'].append(bond_prm)
                
                # Create zero fk angles
                if line[2] == 'angle':
                    angle_prm = f"{line[-3]:11}{line[-2]:11}{line[-1]:11}{0.0: 8.2f}{110.0:>12.7}\n"
                    if angle_prm not in self.zerofk['[angles]']:
                        self.zerofk['[angles]'].append(angle_prm)
                    
                # Create zero fk torsions
                if line[2] == 'torsion':
                    torsion_prm = f"{line[-4]:11}{line[-3]:11}{line[-2]:11}{line[-1]:11}{0.0:<10.3f}{1:2d}.000{'0.000':>10}{'1':>10}\n"
                    if torsion_prm not in self.zerofk['[torsions]']:
                        self.zerofk['[torsions]'].append(torsion_prm)
        
        # Clean up tmp files
        for tmp in glob.glob('tmp*'):
            os.remove(tmp)
        os.chdir('../../')
    
    # Function to generate a parameter file for the hybrid residue if using dual topology
    def merge_prm(self):
        headers =['[options]', '[atom_types]', '[bonds]', '[angles]', '[torsions]', '[impropers]']
        prmfiles = [f"{s.FF_DIR}/{self.forcefield}.prm"]
        if self.cofactors:
            for cofactor in self.cofactors:
                prmfiles.append(f"{cofactor}.prm")
                
        if self.nonAA == True:
            prmfiles.append(f"{self.mutation[2]}.prm")
            
        prms = IO.read_prm(prmfiles)
        prm_merged = f"{self.directory}/inputfiles/{self.forcefield}_merged.prm"
        self.prm_merged = f"{self.forcefield}_merged.prm"

        with open (prm_merged, 'w') as hyb_prm: # Write out hybrid parameter file
            for header in headers:
                hyb_prm.write(f"{header}\n")
                for line in prms[header]:
                    hyb_prm.write(line)
                if self.dual == True:
                    if header in self.zerofk:
                        hyb_prm.write(f'! Zero order {header.strip("[]")} dual topology\n')
                        for zerofk in self.zerofk[header]:
                            hyb_prm.write(zerofk)
                        hyb_prm.write('\n')

    # Function to generate a pdb inputfile for qprep
    def write_pdb(self):
        pdb_file = f"{self.directory}/inputfiles/complex.pdb"
        if self.system == 'protein': # Write out the PDB dict to complex.pdb
            with open(pdb_file, 'w') as pdb:
                for atom in [residue for residues in self.PDB.values() for residue in residues]:
                    try:
                        pdb.write(f"{IO.pdb_parse_out(atom)}\n")
                    except:
                        pdb.write(atom)
                        
        elif self.system == 'water': # Write out only the relevant part for the reference tripeptide to complex.pdb
            ter = 3
            with open(pdb_file, 'w') as pdb: 
                if self.tpt == 'Z': # Write out also the flaning residues of the mutable residue for natural sequence tripeptide
                    tpt_resi = [int(self.PDB2Q[self.chain][self.mutation[1]]) + i for i in range(-1, 2)]
                    resi = [2,3,4]
                else: # Only write out the mutable residue for no, Ala, or Gly flanking residues
                    tpt_resi = [int(self.PDB2Q[self.chain][self.mutation[1]])]
                    resi = [3]
                id = itertools.count(1)
                for _, res in enumerate(tpt_resi):
                    for atom in self.PDB[res]:
                        atom[5], atom[6] = '', resi[_] # Remove chain and reset residue number
                        if resi[_] == 4 and atom[2] == 'CA': # Remember CA id of the post-mutable residue as anchorpoint
                            ter = next(id)
                        else:
                            next(id)
                        pdb.write(f"{IO.pdb_parse_out(atom)}\n")
            
            # Replacements for the template PyMOL file to match tripeptide of use
            replacements = {'PDB': 'complex.pdb',
                            'COMPLEX': 'complex',
                            'TER': str(ter),
                            'TRIPEPTIDE': 'tripeptide.pdb'}
            ala_or_gly = {'A': 'ala', 'G': 'gly'}
            if self.tpt in ala_or_gly:
                aog = ala_or_gly[self.tpt]
                replacements['PRE'] = f"editor.attach_amino_acid('pk1', '{aog}')"
                replacements['POST'] = f"cmd.edit('(COMPLEX`1)',None,None,None,pkresi=0,pkbond=0)\neditor.attach_amino_acid('pk1', '{aog}')"

            # Update the template PyMOL tripeptide.pml script
            with open(f"{s.INPUT_DIR}/tripeptide.pml") as pml_tmplt, \
                 open(f"{self.directory}/inputfiles/tripeptide.pml", 'w') as pml_trgt:
                for line in pml_tmplt:
                    if self.tpt in ['X', 'Z'] and any(x in line for x in ('PRE', 'POST')):
                        continue
                    pml_trgt.write(IO.replace(line, replacements))
            
            # Execute PyMOL to generate a reference tripeptide.pdb
            cwd = self.directory + '/inputfiles'
            pymol = ['pymol', '-c', 'tripeptide.pml']
            result = subprocess.run(pymol, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if result.returncode == 0:
                print("PyMOL executed successfully")
            else:
                print("PymMOL error...")

            # Replacements to update atom and residue names of PyMOL output tripeptide.pdb to match .lib file 
            replace =  [('1HB ', ' HB1'),   ('2HB ', ' HB2'),   ('3HB ', ' HB3'),
                        ('1HH3 ACE', 'HH31 ACE'), ('2HH3 ACE', 'HH32 ACE'), ('3HH3 ACE', 'HH33 ACE'),
                        (' N   NME', ' N   NMA'), (' H   NME', ' H   NMA'),
                        (' CH3 NME', ' CA  NMA'),
                        ('1HH3 NME', ' HA1 NMA'), ('2HH3 NME', ' HA2 NMA'), ('3HH3 NME', ' HA3 NMA'),
                        ('3HA  GLY', ' HA3 GLY'), ('HA  GLY', 'HA2 GLY')]
            
            id = itertools.count(1)
            with open(f"{self.directory}/inputfiles/tripeptide.pdb", 'r') as tpt, \
                open(pdb_file, 'w') as pdb:
                    for atom in tpt:
                        if atom.startswith('ATOM'):
                            for old, new in replace: # replace the atom and residue names
                                atom = re.sub(old, new, atom)
                            pdb.write(atom)
                            next(id)

                    charge = self.charge # Subtract the charge of the mutable residue from the total
                    if any(res in ['ASP', 'GLU', 'LYS', 'ARG', 'HIP'] for res in self.mutation):
                        if self.mutation[0] in ['ASP', 'GLU']:
                            charge += 1
                        if self.mutation[0] in ['ARG', 'LYS', 'HIP']:
                            charge -= 1

                        try: # Add counter ions to the reference water sphere to match the total of the protein sphere
                            xyz = counter_ions.minimize_coulomb_on_sphere(abs(charge), int(float(self.radius))-11, self.sphere)
                            self.ions = len(xyz)
                            offset = next(id)
                            for resi, coords in enumerate(xyz, start=6):
                                x, y, z = coords[0], coords[1], coords[2]
                                ion = 'CLA' if charge <= 0 else 'SOD'
                                pdb.write(f"ATOM   {offset:4}  {ion} {ion}   {resi:3}     {x:7.3f} {y:7.3f} {z:7.3f}  0.00  0.00           \n")
                                offset = next(id)
                            self.systemsize = offset
                        except: pass
                    else: pass
            
        elif self.system == 'vacuum':
            pass
            
    # Removing water molecules clashing with solute atoms or ions
    def select_waters(self):
        cofactor_coordinates = []
        waters_remove = []
        waters = {}
        
        # Add cofactor coordinates
        if self.cofactors:
            for cofactor in self.cofactors:
                with open(f"{cofactor}.pdb") as pdb:
                    for atom in pdb:
                        if atom.startswith(('ATOM','HETATM')):
                            atom = IO.pdb_parse_in(atom)
                            xyz  = [atom[8], atom[9], atom[10]]
                            cofactor_coordinates.append(xyz)
        
        # Store water molecules in dict
        with open ('water.pdb') as water_in, \
            open (f"{self.directory}/inputfiles/water.pdb", 'w') as water_out:
                water_out.write('{} SPHERE\n'.format(self.radius))
                for atom in water_in:
                    atom = IO.pdb_parse_in(atom)
                    resi = atom[6]
                    xyz  = [atom[8], atom[9], atom[10]]
                    waters.setdefault(resi, []).append(atom)

                    # Determine which water molecules overlap with cofactor atoms
                    for cofactor in cofactor_coordinates:
                        if f.euclidian_overlap(xyz, cofactor, 1.6) == True:
                            waters_remove.append(resi)
                
                # Write out only non-clashing water molecules to water.pdb
                for water in waters:
                    if water not in waters_remove:
                        for water in waters[water]:
                            water_out.write(f"{IO.pdb_parse_out(water)}\n")

    # Write the qprep final qprep input file                
    def write_qprep(self):
        self.qprep = f"{self.directory}/inputfiles/qprep.inp"
        replacements = {'PRM':self.prm_merged,
                        'PDB': 'complex.pdb',
                        'CENTER':'{} {} {}'.format(*self.sphere),
                        'SPHERE':self.radius}
        
        if self.system == 'protein':
            replacements['SOLVENT'] = '4 water.pdb'    
        elif self.system == 'water':
            replacements['SOLVENT'] = '1 HOH'
        elif self.system == 'vacuum':
            replacements['solvate']='!solvate'
        
        libraries = [f"{self.forcefield}.lib"]
        if self.cofactors:
            for cofactor in self.cofactors:
                libraries.append(f"{cofactor}.lib")
                
        # Modify the template qprep input file
        with open(f"{s.INPUT_DIR}/qprep_QresFEP.inp") as qprep_tmplt:
            with open(self.qprep, 'w') as qprep_out:
                if self.dual == True:
                    qprep_out.write(f'rl {self.MUTresn}.lib\n')
                
                for line in qprep_tmplt:
                    line = IO.replace(line, replacements)
                    
                    # Add additional library files
                    if line.split()[0] == '!Added':
                        for lib in libraries:
                            qprep_out.write(f'rl {lib}\n')
                        continue
                    
                    # Add sulfide bridges between cysteine residues (SG atoms)
                    if line.split()[0] == '!addbond':
                        for cyx in self.CYX:
                            outline = 'addbond {}:SG {}:SG y\n'.format(*cyx)
                            qprep_out.write(outline)
                        continue
                    
                    qprep_out.write(line)
        
    # Run qprep
    def run_qprep(self):
        os.chdir(self.directory + '/inputfiles/')

        q_commands = getattr(s, self.preplocation)
        qprep = q_commands['QPREP']
        os.system(f'{qprep} < qprep.inp > qprep.out')

        os.chdir('../../')

    # Retrieve the values for the lambda windows
    def get_lambdas(self):
        if self.sampling == 'exponential':
            self.lambdas_1 = IO.get_lambdas(self.windows, 'exponential')
            self.lambdas_2 = IO.get_lambdas(self.windows, 'reverse_exponential')
        else:
            self.lambdas = IO.get_lambdas(self.windows, self.sampling)


    # Write the equilibration (eq) MD input files
    def write_EQ(self):
        # Set anchor restrain on mutable residue CA (atom 19) to the spphere system
        if self.system in ('water', 'vacuum'):
            self.replacements['WATER_RESTRAINT'] = '19 19 1.0 0 0'
        # Set final solute atom to define sequence restraint limit                
        self.replacements['ATOM_END'] = f"{self.systemsize}"
        
        hyb = self.PDB[int(self.PDB2Q[self.chain][self.mutation[1]])]
        topology = self.directory + '/inputfiles/top_p.pdb'
        
        with open(topology, 'r') as top:
            _ = itertools.count()
            for atom in top:
                if atom.startswith('ATOM'):
                    atom = atom.split() # Store variable to improve readability
                    id, atn, resn, resi, x, y, z = atom[1], atom[2], atom[3], atom[4], atom[5], atom[6], atom[7]

                    if self.dual:
                        if '2' in resn: # Hybrid residue name is {WT}2{MUT}
                            i = next(_) # Store the hybrid residue atom informations
                            hyb[i][1], hyb[i][2], hyb[i][8], hyb[i][9], hyb[i][10] = int(id), atn, x, y, z
                            
                            if atn in self.atoms: # Store the atom ids in the global self.atoms dict
                                self.atoms[atn] = id
                    else: # Single topology atom information - for tripeptide or protein
                        match = '2' if self.tpt == 'X' and self.system == 'water' else '3'
                        resi_match = (resi == match) if self.system == 'water' else (resi == self.PDB2Q[self.chain][self.mutation[1]])
                        if resi_match and atn in self.atoms:
                            self.atoms[atn] = id

        # Determine distance restraints between heavy atoms of the wild-type and mutant side chains
        self.dist_rest = []
        wt, mut = IO.restraint_matrix(self.mutation) # Returns heavy atom names for wt and mut residue
        ha_pairs = IO.heavy_atom_match(wt, mut, hyb) # Determines equivalent wt and mut heavy atom pairs to restrain
        for ha_pair in ha_pairs:
            self.dist_rest.append(f"{self.atoms[ha_pair[0]]} {self.atoms[ha_pair[1]]} 0.0 0.5 10.0 0")
        self.replacements['DIST'] = '\n'.join(self.dist_rest)

        # Include wall restraints for ions in reference tripeptide spheres IF mutation involves a charged residue
        try:
            if any(res in ['ASP', 'GLU', 'LYS', 'ARG', 'HIP'] for res in self.mutation) and self.system == "water":
                first = self.systemsize + 1
                last = self.systemsize + self.ions
                self.replacements['WALL'] = f"{first} {last} {int(float(self.radius)) - 5} 1.0 0 0 0"
            else:
                self.replacements['WALL'] = "\n"
        except:
            self.replacements['WALL'] = "\n"

        # Modify the template equilibration input files
        for eq_tmplt in glob.glob(s.INPUT_DIR + '/eq*.inp'):
            eq_file = f"{self.directory}/inputfiles/{eq_tmplt.split('/')[-1]}"
            with open(eq_tmplt) as eq_in, open(eq_file, 'w') as eq_out:
                    for line in eq_in:
                        eq_out.write(IO.replace(line, self.replacements))

    # Write the production (md) MD input files
    def write_MD(self):
        # Set anchor restrain on mutable residue CA (atom 19) to the spphere system
        if self.system in ('water', 'vacuum'):
            self.replacements['WATER_RESTRAINT'] = '19 19 1.0 0 0'
        
        # Determine the lamba values for each md file
        if self.sampling == 'exponential':
            for stage, lambdas in enumerate([self.lambdas_1, self.lambdas_2], start=1):
                for window in range(self.windows):
                    if self.start == '0.5' and stage == 1:
                        value = 1.000 - float(list(reversed(lambdas))[window])
                    else:
                        value = float(lambdas[window])
                    lambda1 = f'{value:.3f}'.replace(".", "")
                    lambda2 = f'{1.000 - value:.3f}'.replace(".", "")

                    # Handle lambda float replacements
                    if self.FromGly or (self.start == '0.5' and stage == 1):
                        self.replacements['FLOAT_LAMBDA1'] = f'{float(lambda2) / 1000:.3f}'
                        self.replacements['FLOAT_LAMBDA2'] = f'{float(lambda1) / 1000:.3f}'
                    else:
                        self.replacements['FLOAT_LAMBDA1'] = f'{float(lambda1) / 1000:.3f}'
                        self.replacements['FLOAT_LAMBDA2'] = f'{float(lambda2) / 1000:.3f}'
                    self.replacements['FILE'] = f"md_{stage}_{lambda1}_{lambda2}"
                    
                    # Handle file section
                    if lambdas[window] == '1.000':
                        md_tmplt = s.INPUT_DIR + '/md_1000_0000.inp'
                        self.replacements['md_1000_0000.dcd'] = f'md_{stage}_1000_0000.dcd'
                        self.replacements['md_1000_0000.en'] = f'md_{stage}_1000_0000.en'
                        self.replacements['md_1000_0000.re'] = f'md_{stage}_1000_0000.re'
                    else:
                        md_tmplt = s.INPUT_DIR + '/md_XXXX_XXXX.inp'
                    md_file = f"{self.directory}/inputfiles/md_{stage}_{lambda1}_{lambda2}.inp"

                    # Write the md file
                    with open(md_tmplt) as md_in, open(md_file, 'w') as md_out:
                        for line in md_in:
                            md_out.write(IO.replace(line, self.replacements))
                    
                    # Store previous file
                    self.replacements['FILE_N'] = f'md_{stage}_{lambda1}_{lambda2}'
        else:
            for window in range(self.windows):
                lambda1 = self.lambdas[window].replace(".", "")
                lambda2 = self.lambdas[::-1][window].replace(".", "")

                self.replacements['FLOAT_LAMBDA1'] = f'{float(lambda1) / 1000:.3f}'
                self.replacements['FLOAT_LAMBDA2'] = f'{float(lambda2) / 1000:.3f}'

                self.replacements['FILE'] = f"md_{lambda1}_{lambda2}"

                md_tmplt = (s.INPUT_DIR + '/md_1000_0000.inp' if window == 0 else s.INPUT_DIR + '/md_XXXX_XXXX.inp')
                md_file = f"{self.directory}/inputfiles/md_{lambda1}_{lambda2}.inp"

                with open(md_tmplt) as md_in, open(md_file, 'w') as md_out:
                    for line in md_in:
                        md_out.write(IO.replace(line, self.replacements))
                
                # Store previous file
                self.replacements['FILE_N'] = f'md_{lambda1}_{lambda2}'

    # Write the bash run script for the HPC cluster (run{HPC}.sh)
    def write_runfile(self):
        # Replace for template definitions 
        replacements = IO.merge_two_dicts(self.replacements, getattr(s, self.cluster))
        replacements['FEPS'] = ' '.join(self.FEPlist)
        replacements['JOBNAME'] = 'QresFEP'
        replacements['RESTART'] = 'md_1_0000_1000.re' if self.start == '1' else 'eq5.re' if self.start == '0.5' else replacements.get('RESTART')

        # Define the template and target run file
        run_tmplt = f"{s.INPUT_DIR}/{'run_exponential.sh' if self.sampling == 'exponential' else 'run_benchmark.sh'}"
        run_file = f"{self.directory}/inputfiles/run{self.cluster}.sh"
        
        eq_files = sorted(glob.glob(self.directory + '/inputfiles/eq*.inp'))

        # Write the runHPC.sh file
        with open(run_tmplt) as run_in, open(run_file, 'w') as run_out:
                for line in run_in: # Remove SLURM account name if not required
                    if line.strip() == '#SBATCH -A ACCOUNT' and 'ACCOUNT' not in replacements:
                        line = ''
                    run_out.write(IO.replace(line, replacements))

                    # Add qdyn lines for equilibration phase
                    if line.strip() == '#EQ_FILES':
                        for eq_file in eq_files:
                            file_base = eq_file.split('/')[-1][:-4]
                            outline = f"time srun $qdyn {file_base}.inp > {file_base}.log\nwait\n"
                            run_out.write(outline)

                    # Add qdyn lines for production phase
                    for tag, pattern in [('#RUN_FILES', 'md_*.inp'), 
                                         ('#RUN_1_FILES', 'md_1*.inp'), 
                                         ('#RUN_2_FILES', 'md_2*.inp')]:
                        if line.strip() == tag:
                            md_files = reversed(sorted(glob.glob(f"{self.directory}/inputfiles/{pattern}")))
                            for md_file in md_files:
                                file_base = md_file.split('/')[-1][:-4]
                                outline = f"echo {file_base}\ntime srun $qdyn {file_base}.inp > {file_base}.log\nwait\n"
                                run_out.write(outline)

    # Define the time step and runtime for final eq step and md steps
    def settimestep(self):
        if self.timestep == '1fs':
            self.replacements['NSTEPS1'] = '500000'
            self.replacements['NSTEPS2'] = '10000'
            self.replacements['STEPSIZE'] = '1.0'
            self.replacements['STEPTOGGLE'] = 'off'
            
        if self.timestep == '2fs':
            self.replacements['NSTEPS1'] = '1250000'
            self.replacements['NSTEPS2'] = '10000'
            self.replacements['STEPSIZE'] = '2.0'
            self.replacements['STEPTOGGLE'] = 'on' 
                        
    # Write the submit file (FEP_submit.sh) for job submission on the HPC cluster
    def write_submitfile(self):
        IO.write_submitfile_benchmark(self.directory, self.replacements)
    
    # Write the FEP (FEP.fep) files for requested topology strategy
    def write_FEPfile(self):
        run.write_dualFEPfile() if self.dual else run.write_singleFEPfile()
        
    # Write the single topology FEP files
    def write_singleFEPfile(self):
        for fep_tmplt in sorted(glob.glob(f"{self.FEPdir}/FEP*.fep")): # Takes the mutation template FEP files
            fep_file = f"{self.directory}/inputfiles/{fep_tmplt.split('/')[-1]}"
            with open (fep_tmplt) as fep_in, open (fep_file, 'w') as fep_out:
                headers = ['[change_charges]', '[FEP]', '[atom_types]', '[change_atoms]', '[softcore]', '[bond_types]']
                fep_atoms = {} # Maps FEP atom id to topology id
                block = 0
                id = itertools.count(1)
                for line in fep_in:
                    if not line.strip(): # Simply write out empty lines
                        fep_out.write('\n')
                        continue
                    if line.strip() == '[atoms]': # Gatekeeper for [atoms] section
                        fep_out.write(line)
                        block = 1
                        continue
                        
                    if line.strip() == '[change_bonds]': # Gatekeeper for [change_bonds] section
                        fep_out.write(line)
                        block = 2
                        continue
                        
                    if line.strip() in headers: # Gatekeeper for other sections from headers list
                        block = 0

                    if block == 0: # Simply write out lines of sections from headers list
                        fep_out.write(line)
                    
                    elif block == 1: # Write out the topology ids for the [atoms] section
                        line = line.split()
                        atom = str(self.atoms[line[3]]) # Finds the topology atom id
                        fep_atoms.setdefault(str(next(id)), atom) # Add the atom ids to the fep_atoms map
                        comment = line[4] if len(line) > 4 else ' ' # Adds the atom name as comment
                        fep_out.write(f"{line[1]:7s}{atom:<7s}{'!':5s}{line[3]:5s}{comment:5s}\n")
                            
                    elif block == 2: # Writes out the topology ids for the [change_bonds] section
                        line = line.split()
                        fep_out.write(f'{fep_atoms[line[0]]:<7s}{fep_atoms[line[1]]:<7s}{line[2]:5s}{line[3]:5s}\n')

                block = None

                            
    # Write the dual topology FEP files
    def write_dualFEPfile(self):
        # The vdW paramers need to be extracted from the reference .prm file
        prmfiles = [f"{s.FF_DIR}/{self.forcefield}.prm"]
        if self.nonAA == True:
            prmfiles.append(f"{self.mutation[2]}.prm")
            
        vdw_prms = IO.read_prm(prmfiles)['[atom_types]']
        
        ## Next, we start to construct the FEP file
        fepfiles = self.FEPlist.copy()
        if self.FromGly: # Make sure to always start from the larger residue for glycine mutations
            fepfiles.reverse()

        for i, fepfile in enumerate(self.FEPlist):
            with open(f"{self.directory}/inputfiles/{fepfiles[i]}", 'w') as fep_out:
                # First create the file header
                fep_out.write(f'! Dual topology QresFEP: {self.MUTresn}\n\n')
                fep_out.write('[FEP]\n')
                fep_out.write('   states 2\n')
                fep_out.write('   softcore_use_max_potential on\n\n')

                # Write out [atoms] section
                fep_out.write('[atoms]\n')
                hyb_resi = int(self.PDB2Q[self.chain][self.mutation[1]])
                for id, atom in enumerate(self.PDB[hyb_resi], start=1):
                    fep_out.write(f"{id:4d} {atom[1]:4d} ! {atom[2]:<4s}\n")

                # Write out [change_charges] section
                fep_out.write('\n[change_charges]\n')
                for id, atom in enumerate(self.merged_lib['[atoms]'], start=1):
                    atom = atom.split()
                    atn, charge = atom[1], atom[3]
                    
                    if atn.islower(): # Mutant side chain atoms
                        if self.FromGly and atn == 'ha': # Only charge mutant HA in FEP1 stage of Gly mutations
                            charges = ('0.000', '0.060') if fepfile == 'FEP1.fep' else ('0.060', '0.060')
                        else: # Charge mutant atoms in the FEP2 stage in all other cases
                            charges = ('0.000', '0.000') if fepfile == 'FEP1.fep' else ('0.000', charge)
                    else: # Wild-type side chain atoms
                        if self.ToGly and atn == 'HA': # Only decharge wild-type HA in FEP2 stage of Gly mutations
                            charges = ('0.060', '0.060') if fepfile == 'FEP1.fep' else ('0.060', '0.000')
                        elif atn in self.backbone:
                            if self.FromGly and atn == 'CA': # Adjust CA charge in FEP1 stage of mutations from Gly
                                charges = ('0.080', '0.140') if fepfile == 'FEP1.fep' else ('0.140', '0.140')
                            elif self.ToGly and atn == 'CA': # ADjust CA charge in FEP2 stage of mutations to Gly
                                charges = ('0.140', '0.140') if fepfile == 'FEP1.fep' else ('0.140', '0.080')
                            else: # Always keep backbone charges (C, O, CA, N, H)
                                charges = (charge, charge) 
                        else: # Decharge wild-type atoms in the FEP1 stage in all other cases
                            charges = (charge, '0.000') if fepfile == 'FEP1.fep' else ('0.000', '0.000')
                            
                    fep_out.write(f"{id:4d}{charges[0]:>10s}{charges[1]:>10s}\n")

                # Write out [atom_types] section
                fep_out.write('\n[atom_types]\n')
                for prm in vdw_prms:
                    prm = prm.split()
                    if len(prm) >= 7:
                        fep_out.write(f"{prm[0]:<8s}{prm[1]:>10s}{prm[3]:>10s}{'  0.0  0.0'}{prm[4]:>10s}{prm[5]:>10s}{prm[6]:>10s}\n")

                # Write out [change_atoms] section
                fep_out.write('\n[change_atoms]\n')
                for id, atom in enumerate(self.merged_lib['[atoms]'], start=1):
                    atom = atom.split()
                    atn, att = atom[1], atom[2]

                    if atn.islower(): # Mutant side chain atoms
                        if self.FromGly and atn == 'ha': # Only charge mutant HA in FEP1 stage of Gly mutations
                            types = ('DUM', 'H140') if fepfile == 'FEP1.fep' else ('H140', 'H140')
                        else: # Charge mutant atoms in the FEP2 stage in all other cases
                            types = ('DUM', att) if fepfile == 'FEP1.fep' else (att, att)
                    else: # Wild-type side chain atoms
                        if self.ToGly and atn == 'HA': # Only decharge wild-type HA in FEP2 stage of Gly mutations
                            types = ('H140', 'H140') if fepfile == 'FEP1.fep' else ('H140', 'DUM')
                        elif atn in self.backbone: # Always keep backbone charges (C, O, CA, N, H)
                                types = (att, att)
                        else: # Decharge wild-type atoms in the FEP1 stage in all other cases
                            types = (att, att) if fepfile == 'FEP1.fep' else (att, 'DUM')
                    
                    fep_out.write(f"{id:4d}{types[0]:>10s}{types[1]:>10s}\n")


                # Write out [softcore] section
                fep_out.write('\n[softcore]\n')
                for id, atom in enumerate(self.merged_lib['[atoms]'], start=1):
                    atom = atom.split()
                    atn = atom[1]
                    if atn in self.backbone: # Don't apply softcore to backbone atoms
                        sc = ('0', '0')
                    else: # Apply softcore on all side-chain atoms in the intermediate state
                        sc = ('0', '20') if fepfile == 'FEP1.fep' else ('20', '0')
                    fep_out.write(f"{id:4d}{sc[0]:>10s}{sc[1]:>10s}\n")
                
                # Write excluded pairs section
                fep_out.write('\n[excluded_pairs]\n')
                atoms = [(atom[1], atom[2]) for atom in self.PDB[hyb_resi] if atom[2] not in self.backbone]
                wt_res = [id for id, atn in atoms if not atn.islower()] # wild-type side-chain atom numbers
                mut_res = [id for id, atn in atoms if atn.islower()] # mutant side-chain atom numbers
                for atom1 in wt_res:
                    for atom2 in mut_res:
                        fep_out.write(f'{atom1:6d}{atom2:6d}{"1":>5s}{"1":>5s}\n')
                
                # Adds bond section - no unphysical bond connections need be set to zero
                fep_out.write('\n[bond_types]\n\n[change_bonds]\n\n')
                
                # Adds angles section - types needn't be defined; only set unphysical angle connections to zero
                fep_out.write('[angle_types]\n\n[change_angles]\n')
                if self.FromGly: # connecting through common backbone CA to mutant side chain atoms for wild-type Gly mutations
                    zero_angles = [('HA3', 'CA', 'cb'), ('HA2', 'CA', 'cb'), ('HA3', 'CA', 'ha'), ('HA2', 'CA', 'ha')]
                elif self.ToGly: # for mutant Gly mutations
                    zero_angles = [('ha3', 'CA', 'CB'), ('ha2', 'CA', 'CB'), ('ha3', 'CA', 'HA'), ('ha2', 'CA', 'HA')]
                else: # for all other mutations
                    zero_angles = [('CB', 'CA', 'cb')]
                
                for atoms in zero_angles: # retrieve the topology ids for the zero-angle atoms and write to FEP file
                    ids = [self.atoms[atoms[0]], self.atoms[atoms[1]], self.atoms[atoms[2]]]
                    fep_out.write(f"{ids[0]:>6s}{ids[1]:>6s}{ids[2]:>6s}{0:>5d}{0:>5d}\n")
                
                # Add torsions section - types needn't be defined; only set unphysical angle connections to zero
                fep_out.write('\n[torsion_types]\n\n[change_torsions]\n')
                wt, mut = self.mutation[0], self.mutation[2] # Define wild-type and mutant residues
                for resi in [wt, mut]: # Define relevant side-chain atoms for each amino acid (templates)
                    if resi in ['ASP', 'GLU', 'HID', 'HIE', 'ARG', 'LYS', 'PHE', 'LEU', 'MET', 'TRP', 'TYR', 'ASN', 'GLN']:
                        atoms = ['HB2', 'HB3', 'CG']
                    elif resi in ['ILE', 'VAL', 'THR']:
                        atoms = ['HB', 'CG2']
                        atoms.append('OG1' if resi == 'THR' else 'CG1')
                    elif resi in ['CYS', 'SER', 'ALA']:
                        atoms = ['HB2', 'HB3']
                        atoms.append('HB1' if resi == 'ALA' else 'SG' if resi == 'CYS' else 'OG')
                    else: # Skip Gly
                        continue

                    # Define the core atoms and the "dummy" atoms (connection set to zero) - including Gly mutations - and set the letter case if wt or mut
                    atoms = [atom.lower() if resi == mut else atom.upper() for atom in atoms]
                    core = ['CB', 'CA'] if resi == wt else ['cb', 'CA']
                    dummies = ['HA2', 'HA3'] if self.FromGly and resi == mut else ['ha2', 'ha3'] if self.ToGly and resi == wt else ['CB'] if resi == mut else ['cb']                    

                    for atom in atoms:
                        for dum in dummies: # retrieve the topology ids for the zero-torsion atoms and write to FEP file
                            ids = [self.atoms[atom], self.atoms[core[0]], self.atoms[core[1]], self.atoms[dum]]
                            fep_out.write(f"{ids[0]:>6s}{ids[1]:>6s}{ids[2]:>6s}{ids[3]:>6s}{0:>4d}{0:>4d}\n")
                    

    # Write the qfep input file
    def write_qfep(self):
        qfep_tmplt = f"{s.ROOT_DIR}/INPUTS/qfep.inp"
        qfep_file = f"{self.directory}/inputfiles/qfep.inp"
    
        # Store values for kT
        replacements = {'kT': f.kT(float(self.temperature)),
                        'WINDOWS': str(self.windows),
                        'TOTAL_L': str(self.windows)}

        with open(qfep_tmplt) as qfep_in, open(qfep_file, 'w') as qfep_out:
            for line in qfep_in:
                line = IO.replace(line, replacements) # Replace the template value
                if line == '!ENERGY_FILES\n':
                    if self.sampling != 'exponential': # Add energy files to read for qfep
                        for lambda1, lambda2 in zip(self.lambdas, reversed(self.lambdas)):
                            qfep_out.write(f"md_{lambda1.replace('.', '')}_{lambda2.replace('.', '')}.en\n")
                    continue # Don't add energy files for exponential sampling - these are different per FEP stage; runHPC.sh takes care of this
                qfep_out.write(line)
                            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='QresFEP',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '===== Generate input files for running residue FEP with Q; First run protprep.py for required input files =====')
    
    required = parser.add_argument_group('required arguments')
    required.add_argument('-m', '--mutation',
                        dest = "mutation",
                        required = True,
                        help = "The intended mutation, input as $WT$RESN$MUT (e.g. ALA123GLY); residue number in original PDB file")

    required.add_argument('-mc', '--mutchain',
                        dest = "mutchain",
                        required = True,
                        help = "Specficy the PDB chain of the intended mutations")

    required.add_argument('-S', '--system',
                        dest = "system",
                        required = True,
                        choices = ['protein', 'water', 'vacuum'],
                        help = "Specify system type; either protein, water or vacuum")

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('-t', '--tripeptide',
                        dest = "tripeptide",
                        required = False,
                        default = 'A',
                        choices = ['A', 'G', 'X', 'Z'],
                        help = "If system is water, use to specify flanking residues of mutable X in reference tripeptide: \
                                A: AXA (Ala)  |  G: GXG  (Gly) |  X: X (None)  |  Z: ZXZ (Natural sequence)")

    optional.add_argument('-d', '--dual',
                        dest = "dual",
                        default = False,
                        action = 'store_true',
                        help = "Turn on for dual topology FEP")

    optional.add_argument('-c', '--cofactors',
                        nargs='*',
                        dest = "cofactors",
                        required = False,
                        help = "If cofactor(s), input name(s) (e.g. ligand); .pdb, .prm and .lib files are required")

    optional.add_argument('-f', '--forcefield',
                        dest = "forcefield",
                        default = 'OPLSAAM',
                        choices = ['OPLSAAM', 'OPLS2015', 'OPLS2005', 'SIDECHAIN', 'AMBER14sb','CHARMM36'],
                        help = "Use to specficy forcefield")
        
    optional.add_argument('-w', '--windows',
                        dest = "windows",
                        default = '50',
                        help = "Use to specify number of lambda windows per FEP stage")

    optional.add_argument('-s', '--sampling',
                        dest = "sampling",
                        default = 'linear',
                        choices = ['linear', 'sigmoidal', 'exponential', 'reverse_exponential'],
                        help = "Use to specify sampling strategy")

    optional.add_argument('-l', '--start',
                        dest = "start",
                        default = '1',
                        choices = ['1', '0.5', '0'],
                        help = "Specify starting point of FEP simulations; either start (1), middle (0.5) or end (0) lambda point")

    optional.add_argument('-ts', '--timestep',
                        dest = "timestep",
                        default = "2fs",
                        choices = ['1fs','2fs'],
                        help = "Use to specify timestep")
    
    optional.add_argument('-T', '--temperature',
                        dest = "temperature",
                        default = '298',
                        help = "Use to specify intended simulation temperature")
    
    optional.add_argument('-r', '--replicates',
                        dest = "replicates",
                        default = '10',
                        help = "Use to specify intended number of simulation replicates")   
    
    optional.add_argument('-C', '--cluster',
                        dest = "cluster",
                        default = s.DEFAULT,
                        help = "Use to specify HPC cluster to run simulations")
    
    optional.add_argument('-P', '--preplocation',
                        dest = "preplocation",
                        default = s.DEFAULT,
                        help = "Use to specify location of Q executables QresFEP ought to use")

    args = parser.parse_args()
    run = Run(mutation = args.mutation,
              mutchain = args.mutchain,
              system = args.system,
              tripeptide = args.tripeptide,
              dual = args.dual,
              cofactors = args.cofactors,
              forcefield = args.forcefield,
              windows = args.windows,
              sampling = args.sampling,
              start = args.start,
              timestep = args.timestep,
              temperature = args.temperature,
              replicates = args.replicates,
              cluster = args.cluster,
              preplocation = args.preplocation,)
    
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
