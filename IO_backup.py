import re
import shlex
from subprocess import check_output
import os
import stat
import numpy as np
import pandas as pd

import functions as f
import settings as s

## Some useful objects TO DO add GLH etc.
charged_res = {'HIS': {'HD1' : 'HID',
                       'HE2' : 'HIE'},
               'GLU': {'HE2' : 'GLH'},
               'ASP': {'HD2' : 'ASH'}
              }    

def pdb_parse_in(line, include=('ATOM','HETATM')):
    """
    Takes a pdb file line and parses it into a list, according to Atomic Coordinate Entry Format 
    v3.3
    """
    at_entry = []
    line = line.strip('\n')
    if line.startswith(include):
        at_entry.append(line[0:6])              #  0 ATOM/HETATM
        at_entry.append(int(line[6:11]))        #  1 ATOM serial number
        at_entry.append(line[12:16].strip())    #  2 ATOM name
        at_entry.append(line[16:17])            #  3 Alternate location indicator
        at_entry.append(line[17:21].strip())    #  4 Residue name
        at_entry.append(line[21:22])            #  5 Chain identifier
        at_entry.append(int(line[22:26]))       #  6 Residue sequence number
        at_entry.append(line[26:27])            #  7 Code for insertion of residue
        at_entry.append(float(line[30:38]))     #  8 Orthogonal coordinates for X
        at_entry.append(float(line[38:46]))     #  9 Orthogonal coordinates for Y
        at_entry.append(float(line[46:54]))     # 10 Orthogonal coordinates for Z
        # These entries can be empty
        try:
            at_entry.append(float(line[54:60])) # 11 Occupancy
            
        except:
            at_entry.append(0.0)                # 11 Empty Occupancy
            
        try:
            at_entry.append(float(line[60:66])) # 12 Temperature factor
            
        except:
            at_entry.append(0.0)                # 12 Empty Temperature factor
            
        try:
            at_entry.append(line[76:78])        # 13 Element symbol
            
        except:
            at_entry.append('  ')               # 13 Empty Element symbol
            
        try:
            at_entry.append(line[78:80])        # 14 Charge on atom
            
        except:
            at_entry.append('  ')               # 14 Empty charge
        
    else:
        at_entry = line
    
    return at_entry
    
def pdb_parse_out(line):
    """
    Takes a list and parses it into a pdb writeable line
    """
    if len(line[2]) <= 3: 
        line = '{:6s}{:5d}  {:3s}{:1s}{:4s}{:1s}{:4d}{:1s}   '\
               '{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format(*line)
            
    elif len(line[2]) == 4: 
        line = '{:6s}{:5d} {:4s}{:1s}{:4s}{:1s}{:4d}{:1s}   '\
               '{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}'.format(*line)
    return line

def replace(string, replacements):
    pattern = re.compile(r'\b(' + '|'.join(replacements.keys()) + r')\b')
    replaced_string = pattern.sub(lambda x: replacements[x.group()], string)
    return replaced_string

def run_command(executable, options, string = False):
    """
    Takes three variables, the executable location and its options as strings and a tag if the
    options need to be split or not (e.g. Q runs with one string), and runs the program.
    Returns the output of that program as an unformatted string.
    """
    if string == False:
        args = shlex.split(executable + options)
        out = check_output(args)
        print(' '.join(args))
    else:
        os.system(executable + options)
        out = None

    return out

def AA(AA):
    """
    Handy dictionary to convert 3 letter AA code to one and vice versa
    """
    threeAA = {'CYS': 'C', 'CYX': 'C', 'ASH': 'D', 'ASP': 'D', 'SER': 'S', 
               'GLN': 'Q', 'LYN': 'K', 'LYS': 'K', 'ILE': 'I', 'PRO': 'P', 
               'THR': 'T', 'PHE': 'F', 'ASN': 'N', 'GLY': 'G', 'HID': 'H', 
               'HIP': 'H', 'HIE': 'H', 'HIS': 'H', 'LEU': 'L', 'ARN': 'R', 
               'ARG': 'R', 'TRP': 'W', 'ALA': 'A', 'VAL': 'V', 'GLH': 'E', 
               'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
              }
    
    fourAA = { 'CCYS': 'C', 'CASP': 'D', 'CASH': 'H', 'CSER': 'S', 
               'CGLN': 'Q', 'CLYN': 'K', 'CLYS': 'K', 'CILE': 'I', 
               'CPRO': 'P', 'CTHR': 'T', 'CPHE': 'F', 'CASN': 'N', 
               'CGLY': 'G', 'CHIE': 'H', 'CHID': 'H', 'CHIP': 'H', 
               'CLEU': 'L', 'CARG': 'R', 'CARN': 'R', 'CTRP': 'W', 
               'CALA': 'A', 'CVAL': 'V', 'CGLU': 'E', 'CGLH': 'E',
               'CTYR': 'Y', 'CMET': 'M'
             }
    
    oneAA = {  'C' : 'CYS', 'D' : 'ASP', 'S' : 'SER', 'Q' : 'GLN',
               'K' : 'LYS', 'I' : 'ILE', 'P' : 'PRO', 'T' : 'THR', 
               'F' : 'PHE', 'N' : 'ASN', 'G' : 'GLY', 'H' : 'HID',
               'L' : 'LEU', 'R' : 'ARG', 'W' : 'TRP', 'A' : 'ALA',
               'V' : 'VAL', 'E' : 'GLU', 'Y' : 'TYR', 'M' : 'MET'
            }
    
    if len(AA) == 4:
        AA = fourAA[AA]
        return AA

    if len(AA) == 3:
        AA = threeAA[AA]
        return AA
        
    if len(AA) == 1:
        AA = oneAA[AA]
        return AA
    
def restraint_matrix(mutation):
    """
    Matrix stroring appropriate restraints for side chain overlap of given mutation
    """
    wt, mut = AA(AA(mutation[0])), AA(AA(mutation[2]))
    
    sidechains = {
        'ASP': ['CB', 'CG', 'OD1', 'OD2'],
        'GLU': ['CB', 'CG', 'CD', 'OE1', 'OE2'],
        'HID': ['CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'],
        'ARG': ['CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2'],
        'LYS': ['CB', 'CG', 'CD', 'CE', 'NZ'],
        'ALA': ['CB'],
        'CYS': ['CB', 'SG'],
        'PHE': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
        'ILE': ['CB', 'CG1', 'CG2', 'CD1'],
        'LEU': ['CB', 'CG', 'CD1', 'CD2'],
        'MET': ['CB', 'CG', 'SD', 'CE'],
        'VAL': ['CB', 'CG1', 'CG2'],
        'TRP': ['CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
        'TYR': ['CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH'],
        'ASN': ['CB', 'CG', 'OD1', 'ND2'],
        'GLN': ['CB', 'CG', 'CD', 'OE1', 'NE2'],
        'SER': ['CB', 'OG'],
        'THR': ['CB', 'CG2', 'OG1'],
        'GLY': []
    }
    
    residues = ['ASP', 'GLU', 'HID', 'ARG', 'LYS', 'ALA', 'CYS', 'PHE', 'ILE', 'LEU', 'MET', 'VAL', 'TRP', 'TYR', 'ASN', 'GLN', 'SER', 'THR', 'GLY']
    matrix = [[0, 'CG', 'CG', 'CG', 'CG', 'CB', 'CB', 'CG', 'CG', 'CG', 'CG', 'CG', 'CG', 'CG', ['CG','OD1'], 'CG', 'CB', 'CG', 0],
             ['CG', 0, 'CG', 'CD', 'CD', 'CB', 'CB', 'CG', 'CD', 'CD', 'CG',	'CG', 'CG', 'CG', 'CG',	['CD','OE1'], 'CB', 'CG', 0],
             ['CG', 'CG', 0, 'CG', 'CG', 'CB', 'CB', 'CG', 'CG', 'CG', 'CG', 'CG', 'CG', 'CG', 'CG', 'CG', 'CB', 'CG', 0],
             ['CG', 'CD', 'CG', 0, 'CD', 'CB', 'CB',	'CG', 'CD', 'CD', 'CG', 'CG', 'CG', 'CG', 'CG', 'CD', 'CB', 'CG', 0],
             ['CG', 'CD', 'CG', 'CD', 0, 'CB', 'CB', 'CG', 'CD', 'CD', 'CG', 'CG', 'CG', 'CG', 'CG', 'CD', 'CB', 'CG', 0],
             ['CB', 'CB', 'CB', 'CB', 'CB', 0, 'CB', 'CB', 'CB', 'CB', 'CB', 'CB', 'CB', 'CB', 'CB', 'CB', 'CB', 'CB', 0],
             ['CB', 'CB', 'CB', 'CB', 'CB', 'CB', 0, 'CB', 'CB', 'CB', 'CB', 'CB', 'CB', 'CB', 'CB', 'CB', ['CB','SG'], ['CB','SG'], 0],
             ['CG', 'CG', 'CG', 'CG', 'CG', 'CB', 'CB', 0, 'CG', 'CG', 'CG', 'CG', 'CG', 'CZ', 'CG', 'CG', 'CB', 'CG', 0],
             ['CG1', 'CD1', 'CG1', 'CD1', 'CD1', 'CB', 'CB', 'CG1', 0, 'CD1', 'CG1', ['CG1', 'CG2'], 'CG1', 'CG1', 'CG1', 'CD1', 'CB', 'CG2', 0],
             ['CG', 'CD1', 'CG', 'CD1', 'CD1', 'CB', 'CB', 'CG', 'CD1', 0, 'CG', 'CG', 'CG', 'CG', 'CG', 'CD1', 'CB', 'CG', 0],
             ['CG', 'CG', 'CG', 'CG', 'CG', 'CB', 'CB', 'CG', 'CG', 'CG', 0, 'CG', 'CG', 'CG', 'CG', 'CG', 'CB', 'CG', 0],
             ['CG1', 'CG1', 'CG1', 'CG1', 'CG1', 'CB', 'CB', 'CG1', ['CG1', 'CG2'], 'CG1', 'CG1', 0, 'CG1', 'CG1', 'CG1', 'CG1', 'CB', 'CG1', 0],
             ['CG', 'CG', 'CG', 'CG', 'CG', 'CB', 'CB', 'CG', 'CG', 'CG', 'CG', 'CG', 0, 'CG', 'CG', 'CG', 'CB', 'CG', 0],
             ['CG', 'CG', 'CG', 'CG', 'CG', 'CB', 'CB', 'CZ', 'CG', 'CG', 'CG', 'CG', 'CG', 0, 'CG', 'CG', 'CB', 'CG', 0],
             [['CG', 'OD1'], 'CG', 'CG', 'CG', 'CG', 'CB', 'CB', 'CG', 'CG', 'CG', 'CG', 'CG', 'CG', 'CG', 0, 'CG', 'CB', 'CG', 0],
             ['CG', ['CD', 'OE1'], 'CG', 'CD', 'CD', 'CB', 'CB', 'CG', 'CD', 'CD', 'CG', 'CG', 'CG', 'CG', 'CG', 0, 'CB', 'CG', 0],
             ['CB', 'CB', 'CB', 'CB', 'CB', 'CB', ['CB', 'OG'], 'CB', 'CB', 'CB', 'CB', 'CB', 'CB', 'CB', 'CB', 'CB', 0, ['CB', 'OG'], 0],
             ['CG2', 'CG2', 'CG2', 'CG2', 'CG2', 'CB', ['CB', 'OG1'], 'CG2', 'CG2', 'CG2', 'CG2', 'CG2', 'CG2', 'CG2', 'CG2', 'CG2', ['CB', 'OG1'], 0, 0],
             [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]

    restraints = pd.DataFrame(matrix, index=residues, columns=residues)
    
    wt_switch  = restraints.loc[wt, mut]
    mut_switch = restraints.loc[mut, wt]
    
    if type(wt_switch) == str:
        wt_chain  = sidechains[wt][:sidechains[wt].index(wt_switch) + 1]
        mut_chain = [x.lower() for x in sidechains[mut][:sidechains[mut].index(mut_switch) + 1]]
        if wt == 'THR' and mut == 'ILE':
            mut_chain[mut_chain.index('cg1')] = 'cg2'
        elif wt == 'ILE' and mut == 'THR':
            wt_chain[wt_chain.index('CG1')] = 'CG2'
        elif wt in ['ARG', 'GLN', 'GLU', 'LEU', 'LYS'] and mut == 'ILE':
            mut_chain.remove('cg2')
        elif wt == 'ILE' and mut in ['ARG', 'GLN', 'GLU', 'LEU', 'LYS']:
            wt_chain.remove('CG2')
            
    elif type(wt_switch) == list:
        wt_index = -1
        for atom in wt_switch:
            if atom in sidechains[wt]:
                index = sidechains[wt].index(atom)
                if index > wt_index:
                    wt_index = index
        wt_chain = sidechains[wt][:wt_index + 1]
        if wt == 'THR' and mut == 'SER':
            wt_chain.remove('CG2')

        mut_index = -1
        for atom in mut_switch:
            if atom in sidechains[mut]:
                index = sidechains[mut].index(atom)
                if index > mut_index:
                    mut_index = index
        mut_chain = [x.lower() for x in sidechains[mut][:mut_index + 1]]
        if wt == 'SER' and mut == 'THR':
            wt_chain.remove('cg2')
        
    else:
        wt_chain  = []
        mut_chain = []
    
    return wt_chain, mut_chain


def read_prm(prmfiles):
    """
    Takes a list of Q .prm files and merges them, first file is the referene .prm file
    Returns a dicitonary of the merged .prm files
    """    
    block = 0
    prm = {'[options]':[],
            '[atom_types]':[],
            '[bonds]':[],
            '[angles]':[],
            '[torsions]':[],
            '[impropers]':[]}
    
    for filename in prmfiles:
        with open(filename) as infile:
            for line in infile:
                if line == '[options]\n':
                    block = 1
                    continue                
                elif line == '[atom_types]\n':
                    block = 2
                    continue
                elif line == '[bonds]\n':
                    block = 3
                    continue
                elif line == '[angles]\n':
                    block = 4
                    continue
                elif line == '[torsions]\n':
                    block = 5
                    continue
                if line == '[impropers]\n':
                    block = 6
                    continue

                if block == 1:
                    prm['[options]'].append(line)

                if block == 2:
                    prm['[atom_types]'].append(line)

                elif block == 3:
                    prm['[bonds]'].append(line)                

                elif block == 4:
                    prm['[angles]'].append(line)                

                elif block == 5:
                    prm['[torsions]'].append(line)                

                elif block == 6:
                    prm['[impropers]'].append(line)
                
    return prm

def get_lambdas(windows, sampling):
    windows = int(windows)
    step = int(windows/2)
    lambdas = []
    lmbda_1 = []
    lmbda_2 = []
    k_dic = {'sigmoidal':-1.1, 
             'linear':1000,
             'exponential':-1.1,
             'reverse_exponential':1.1
            }
    k = k_dic[sampling]

    if sampling == 'sigmoidal': 
        for i in range(0, step + 1):
            lmbda1 = '{:.3f}'.format(0.5 * (f.sigmoid(float(i)/float(step), k) + 1))
            lmbda2 = '{:.3f}'.format(0.5 * (-f.sigmoid(float(i)/float(step), k) + 1))
            lmbda_1.append(lmbda1)
            lmbda_2.append(lmbda2)

        lmbda_2 = lmbda_2[1:]

        for i in reversed(lmbda_2):
            lambdas.append(i)

        for i in lmbda_1:
            lambdas.append(i)

    else:
        for i in range(0, windows + 1):
            lmbda = '{:.3f}'.format(f.sigmoid(float(i)/float(windows), k))
            lambdas.append(lmbda)

    lambdas = lambdas[::-1]
    return lambdas   

def write_submitfile(writedir, replacements):
    submit_in = s.ROOT_DIR + '/INPUTS/FEP_submit.sh'
    submit_out = writedir + ('/FEP_submit.sh')
    with open(submit_in) as infile, open (submit_out, 'w') as outfile:
        for line in infile:
            line = replace(line, replacements)
            outfile.write(line)

    try:
        st = os.stat(submit_out)
        os.chmod(submit_out, st.st_mode | stat.S_IEXEC)

    except:
        print("WARNING: Could not change permission for " + submit_out)

def write_submitfile_benchmark(writedir, replacements):
    submit_in = s.ROOT_DIR + '/INPUTS/FEP_submit_benchmark.sh'
    submit_out = writedir + ('/FEP_submit.sh')
    with open(submit_in) as infile, open (submit_out, 'w') as outfile:
        for line in infile:
            line = replace(line, replacements)
            outfile.write(line)

    try:
        st = os.stat(submit_out)
        os.chmod(submit_out, st.st_mode | stat.S_IEXEC)

    except:
        print("WARNING: Could not change permission for " + submit_out)

def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z

def regex_str_int(line):
    a = re.split("(\d+)", line)
    return a
    
def read_qfep(qfep):
    """
    Reads a given qfep.out file.

    returns [Zwanzig, dGfr, dGr, TI, OS, BAR]
    """
    with open(qfep) as infile:
        block = 0
        for line in infile:
            line = line.split()
            if len(line) > 3:
                if line[0] == 'ERROR:' or line[1] == 'ERROR:':
                    ERROR = True

                if line[3] == 'Free':
                    block = 1

                if line[3] == 'Termodynamic':
                    #continue
                    block = 2

                if line[3] == 'Overlap':
                    block = 3

                if line[3] == 'BAR':
                    block = 4

                if line[3] == 'Reaction':
                    block = 0

            if len(line) > 1:
                if block == 1:
                    if line[0] == '1.000000':
                        Zwanzig_r = float(line[4])

                    elif line[0] == '0.000000':
                        Zwanzig_f = float(line[2])

                        if line[5] == '-Infinity':
                            Zwanzig = np.nan

                        else:
                            Zwanzig = float(line[5])

                if block == 2 and line[0] == '0.000000':
                    try:
                        TI = line[2]
                        if line[2] == '-Infinity':
                            TI = np.nan
                    except:
                        TI = np.nan

                if block == 3 and line[0] == '0.000000':
                    if line[2] == '-Infinity':
                        OS = np.nan

                    else:
                        OS = float(line[2])

                if block == 4 and line[0] == '0.000000':
                    if line[2] == '-Infinity':
                        BAR = np.nan
                    else:
                        BAR = float(line[2])
    try: OS
    except: OS = np.nan
    try: BAR
    except: BAR = np.nan
                        
    return [Zwanzig, Zwanzig_f, Zwanzig_r, OS, BAR]

def read_qfep_esc(qfep):
    """
    Reads a given qfep.out file.

    returns [Zwanzig, dGfr, dGr, TI, OS, BAR]
    """
    with open(qfep) as infile:
        block = 0
        for line in infile:
            line = line.split()
            if len(line) > 3:
                if line[0] == 'ERROR:' or line[1] == 'ERROR:':
                    ERROR = True

                if line[3] == 'Free':
                    block = 1

                if line[3] == 'Termodynamic':
                    #continue
                    block = 2

                if line[3] == 'Overlap':
                    block = 3

                if line[3] == 'BAR':
                    block = 4

                if line[3] == 'Reaction':
                    block = 0

            if len(line) > 1:
                if block == 1:
                    if line[0] == '0.980000':
                        Zwanzig_r = float(line[4])

                    elif line[0] == '0.020000':
                        Zwanzig_f = float(line[2])

                        if line[5] == '-Infinity':
                            Zwanzig = np.nan

                        else:
                            Zwanzig = float(line[5])

                if block == 2 and line[0] == '0.020000':
                    try:
                        TI = line[2]
                        if line[2] == '-Infinity':
                            TI = np.nan
                    except:
                        TI = np.nan

                if block == 3 and line[0] == '0.020000':
                    if line[2] == '-Infinity':
                        OS = np.nan

                    else:
                        OS = float(line[2])

                if block == 4 and line[0] == '0.020000':
                    if line[2] == '-Infinity':
                        BAR = np.nan
                    else:
                        BAR = float(line[2])
    try: OS
    except: OS = np.nan
    try: BAR
    except: BAR = np.nan

    return [Zwanzig, Zwanzig_f, Zwanzig_r, OS, BAR]

def read_qfep_verbose(qfep):
    """
    Reads a given qfep.out file.

    returns [[Zwanzig, dGfr, dGr, OS, BAR]   lambda 1
                          ....               lambda ..
             [Zwanzig, dGfr, dGr, OS, BAR]]  lambda n
    """
    # Merge this and the above function?
    # Zwanzig, frwd, rv, OS, BAR
    array = [[],[],[],[],[]]
    with open(qfep) as infile:
        block = 0
        for line in infile:
            line = line.split()
            if len(line) > 3:
                if line[0] == 'ERROR:' or line[1] == 'ERROR:':
                    ERROR = True

                if line[3] == 'Free':
                    block = 1

                if line[3] == 'Termodynamic':
                    #continue
                    block = 2

                if line[3] == 'Overlap':
                    block = 3

                if line[3] == 'BAR':
                    block = 4

                if line[3] == 'Reaction':
                    block = 0
                    
            if len(line) > 1:
                if line[0] == '#':
                    continue
                
                if block == 1:
                    try:
                        array[0].append(np.float(line[5]))
                    except:
                        array[0].append(np.nan)
                    try:
                        array[1].append(np.float(line[4]))
                    except:
                        array[1].append(np.nan)
                    try:
                        array[2].append(np.float(line[2]))
                    except:
                        array[2].append(np.nan)

                if block == 3:
                    try:
                        array[3].append(np.float(line[2]))
                    except:
                        array[3].append(np.nan)

                if block == 4:
                    try:
                        array[4].append(np.float(line[2]))
                    except:
                        array[4].append(np.nan)
    
    return array   
