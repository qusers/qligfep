import re
import shlex
from subprocess import check_output
import os
import stat
import numpy as np

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
    
    if len(AA) == 3:
        AA = threeAA[AA]
        
    if len(AA) == 1:
        AA = oneAA[AA]
    return(AA)

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
    step = windows/2
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
        print "WARNING: Could not change permission for " + submit_out

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
                        Zwanzig_r = line[4]

                    elif line[0] == '0.000000':
                        Zwanzig_f = line[2]

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
    
    return [Zwanzig, Zwanzig_f, Zwanzig_r, OS, BAR]
