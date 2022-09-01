import os
# Root directory of the setup FEP modules
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

# The directories to the input FF and run related input files are given here
FF_DIR = os.path.join(ROOT_DIR, "FF")
INPUT_DIR = os.path.join(ROOT_DIR, "INPUTS")
# Dicionary of locations of Q executables
Q_DIR = {'CSB':'/Users/pantelismaragoudakis/software/Q/bin/q6/',
         'LOCAL':'/Users/pantelismaragoudakis/software/Q/bin/q6'
        }
#BIN = os.path.join(ROOT_DIR, "bin")

# some example schrodinger directory
SCHROD_DIR = '/opt/schrodinger/suites2022-1/'

# CLUSTER INPUTS. To add your own cluster, use the same input as below
ALICE = {'MAINDIR'      : '/home/jespersw/software/q6/',
         'NODES'        : '1',
         'NTASKS'       : '24',
         'TIME'         : '0-3:00:00',  # d-hh:mm:ss
         'MODULES'      : 'module load OpenMPI/3.1.3-GCC-8.2.0-2.31.1',
         'QDYN'         : 'qdyn=/home/jespersw/software/q6/bin/qdynp',
         'QPREP'        : Q_DIR['CSB'] + 'qprep',
         'QFEP'         : '/home/jespersw/software/q6/bin/qfep',
         'QCALC'        : '/home/jespersw/software/q6/bin/qcalc'
        }
