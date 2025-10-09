import os
# Root directory of the setup FEP modules
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

# The directories to the input FF and run related input files are given here
FF_DIR = os.path.join(ROOT_DIR, "FF")
INPUT_DIR = os.path.join(ROOT_DIR, "INPUTS")
# Dicionary of locations of Q executables
Q_DIR = {'CSB':'/home/jespers/software/q6/bin/',
         'LOCAL':'/home/jespers/software/q6/bin/'
        }
BIN = os.path.join(ROOT_DIR, "bin")

# some example schrodinger directory
# SCHROD_DIR = '/home/apps/apps/schrodinger2024-4/'

# CLUSTER INPUTS. To add your own cluster, use the same input as below
CSB = {'NODES'        : '1',
       'NTASKS'       : '8',
       'TIME'         : '1-00:00:00',  # d-hh:mm:ss
       'PARTITION'    : 'CLUSTER',
       'EXCLUDE'      : 'compute-0-33,compute-0-15',
       'MODULES'      : '\n',
       'QDYN'         : 'qdyn=' + Q_DIR['CSB'] + 'qdynp',
       'QPREP'        : Q_DIR['CSB'] + 'qprep',
       'QFEP'         : Q_DIR['CSB'] + 'qfep',
       'QCALC'        : Q_DIR['CSB'] + 'qcalc'
      }

LOCAL = {'NODES'      : '',
         'NTASKS'     : '',
         'TIME'       : '',  # d-hh:mm:ss
         'MODULES'    : '\n', # Add a \n for every added module
         'QDYN'       : 'qdyn=qdyn5p', #fix qdyn= !!!!!
         'QPREP'      : 'qprep5', # NOTE: change to where you are setting up, not where you are running!
         'QFEP'       : 'qfep5',
         'QCALC'      : 'qcalc5',           
         'ACCOUNT'    : ''
        }
