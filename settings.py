import os
# Root directory of the setup FEP modules
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

# The directories to the input FF and run related input files are given here
FF_DIR = os.path.join(ROOT_DIR, "FF")
INPUT_DIR = os.path.join(ROOT_DIR, "INPUTS")
# Dicionary of locations of Q executables
Q_DIR = {'CSB':     '/home/koenekoop/software/q6/bin/',
         'TETRA':   '/proj/alexandria/users/x_lucko/software/q6/bin/',
         'DARDEL':  '/cfs/klemming/home/l/lucko/q6/bin/',
         'LUMI':    '/scratch/project_465000920/lucien/software/q6/bin/',
         'LOCAL':   '/home/nvandebrug/software/q6/bin/'
        }
BIN = os.path.join(ROOT_DIR, "bin")
SCHROD_DIR = '/home/apps/apps/schrodinger2024-4/'

DEFAULT = 'CSB'

# CLUSTER INPUTS. To add your own cluster, use the same input as below
CSB = {'NODES'        : '1',
       'NTASKS'       : '8',
       'TIME'         : '1-00:00:00',  # d-hh:mm:ss
       'PARTITION'    : 'CLUSTER-AMD',
       'EXCLUDE'      : 'compute-0-33,compute-0-15',
       'MODULES'      : '\n',
       'QDYN'         : 'qdyn=' + Q_DIR['CSB'] + 'qdynp',
       'QPREP'        : Q_DIR['CSB'] + 'qprep',
       'QFEP'         : Q_DIR['CSB'] + 'qfep',
       'QCALC'        : Q_DIR['CSB'] + 'qcalc'
      }

TETRA  = {'NODES'      : '1',
          'NTASKS'     : '8',
          'TIME'       : '1-00:00:00',  # d-hh:mm:ss
          'PARTITION'  : 'tetralith',
          'MODULES'    : '\n', # Add a \n for every added module
          'QDYN'       : 'qdyn=' + Q_DIR['TETRA'] + 'qdynp',
          'QPREP'      : Q_DIR['TETRA'] + 'qprep',
          'QFEP'       : Q_DIR['TETRA'] + 'qfep',
          'ACCOUNT'    : 'naiss2025-3-3'
        }

DARDEL = {'NODES'      : '1',
          'NTASKS'     : '8',
          'TIME'       : '1-00:00:00',  # d-hh:mm:ss
          'PARTITION'  : 'shared',
          'MODULES'    : '\n', # Add a \n for every added module
          'QDYN'       : 'qdyn=' + Q_DIR['DARDEL'] + 'qdynp',
          'QPREP'      : Q_DIR['DARDEL'] + 'qprep',
          'QFEP'       : Q_DIR['DARDEL'] + 'qfep',
          'ACCOUNT'    : 'naiss2024-3-13'
          }

LUMI   = {'NODES'      : '1',
          'NTASKS'     : '8',
          'TIME'       : '0-05:00:00',  # d-hh:mm:ss
          'PARTITION'  : 'small',
          'MODULES'    : '\n', # Add a \n for every added module
          'QDYN'       : 'qdyn=' + Q_DIR['LUMI'] + 'qdynp',
          'QPREP'      : Q_DIR['LUMI'] + 'qprep',
          'QFEP'       : Q_DIR['LUMI'] + 'qfep',
          'ACCOUNT'    : 'project_465000920' 
        }

KEBNE = {'NODES'      : '1',
         'NTASKS'     : '8',
         'TIME'       : '0-24:00:00',  # d-hh:mm:ss
         'MODULES'    : '\n', # Add a \n for every added module
         'QDYN'       : 'qdyn=qdyn5p', #fix qdyn= !!!!!
         'QPREP'      : 'qprep5', # NOTE: change to where you are setting up, not where you are running!
         'QFEP'       : 'qfep5',
         'QCALC'      : 'qcalc5',
         'ACCOUNT'    : 'naiss2023-3-5'
        }

RACKHAM = {'NODES'      : '1',
           'NTASKS'     : '10',
           'TIME'       : '0-10:00:00',  # d-hh:mm:ss
           'MODULES'    : '\n', # Add a \n for every added module
           'QDYN'       : 'qdyn=qdyn5p', #fix qdyn= !!!!!
           'QPREP'      : 'qprep5', # NOTE: change to where you are setting up, not where you are running!
           'QFEP'       : 'qfep5',
           'QCALC'      : 'qcalc5',
           'ACCOUNT'    : 'SNIC2021-3-1'
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
