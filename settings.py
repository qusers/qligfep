import os
# Root directory of the setup FEP modules
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

# The directories to the input FF and run related input files are given here
FF_DIR = os.path.join(ROOT_DIR, "FF")
INPUT_DIR = os.path.join(ROOT_DIR, "INPUTS")
# Dicionary of locations of Q executables
Q_DIR = {'CSB':'/home/apps/apps/Q/5.10/bin/',
         'LOCAL':'/home/apps/apps/Q/5.10/bin/'
        }
BIN = os.path.join(ROOT_DIR, "bin")
SCHROD_DIR = '/home/apps/apps/schrodinger2020-3/'

# CLUSTER INPUTS. To add your own cluster, use the same input as below
CSB = {'NODES'        : '1',
       'NTASKS'       : '8',
       'TIME'         : '0-08:00:00',  # d-hh:mm:ss
       'PARTITION'    : 'CLUSTER-AMD',
       'MODULES'      : '\n',
       'QDYN'         : 'qdyn=' + Q_DIR['CSB'] + 'qdyn5p',
       'QPREP'        : Q_DIR['CSB'] + 'qprep5',
       'QFEP'         : Q_DIR['CSB'] + 'qfep5',
       'QCALC'        : Q_DIR['CSB'] + 'qcalc5'
      }

TETRA  = {'NODES'      : '1',
          'NTASKS'     : '8',
          'TIME'       : '2-00:00:00',  # d-hh:mm:ss
          'PARTITION'  : 'tetralith',
          'MODULES'    : '\n', # Add a \n for every added module
          'QDYN'       : 'qdyn=qdyn5p', #fix qdyn= !!!!!
          'QPREP'      : 'qprep', # NOTE: change to where you are setting up, not where you are running!
          'QFEP'       : 'qfep5',
          'ACCOUNT'    : 'naiss2023-3-5'
        }

DARDEL = {'NODES'      : '1',
          'NTASKS'     : '8',
          'TIME'       : '2-00:00:00',  # d-hh:mm:ss
          'PARTITION'  : 'shared',
          'MODULES'    : '\n', # Add a \n for every added module
          'QDYN'       : 'qdyn=qdyn5p', #fix qdyn= !!!!!
          'QPREP'      : 'qprep', # NOTE: change to where you are setting up, not where you are running!
          'QFEP'       : 'qfep5',
          'ACCOUNT'    : 'naiss2023-3-5'
        }

KEBNE = {'NODES'      : '1',
         'NTASKS'     : '8',
         'TIME'       : '0-24:00:00',  # d-hh:mm:ss
         'MODULES'    : '\n', # Add a \n for every added module
         'QDYN'       : 'qdyn=qdyn5p', #fix qdyn= !!!!!
         'QPREP'      : 'qprep5', # NOTE: change to where you are setting up, not where you are running!
         'QFEP'       : 'qfep5',
         'QCALC'      : 'qcalc5',
         'ACCOUNT'    : 'SNIC2021-3-1'
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
