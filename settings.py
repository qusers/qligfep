import os
# Root directory of the setup FEP modules
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

# The directories to the input FF and run related input files are given here
FF_DIR = os.path.join(ROOT_DIR, "FF")
INPUT_DIR = os.path.join(ROOT_DIR, "INPUTS")
# Dicionary of locations of Q executables
Q_DIR = {'CSB':'/home/jespers/software/q6/bin/',
         'LOCAL':'/home/jespers/software/q6/bin/'
         #'LOCAL':'/Users/willemjespers/Software/Q6/bin/'
        }
BIN = os.path.join(ROOT_DIR, "bin")
#SCHROD_DIR = '/opt/schrodinger/suites2017-3/'
#SCHROD_DIR = '/home/apps/schrodinger2017/'
SCHROD_DIR = '/home/jespers/software/schrodinger/'

# CLUSTER INPUTS. To add your own cluster, use the same input as below
CSB = {'NODES'        : '1',
       'NTASKS'       : '16',
       'TIME'         : '0-12:00:00',  # d-hh:mm:ss
       'MODULES'      : 'module load gcc/6.2.0\n module load openmpi/2.1.0',
       'QDYN'         : 'qdyn=' + Q_DIR['CSB'] + 'qdynp',
       'QPREP'        : Q_DIR['CSB'] + 'qprep',
       'QFEP'         : Q_DIR['CSB'] + 'qfep',
       'QCALC'        : Q_DIR['CSB'] + 'qcalc'
      }

HEBBE = {'NODES'      : '1',
         'NTASKS'     : '20',
         'TIME'       : '0-02:00:00',  # d-hh:mm:ss
         'MODULES'    : 'module load GCC/5.4.0-2.26\nmodule load OpenMPI/1.10.3\n', # Add a \n for every added module
         'QDYN'       : 'qdyn=/c3se/users/jwillem/Hebbe/software/qsource/bin/qdyn5p', #fix qdyn= !!!!
         'QFEP'       : 'qdyn=/c3se/users/jwillem/Hebbe/software/qsource/bin/qfep5', #fix qdyn= !!!!!!
         'QPREP'      : '/home/jespers/software/q6/bin/qprep', # NOTE: change to where you are setting up, not where you are running!
         'ACCOUNT'    : 'SNIC2018-2-3'
        }

KEBNE = {'NODES'      : '1',
         'NTASKS'     : '28',
         'TIME'       : '0-04:00:00',  # d-hh:mm:ss
         'MODULES'    : 'module load gompi/2017b\n', # Add a \n for every added module
         'QDYN'       : 'qdyn=/home/w/wije/pfs/software/q/bin/qdynp', #fix qdyn= !!!!!
         'QPREP'      : '/home/jespers/software/q6/bin/qprep', # NOTE: change to where you are setting up, not where you are running!
         'QFEP'       : '/home/w/wije/pfs/software/q/bin/qfep',
         'QCALC'      : '/home/w/wije/pfs/software/q/bin/qcalc',
         'ACCOUNT'    : 'SNIC2018-2-3'
        }

STALLO = {'NODES'      : '1',
         'NTASKS'     : '20',
         'TIME'       : '0-12:00:00',  # d-hh:mm:ss
         'MODULES'    : 'module load impi/2018.1.163-iccifort-2018.1.163-GCC-6.4.0-2.28\n', # Add a \n for every added module
         'QDYN'       : 'qdyn=/home/jespersw/software/Q6/bin/qdynp', #fix qdyn= !!!!!
         'QPREP'      : '/home/apps/q-5.06/qprep', # NOTE: change to where you are setting up, not where you are running!
         'ACCOUNT'    : 'nn4654K'
        }

TETRA = {'NODES'      : '1',
         'NTASKS'     : '32',
         'TIME'       : '0-24:00:00',  # d-hh:mm:ss
         'MODULES'    : '\n', # Add a \n for every added module
         'QDYN'       : 'qdyn=/home/x_wilje/Software/q6/bin/qdynp', #fix qdyn= !!!!!
         'QPREP'      : '/home/apps/q-5.06/qprep', # NOTE: change to where you are setting up, not where you are running!
         'QFEP'       : '/home/x_wilje/Software/q6/bin/qfep',
         'ACCOUNT'    : 'snic2018-2-3'
        }

LOCAL = {'NODES'      : '',
         'NTASKS'     : '',
         'TIME'       : '',  # d-hh:mm:ss
         'MODULES'    : '\n', # Add a \n for every added module
         'QDYN'       : 'qdyn=/Users/willemjespers/Software/Q6/bin/qdyn', #fix qdyn= !!!!!
         'QPREP'      : '/Users/willemjespers/Software/Q6/bin/qprep', # NOTE: change to where you are setting up, not where you are running!
         'QFEP'       : '/Users/willemjespers/Software/Q6/bin/qfep',
         'ACCOUNT'    : ''
        }
