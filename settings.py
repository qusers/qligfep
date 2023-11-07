import os
# Root directory of the setup FEP modules
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

# The directories to the input FF and run related input files are given here
FF_DIR = os.path.join(ROOT_DIR, "FF")
INPUT_DIR = os.path.join(ROOT_DIR, "INPUTS")
# Dicionary of locations of Q executables
Q_DIR = {'CSB':'/home/willem/software/Q/bin/q6/',
         'LOCAL':'/home/yannickrvd/software/Q/bin/q6/',
         'ALICE':'/home/s2904160/data1/projects/pi-gerard/Q/bin/qfep',
        #  'SNELLIUS':'/home/wjespers/software/Q/bin/q6/'
         #'LOCAL':'/Users/willemjespers/Software/Q6/bin/'
        }
BIN = os.path.join(ROOT_DIR, "bin")

# some example schrodinger directory
SCHROD_DIR = '/mnt/c/Program\ Files/Schrodinger2022-2/'

# quick fix to run .exe files on wsl
EXE = '.exe'

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

# CLUSTER INPUTS. To add your own cluster, use the same input as below
# SNELLIUS = {'NODES'        : '1',
#             'NTASKS'       : '16',
#             'TIME'         : '0-12:00:00',  # d-hh:mm:ss
#             'MODULES'      : 'module load 2021\n module load gompi/2021a',
#             'QDYN'         : 'qdyn=' + Q_DIR['SNELLIUS'] + 'qdynp',
#             'QPREP'        : Q_DIR['LOCAL'] + 'qprep',
#             'QFEP'         : Q_DIR['SNELLIUS'] + 'qfep',
#             'QCALC'        : Q_DIR['SNELLIUS'] + 'qcalc'
#       }

ALICE = {'MAINDIR'      : '/home/jespersw/software/q6/',
         'NODES'        : '1',
         'NTASKS'       : '24',
         'TIME'         : '0-3:00:00',  # d-hh:mm:ss
         'MODULES'      : 'module load OpenMPI/3.1.3-GCC-8.2.0-2.31.1',
         'QDYN'         : 'qdyn=/home/jespersw/software/q6/bin/qdynp',
         'QPREP'        : '/home/yannickrvd/software/Q/bin/q6/qprep',
         'QFEP'         : '/home/jespersw/software/q6/bin/qfep',
         'QCALC'        : '/home/jespersw/software/q6/bin/qcalc'
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

UPPMAX = {'NODES'      : '1',
         'NTASKS'     : '20',
         'TIME'       : '0-24:00:00',  # d-hh:mm:ss
         'MODULES'    : 'gcc/9.2.0\nopenmpi/4.0.2\n', # Add a \n for every added module
         'QDYN'       : 'qdyn=/domus/h1/willem/software/q6/bin/qdynp', #fix qdyn= !!!!!
         'QPREP'      : '/home/apps/q-5.06/qprep', # NOTE: change to where you are setting up, not where you are running!
         'QFEP'       : '/domus/h1/willem/software/q6/bin/qfep',
          'ACCOUNT'    : 'snic2018-2-3'
        }

TETRA  = {'NODES'      : '1',
          'NTASKS'     : '32',
          'TIME'       : '0-24:00:00',  # d-hh:mm:ss
          'MODULES'    : '\n', # Add a \n for every added module
          'QDYN'       : 'qdyn=/home/x_wilje/Software/q6/bin/qdynp', #fix qdyn= !!!!!
          'QPREP'      : '/home/jespers/software/q6/bin/qprep', # NOTE: change to where you are setting up, not where you are running!
          'QFEP'       : '/home/x_wilje/Software/q6/bin/qfep',
          'ACCOUNT'    : 'snic2019-2-1'
        }

LOCAL = {'NODES'      : '',
         'NTASKS'     : '',
         'TIME'       : '',  # d-hh:mm:ss
         'MODULES'    : '\n', # Add a \n for every added module
         'QDYN'       : 'qdyn=/home/yannickrvd/software/Q/bin/q6/qdyn', #fix qdyn= !!!!!
         'QPREP'      : '/home/yannickrvd/software/Q/bin/q6/qprep', # NOTE: change to where you are setting up, not where you are running!
         'QFEP'       : '/home/yannickrvd/software/Q/bin/q6/qfep',
         'ACCOUNT'    : ''
        }
