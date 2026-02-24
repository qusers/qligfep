import os
# Root directory of the setup FEP modules
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

# The directories to the input FF and run related input files are given here
FF_DIR = os.path.join(ROOT_DIR, "FF")
INPUT_DIR = os.path.join(ROOT_DIR, "INPUTS")
# Dicionary of locations of Q executables
Q_DIR = {'$DEFAULT': '$PATH/bin/',
        }
BIN = os.path.join(ROOT_DIR, "bin")
SCHROD_DIR = '$SCHROD/'

DEFAULT = '$DEFAULT'

# CLUSTER INPUTS. To add your own cluster, use the same input as below
$DEFAULT = {'NODES'        : '1',
       'NTASKS'       : '8',
       'TIME'         : '1-00:00:00',  # d-hh:mm:ss
       'PARTITION'    : '',
       'EXCLUDE'      : '',
       'MODULES'      : '\n',
       'QDYN'         : 'qdyn=' + Q_DIR['$DEFAULT'] + 'qdynp',
       'QPREP'        : Q_DIR['$DEFAULT'] + 'qprep',
       'QFEP'         : Q_DIR['$DEFAULT'] + 'qfep',
       'QCALC'        : Q_DIR['$DEFAULT'] + 'qcalc'
      }