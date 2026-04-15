import os
# Root directory of the setup FEP modules
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

# The directories to the input FF and run related input files are given here
FF_DIR = os.path.join(ROOT_DIR, "FF")
INPUT_DIR = os.path.join(ROOT_DIR, "INPUTS")

# Dicionary of locations of Q executables
# NOTE: Q_PATH should be the directory containing the Q binaries (qprep, qdyn, qfep, qcalc)
# Examples:
#   /usr/local/q/bin (for system-wide Q installation)
#   /path/to/qligfep/Q/bin (for local Q in repository)
Q_DIR = {'default': '/workspaces/qligfep/Q/bin/',
        }

BIN = os.path.join(ROOT_DIR, "bin")

# Schrodinger directory (optional for certain forcefields)
SCHROD_DIR = ''

# Default cluster to use when not specified
DEFAULT = 'default'

# CLUSTER INPUTS. To add your own cluster, use the same input as below.
# Each cluster needs an entry in Q_DIR with the path to its Q binaries.
default = {'NODES'        : '1',
       'NTASKS'       : '8',
       'TIME'         : '1-00:00:00',  # d-hh:mm:ss
       'PARTITION'    : '',
       'EXCLUDE'      : '',
       'MODULES'      : '\n',
       'QDYN'         : 'qdyn=' + Q_DIR['default'] + 'qdyn',
       'QPREP'        : Q_DIR['default'] + 'qprep',
       'QFEP'         : Q_DIR['default'] + 'qfep',
       'QCALC'        : Q_DIR['default'] + 'qcalc'
      }
