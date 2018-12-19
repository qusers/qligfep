import os
# Root directory of the setup FEP modules
ROOT_DIR = os.path.dirname(os.path.realpath(__file__))

# The directories to the input FF and run related input files are given here
FF_DIR = os.path.join(ROOT_DIR, "FF")
INPUT_DIR = os.path.join(ROOT_DIR, "INPUTS")
# Dicionary of locations of Q executables
Q_DIR = {'CLUSTER_NAME':'$QDIR',
         'LOCAL':'$QDIR_LOCAL'
        }
BIN = os.path.join(ROOT_DIR, "bin")
#SCHROD_DIR = '/opt/schrodinger/suites2017-3/'
SCHROD_DIR = '/home/apps/schrodinger2017/'


# CLUSTER INPUTS. To add your own cluster, use the same input as below
CLUSTER_NAME = {'NODES'      : '1',
                'NTASKS'     : '20',
                'TIME'       : '0-03:00:00',  # d-hh:mm:ss
                'MODULES'    : 'module load $SOME_LOCAL_MODULE\n', # Add a \n for every added module
                'QDYN'       : 'qdyn=$QDIR_CLUSTER/Qdyn5p', #fix qdyn= !!!!!
                'QPREP'      : '$QDIR_LOCAL/qprep', # NOTE: change to where you are setting up, not where you are running!
                'ACCOUNT'    : '$SNIC_ACCOUNT'
               }
