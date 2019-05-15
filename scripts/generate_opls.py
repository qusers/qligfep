import glob
from subprocess import check_output
import shlex
import os

curdir = os.getcwd()
os.chdir(curdir)
opls2Q = '/home/jespers/software/qligfep-p3/qligfep/opls2Q.py'
generate = 'python ' + opls2Q


for pdb in glob.glob('*.pdb'):
    name = pdb.split('.')[0]
    options = ' -l ' + name + ' -FF OPLS2015 -o Q -m -vs'
    args = shlex.split(generate + options)
    out = check_output(args)
