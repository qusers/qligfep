import glob
from subprocess import check_output
import shlex
import os

cur_dir = os.getcwd()
script_dir = os.path.dirname(os.path.dirname(__file__))

os.chdir(cur_dir)
opls2Q = f'{script_dir}/opls2Q.py'

generate = 'python3 ' + opls2Q

for pdb in glob.glob('*.pdb'):
    name = pdb.split('.')[0]
    options = ' -l ' + name + ' -FF OPLS2015 -m'
    args = shlex.split(generate + options)
    out = check_output(args)
