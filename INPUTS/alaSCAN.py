import shutil
import shlex
import os
from subprocess import check_output

mutations = [
                MUTATIONS
            ]
SYSTEMS
systems.append('apo')
QLIGFEP
executable = 'python {}/QresFEP.py '.format(qligfep)
options = {
           '-f'    : 'OPLS2015',       # Forcefield
           '-l'     : '1',             # starting point
           '-T'     : '298',           # Temperature
           '-r'     : '10',            # Total repeats
           '-s'     : 'linear',        # Sampling type
           '-w'     : '20',            # Number of windows
           '-C'     : 'KEBNE',         # Cluster to run on
           '-mc'    : 'A',             # Chain of the mutations
          }

rootdir = os.getcwd()
aladir = rootdir + '/alaSCAN'

for system in systems:
    for mutation in mutations:
        options['-S'] = 'protein'
        options['-m'] = mutation
        if system != 'apo':
            options['-c'] = system
        # Run commands and save to directory
        option_list = ' '.join(['{} {}'.format(k,v) for k,v in options.items()])
        args = shlex.split(executable + option_list)
        out = check_output(args)
        fepdir = 'FEP_{}'.format(mutation)
        tgtdir = aladir + '/' + system + '/' + fepdir
        print(fepdir, tgtdir)
        shutil.move(fepdir, tgtdir)
