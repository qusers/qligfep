import os
import sys
import re

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
import IO

path = os.getcwd()
pymol_file = 'mutagenesis.pml'

try:
    protein = sys.argv[1]
    mutation_file = sys.argv[2]
except (IndexError):
    print("Usage: {} <structure.pdb> <mutations.txt>".format(sys.argv[0]))
    sys.exit(1)

mutations = {}

with open(mutation_file, 'r') as infile:
    for fep in infile:
        fep = (re.split('(\d+)', fep.strip('\n')))
        if len(fep) == 3:
            if not fep[1] in mutations:
                if len(fep[0]) == 1:
                    mutations[fep[1]] = [[IO.AA(fep[0]), IO.AA(fep[2])]]
                else:
                    mutations[fep[1]] = [[fep[0], fep[2]]]
            else:
                if len(fep[0]) == 1:
                    mutations[fep[1]].append([IO.AA(fep[0]), IO.AA(fep[2])])
                else:
                    mutations[fep[1]].append([fep[0], fep[2]])

with open(pymol_file, 'w') as pml:
    for feps in mutations:
        for i,fep in enumerate(mutations[feps]):
            pml.write("""
reinitialize
cmd.load('{}')
        
cmd.wizard('mutagenesis')
cmd.do('refresh_wizard')
        
cmd.get_wizard().set_mode('{}')
cmd.get_wizard().do_select('resi {}')
cmd.get_wizard().apply()

cd {}
save {}{}.pdb, resi {}
            """.format(protein, \
                    mutations[feps][i][1], \
                    feps, \
                    path, \
                    mutations[feps][i][1], \
                    feps, \
                    feps
                   ))
    pml.write("cmd.quit()")
sys.exit(0)
