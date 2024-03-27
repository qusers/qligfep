import os
import IO
import functions
import settings

# can just be done with os.system(executable + options) probably
qprep = '/home/yannickrvd/software/q6/bin/qprep'
options = ' < qprep.inp > qprep.out'
IO.run_command(qprep, options, string = True)