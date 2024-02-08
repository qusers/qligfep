import glob
import os, sys

path = os.getcwd()
software= path+'/../../../software'

files = glob.glob('HM_????.pdb') + glob.glob('hoh.pdb') + glob.glob('ions_local.pdb')
if len(files) == 2:
    os.system(software+'/pymemdyn/pymemdyn_2.5ns -p '+files[0]+' -w hoh')
elif len(files) == 3:
    os.system(software+'/pymemdyn/pymemdyn_2.5ns -p '+files[0]+' -w hoh -i ions_local')
else:
    print('Error: insufficent amount of files in folder '+XXXX)