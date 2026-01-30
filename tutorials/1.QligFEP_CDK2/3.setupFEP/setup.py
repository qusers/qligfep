import os
import shutil
import glob

systems = ['protein', 'water']
cnt = 0
# Change this to where you installed QligFEP
setupFEP = 'python $QligFEPDIR/QligFEP.py'
cysbond = ' '

for system in systems:
    cnt += 1
    directory = str(cnt) + '.' + system
    os.mkdir(directory)
    with open('pairs.txt') as infile:

        for line in infile:
            line = line.split()
            mol1 = line[0]
            mol2 = line[1]
            if system == 'water':

                call = setupFEP + ' -l1 ' + mol1 + ' -l2 ' + mol2 + ' -FF OPLS2015 -s water -c CSB -r 22 -l 0.5'

            if system == 'protein':
                call = setupFEP + ' -l1 ' + mol1 + ' -l2 ' + mol2 + ' -FF OPLS2015 -s protein -c CSB -r 22 -l 0.5' 

            src = 'FEP_' + mol1 + '-' + mol2
            dst = directory + '/FEP_' + mol1 + '-' + mol2
            os.system(call)
            shutil.move(src, dst)
