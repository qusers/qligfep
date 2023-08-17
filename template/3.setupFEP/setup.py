import os
import shutil
import glob
import re

# Get the path to the '1.ligprep' directory
current_folder = os.path.dirname(os.path.abspath(__file__))
ligandprep_folder = os.path.join(current_folder, '..', '1.ligprep')

# Check if the '1.ligprep' directory exists
if not os.path.exists(ligandprep_folder):
    print("Error: '1.ligprep' directory not found.")
else:
    # Import ligand files needed for QligFEP
    lig_files = glob.glob(os.path.join(ligandprep_folder, '*.pdb')) + \
                glob.glob(os.path.join(ligandprep_folder, '*.lib')) + \
                glob.glob(os.path.join(ligandprep_folder, '*.prm'))

    for file in lig_files:
        shutil.copy(file, '.')

# Get the path to the '2.protprep' directory
protprep_folder = os.path.join(current_folder, '..', '2.protprep')

# Check if the '2.protprep' directory exists
if not os.path.exists(protprep_folder):
    print("Error: '2.protprep' directory not found.")
    print("Have you used --noclean when calling qligfep/protprep.py?")
else:
    # Import 'protein.pdb'
    protein_pdb = os.path.join(protprep_folder, 'protein.pdb')
    if os.path.exists(protein_pdb):
        shutil.copy(protein_pdb, '.')

    # Import 'water.pdb'
    water_pdb = os.path.join(protprep_folder, 'water.pdb')
    if os.path.exists(water_pdb):
        shutil.copy(water_pdb, '.')

# Name for folders. Can be changed.
systems = ['protein-TETRA', 'water-TETRA']
cnt = 0

# Change this to where you installed QligFEP
setupFEP = 'python /home/USER/QLIGFEP/qligfep/QligFEP.py'
windows = '100'

# Open qprep.inp where cysbond to add are specified.
if os.path.isfile(protprep_folder+'/qprep.inp'):
    with open(protprep_folder+'/qprep.inp', 'r') as file:
        file_content = file.read()

    # Use regular expressions to find lines with 'addbond' and extract the numbers.
    bond_lines = re.findall(r'addbond (\d+):SG (\d+):SG', file_content)
    bonds = {}
    for idx, (atom1, atom2) in enumerate(bond_lines, 1):
        bond_name = f"bond-{idx}"
        bonds[bond_name] = [int(atom1), int(atom2)]

    cys = []  
    for x, i in bonds.values(): 
        with open('protein.pdb', 'r') as file:
            for line in file:
                column = line.split()
                if column[0] == 'ATOM' and column[2] == 'SG':
                    if column[4] == str(x):
                        atom_1= column[1]
                        cys.append(atom_1)
                    if column[4] == str(i):
                        atom_2 = column[1]
                        cys.append(atom_2)
    pairs = [cys[i] + '_' + cys[i+1] for i in range(0, len(cys), 2)] # Add '_' between atoms and match style required by QligFEP
    cysbond =  ','.join(pairs)
    print('cysbond added: ', cysbond) 
    for system in systems:
        cnt += 1
        directory = str(cnt) + '.' + system

        # Check if the directory already exists and create a new one with a different name if needed
        new_directory = directory
        suffix = 2
        while os.path.exists(new_directory):
            new_directory = f"{directory}-{suffix}"
            suffix += 1

        os.mkdir(new_directory)

        # Call QligFEP 
        with open('pairs.txt') as infile:
            for line in infile:
                line = line.split()
                mol1 = line[0]
                mol2 = line[1]
                if system == systems[1]:
                    call = setupFEP + ' -l1 ' + mol1 + ' -l2 ' + mol2 + ' -FF OPLS2015 -s water -S sigmoidal -c TETRA -r 25 -l 1  -w' + windows
                if system == systems[0]:
                    call = setupFEP + ' -l1 ' + mol1 + ' -l2 ' + mol2 + ' -FF OPLS2015 -b' + cysbond + ' -S sigmoidal -s protein -c TETRA -r 25 -l 1 -w' + windows
                src = 'FEP_' + mol1 + '-' + mol2
                dst = os.path.join(new_directory, 'FEP_' + mol1 + '-' + mol2)
                os.system(call)
                shutil.move(src, dst)
else:
    print('Error: qprep.inp file not found in 2.protprep folder. Have you used --noclean when calling protprep.py?')
