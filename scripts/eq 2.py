import mdtraj as md
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from glob import glob
from scipy import stats

from scipy.optimize import curve_fit
import glob
import os
import sys

inputfiles = sys.argv[1]

top_pdb = f'{inputfiles}/complexnotexcluded.pdb'

dcds = sorted(glob.glob('eq*.dcd'))[4:]
print(dcds)
print(dcds)
print(top_pdb)

traj = []
traj += [md.load(dcd, top=top_pdb) for dcd in dcds]
traj = md.join(traj)

top = traj[0].topology
bb = top.select('protein and name CA')

rmsd = md.rmsd(traj, traj, atom_indices=bb, frame=0) * 10
print('rmsd:', np.mean(rmsd))
print('std:', np.std(rmsd))

U_pot = []
logs = sorted(glob.glob('eq*.log'))[4:]
print(logs)
for log in logs:
    with open(log, 'r') as file:
        for line in file:
            if line.startswith('SUM'):
                U_pot.append(float(line.split()[2]))
        file.close()

print('Upot:', np.mean(U_pot))
print('std:', np.std(U_pot))

U_pot = np.array(U_pot)
N = len(U_pot)
with open('block_averaging.txt', 'w') as f:
    for n in range(1,N+1):
        M = N / n
        blocks = [U_pot[i:i+n] for i in range(0, N, n)]
        block_means = [np.mean(block) for block in blocks]
        if len(block_means) <= 1:
            continue
        blocks_std = np.std(block_means)
        blocks_sem = blocks_std / np.sqrt(M)
        f.write(f'{n} {blocks_sem}\n')
    f.close()