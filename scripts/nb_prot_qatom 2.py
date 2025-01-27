import glob
import sys
import os
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
import IO
import settings as s

qcalc =  s.KEBNE['QCALC']
q_res = ['258', '259']

replicates = 10
replacements = {'PROT_START':'1',
                'PROT_END':'257'
               }

def write_qcalc(dcd):
    out = 'qcalc_{}.inp'.format(q_res[i])
    with open(s.INPUT_DIR+ '/qcalc.inp') as infile,          \
         open(out, 'w') as outfile:
        for line in infile:
            line = IO.replace(line, replacements)
            if line.rstrip() == 'TRAJECTORIES':
                for trajectory in dcd:
                    trajectory = trajectory + '\n'
                    outfile.write(trajectory)
                continue
                    
            outfile.write(line)

def run_qcalc(q_resid):
    for i in range(1, replicates + 1):
        curdir = os.getcwd()
        os.chdir('FEP1/298/{}'.format(i))
        if q_resid == 0:
            dcd = glob.glob('md_1000_0000.dcd')
        elif q_resid == 1:
            dcd = glob.glob('md_0000_1000.dcd')
        write_qcalc(dcd)
        os.system('{} < qcalc_{}.inp > qcalc_{}.log'.format(qcalc,
                                                            q_res[q_resid],
                                                            q_resid))
        os.chdir(curdir)

def get_results(q_resid):
    data = {'vdw':[],'el':[]}
    for i in range(1, replicates + 1):
        j = i - 1
        data['vdw'].append([])
        data['el'].append([])
        with open('FEP1/298/{}/qcalc_{}.log'.format(i,q_resid)) as infile:
            block = 0            
            for line in infile:
                line = line.split()
                if len(line) < 2:
                    continue
                
                if line[0] == 'Residue':
                    block = 1
                    continue
                
                if block == 1:
                    data['vdw'][j].append(float(line[1]))
                    data['el'][j].append(float(line[2]))

        if len(data['vdw'][j]) == 0:
            for i in range(0, int(replacements['PROT_END'])):
                data['vdw'][j].append(np.nan)
                data['el'][j].append(np.nan)
                        
    vdw = np.array(data['vdw'])
    el = np.array(data['el'])
    
    vdw_mean = np.nanmean(vdw, axis=0)
    el_mean = np.nanmean(el, axis=0)
    
    vdw_err = np.nanstd(vdw, axis=0)/np.sqrt(10.0)
    el_err = np.nanstd(el, axis=0)/np.sqrt(10.0)
    
    
    return vdw_mean, vdw_err, el_mean, el_err

#replacements['Q_ATOMS'] = ' ' + q_res[0]
#run_qcalc(0)
data_1 = get_results(0)

#replacements['Q_ATOMS'] = ' ' + q_res[1]
#run_qcalc(1)
data_2 = get_results(1)

el_diff = data_1[2] - data_2[2]
for i in range(0, len(el_diff)):
    if el_diff[i] > 1.0:
        print i, el_diff[i]
