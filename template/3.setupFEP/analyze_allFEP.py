import glob
import numpy as np
import math

def get_results():
    results = {}
    residues = []
        
    for filename in glob.glob('1.*/*/FEP_*/FEP*/*/*'):
        line = filename.split('/')
        run = line[0]
        FEP = line[1]
        name = line[2] + '_' + line[1].split('.')[1]
        replicate = line[-1]
        ERROR = False
     
        with open(filename + '/qfep.out') as infile:
            block = 0
            for line in infile:
                line = line.split()
                if len(line) > 3:
                    if line[0] == 'ERROR:' or line[1] == 'ERROR:' or line[1] == 'Failed':
                        ERROR = True
                    
                    if line[3] == 'Free':
                        block = 1

                    if line[3] == 'Termodynamic':
                        block = 2

                    if line[3] == 'Overlap':
                        block = 3

                    if line[3] == 'BAR':
                        block = 4

                    if line[3] == 'Reaction':
                        block = 0

                if len(line) > 1:
                    if block == 1:
                        if line[0] == '1.000000':
                            dGr = line[4]

                        elif line[0] == '0.000000':
                            dGf = line[2]
                            
                            if line[5] == '-Infinity':
                                dG = np.nan
                                
                            else:
                                dG = float(line[5])

                    if block == 2 and line[0] == '0.000000':
                        #dG_TI = line[2]
                        dG_TI = np.nan

                    if block == 3 and line[0] == '0.000000':
                        if line[2] == '-Infinity':
                            dG_overlap = np.nan
                            
                        else:
                            dG_overlap = float(line[2])

                    if block == 4 and line[0] == '0.000000':
                        if line[2] == '-Infinity':
                            dG_BAR = np.nan
                        else:
                            dG_BAR = float(line[2])
        
        if ERROR != True:
            data = [name, replicate, dG, dG_overlap, dG_BAR]
            
        else:
            data = [name, replicate, np.nan, np.nan, np.nan]

        if name in results:
            #if len(results[name]) < 10:    # temp fix
            results[name].append(data)
        else:
            results[name] = [data]
    return results

def calc_sum_error(data):
    cnt = 0
    for line in data:
        if np.isnan(line) == True:
            continue
        else:
            cnt = cnt + 1
    #data = np.array(data)
    try:
        data = np.array(data)
        mean = np.nanmean(data)
        sem = np.nanstd(np.array(data), ddof =1)/np.sqrt(cnt) #Check if this is correct!
        return mean, sem
    except:
        return None, None
    #return mean, sem

def calc_ddG(raw_data):
    with open('results.txt', 'w') as outfile:
        outfile.write('FEP          Zwanzig     error    OS    error    BAR      error\n')
        dG = {}
        for name in raw_data:
            scores = []
            a = name.split('_')
            FEP = a[0] + '_' + a[1] 
            print(FEP)
            system = a[-1]
            scoring = [[], [], []]
            for line in raw_data[name]:
                scoring[0].append(line[2])
                scoring[1].append(line[3])
                scoring[2].append(line[4])

            for data in scoring:
                scores.append(calc_sum_error(data))

            scores = [system,scores]
            if FEP not in dG:
                dG[FEP] = [scores]

            else:            
                dG[FEP].append(scores)

        for key in dG:
            data = dG[key]
            for system in data:
                if system[0] == 'protein':
                    Zwanzig_avg_prot= system[1][0][0]
                    Zwanzig_err_prot= system[1][0][1]
                    OS_avg_prot= system[1][1][0]
                    OS_err_prot= system[1][1][1]
                    BAR_avg_prot= system[1][2][0]
                    BAR_err_prot= system[1][2][1]

                if system[0] == 'water':
                    Zwanzig_avg_wat= system[1][0][0]
                    Zwanzig_err_wat= system[1][0][1]
                    OS_avg_wat= system[1][1][0]
                    OS_err_wat= system[1][1][1]
                    BAR_avg_wat= system[1][2][0]
                    BAR_err_wat= system[1][2][1]

            try:
                ddG_Zwanzig = Zwanzig_avg_prot - Zwanzig_avg_wat
                ddG_OS = OS_avg_prot - OS_avg_wat
                ddG_BAR = BAR_avg_prot - BAR_avg_wat
                ddG_Zwanzig_error = (Zwanzig_err_prot + Zwanzig_err_wat)/math.sqrt(2)
                ddG_OS_error = (OS_err_prot + OS_err_wat)/math.sqrt(2)
                ddG_BAR_error = (BAR_err_prot + BAR_err_wat)/math.sqrt(2)
            except:
                ddG_Zwanzig = 0
                ddG_Zwanzig_error = 0
                ddG_OS = 0
                ddG_OS_error = 0
                ddG_BAR = 0
                ddG_BAR_error = 0 
            outfile.write('{:30}{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}{:8.2f}\n'.format(key, 
                                                                               ddG_Zwanzig, 
                                                                                ddG_Zwanzig_error,
                                                                                ddG_OS,
                                                                                ddG_OS_error,
                                                                                ddG_BAR,
                                                                                ddG_BAR_error
                                                                                )) 
            
            

raw_data = get_results()
calc_ddG(raw_data)
