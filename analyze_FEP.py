import argparse
import glob
import numpy as np
import os

import functions as f
import settings as s
import IO

try:
    import matplotlib
    matplotlib.use('Agg')    
    import matplotlib.pyplot as plt
    plot = True
except:
    print('cannot import matplotlib, skipping plot generation')
    plot = False

class Run(object):
    """
    """
    def __init__(self, FEP, color, *args, **kwargs):
        self.FEP = FEP.strip('/')
        self.energies = {}
        self.FEPstages = []
        FEPfiles = glob.glob(self.FEP + '/inputfiles/FEP*.fep')
        inputs = glob.glob(self.FEP + '/inputfiles/md*.inp')
        FEPfiles.sort()
        self.failed = []
        for FEPfile in FEPfiles:
            FEPstage = FEPfile.split('/')[-1]
            FEPstage = FEPstage.split('.')[0]
            self.FEPstages.append(FEPstage)
            
        self.lambda_sum = len(FEPfiles) * (len(inputs)-1)
        
        colors = {'blue':['navy','lightblue'],
                  'red' :['darkred','mistyrose']
                 }
        
        self.color = colors[color]
        
    def create_environment(self):
        self.analysisdir = self.FEP + '/analysis'
        # Add overwrite function?
        if os.path.isdir(self.analysisdir) != True:
            os.mkdir(self.analysisdir)
    
    def read_FEPs(self):
        methods_list = ['dG', 'dGf', 'dGr', 'dGos', 'dGbar']
        methods = {'dG'     : {},
                   'dGf'    : {},
                   'dGr'    : {},
                   'dGos'   : {},
                   'dGbar' :  {}
                  }
        results = {}
        out = []
        
        FEPs = sorted(glob.glob(self.FEP + '/*/*/*/qfep.out'))
        for filename in FEPs:
            i = -1
            file_parse = filename.split('/')
            FEP = file_parse[1]
            temperature = file_parse[2]
            replicate = file_parse[3]
            
            try:
                energies = IO.read_qfep(filename)
            except:
                print("Could not retrieve energies for: " + filename)
                energies = [np.nan, np.nan, np.nan, np.nan, np.nan]
                self.failed.append(replicate)
            #try:
            #    energies = IO.read_qfep(filename)
            #except:
            #    print "Could not retrieve energies for: " + filename
            #    energies = [np.nan, np.nan, np.nan, np.nan, np.nan]

            for key in methods_list:
                i += 1
                try:
                    methods[key][FEP].append(energies[i])
                except:
                    methods[key][FEP] = [energies[i]]
                    
            # Construct for the energy figure
            if not replicate in self.energies:
                self.energies[replicate] = {}
                
            self.energies[replicate][FEP] = IO.read_qfep_verbose(filename)
            
        for method in methods:
            dG_array = []
            for key in methods[method]:
                #print method, key, methods[method][key]
                dG_array.append(methods[method][key])
            dG_array = np.array(dG_array)
            dG_array = dG_array.astype(np.float)
            dG = f.avg_sem(dG_array)
            results[method]='{:6.2f}{:6.2f}'.format(*dG)
            
        for method in methods_list:
            out.append(results[method])

        print(self.FEP, '{} {} {} {} {}'.format(*out))
        
    def read_mdlog(self):
        mapping = {}
        cnt = -1
        # Add temperature variable later
        md_files = glob.glob(self.FEP + '/FEP*/*/*/md*.log')        
        md_files.sort()
        md_ref = glob.glob(self.FEP + '/inputfiles/md*.inp')
        windows = len(glob.glob(self.FEP + '/inputfiles/md*.inp')) - 1
        stages = len(glob.glob(self.FEP + '/inputfiles/FEP*.fep'))
        for ref in md_ref:
            w = ref.split('/')[-1].split('_')[2].split('.')[0]
            cnt += 1
            mapping[w]=cnt

        for md_file in md_files:
            stage = md_file.split('/')[1][-1]
            l = md_file.split('/')[-1].split('_')[2].split('.')[0]
            offset = (int(stage) - 1) * int(windows)
            cumulative_l = (mapping[l] + offset)
            (cumulative_l)
            
    def plot_data(self):
        y_axis = {}
        x_axis = range(0,self.lambda_sum+1)
        avg = []
        for replicate in self.failed:
            del self.energies[replicate]
        for replicate in self.energies:
            y_axis[replicate] = [0]
            dG = 0
            for FEPstage in self.FEPstages:
                for energy in self.energies[replicate][FEPstage][0][1:]:
                    energy = dG + energy
                    y_axis[replicate].append(energy)
                dG+=self.energies[replicate][FEPstage][0][-1]
        
        for y in y_axis:
            for i,energy in enumerate(y_axis[y]):
                if len(avg) < self.lambda_sum + 1:
                    avg.append(energy)
                else:
                    avg[i] += energy

            plt.plot(x_axis,y_axis[y],color=self.color[1])
        y_avg = [x / len(y_axis) for x in avg]
        plt.plot(x_axis,y_avg,color=self.color[0])
        axes = plt.gca()
        axes.set_xlim([0,self.lambda_sum])
        plt.xlabel(r'cumulative $\lambda$', fontsize=18)
        plt.ylabel(r'$\Delta$G (kcal/mol)', fontsize=16)        
        plt.savefig(self.analysisdir+'/dG.png',dpi=300,transparent=True)
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='protPREP',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Analyse FEP == ')

    
    parser.add_argument('-F', '--FEP',
                        dest = "FEP",
                        required = True,
                        help = "name of FEP directory (FEP_$)")
    
    parser.add_argument('-c', '--color',
                        dest = "color",
                        required = False,
                        default = 'blue',
                        choices = ['blue', 'red'],
                        help = "name of FEP directory (FEP_$)")
    
    
    args = parser.parse_args()
    run = Run(FEP = args.FEP,
              color = args.color
             )
    
    run.create_environment()
    run.read_FEPs()
    run.read_mdlog()
    if plot == True:
        run.plot_data()
