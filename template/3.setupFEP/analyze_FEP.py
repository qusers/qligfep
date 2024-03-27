#!/usr/bin/env python

import argparse
import glob
import numpy as np
import pandas as pd
import os
import csv

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
    def __init__(self, FEP, color, PDB, cluster, esc, *args, **kwargs):
        self.esc = esc
        self.cluster=cluster
        self.FEP = FEP.strip('/')
        self.energies = {}
        self.FEPstages = []
        FEPfiles = glob.glob(self.FEP + '/inputfiles/FEP*.fep')
        inputs = glob.glob(self.FEP + '/inputfiles/md*.inp')
        inputs = [x for x in inputs if not '_F.inp' in x]
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
                if self.esc:
                    energies = IO.read_qfep_esc(filename)
                else:
                    energies = IO.read_qfep(filename)
            except:
                print("Could not retrieve energies for: " + filename)
                energies = [np.nan, np.nan, np.nan, np.nan, np.nan]
                if not replicate in self.failed:
                    self.failed.append(replicate)
            
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
            for key in sorted(methods[method]):
               # print(method, key, methods[method][key])
                dG_array.append(methods[method][key])
            dG_array = [[float(i) for i in fep_stage] for fep_stage in dG_array]
            dG_array = np.array(dG_array)
            dG = f.avg_sem(dG_array)
            results[method]='{:6.2f}{:6.2f}'.format(*dG)
            
        self.methods = methods
        print(methods)
        
        key_method = []
        dG_values = []

        # Iterate through the dictionary
        for key, value in methods.items():
            key_method.append(key)
            dG_values.append(value['1.QligFEP_CDK2'])
    
        # for method in methods:
            # print(method) 
        # print('crashes  {}'.format(len(self.failed)))
    
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
    
 

    def SEM_calculator(self, output_file, FEP):
        methods = self.methods
        existing_file = os.path.isfile(output_file)
        if existing_file:
            filename_base, file_extension = os.path.splitext(output_file)
            index = 1
            while existing_file:
                new_filename = f"{filename_base}-{index}{file_extension}"
                existing_file = os.path.isfile(new_filename)  # Update the existing_file variable
                index += 1
            output_file = new_filename
            
        with open(output_file, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)

            # Add the FEP directory path as a row at the top (header)
            csv_writer.writerow(['FEP Directory Path'])
            csv_writer.writerow([FEP])

            csv_writer.writerow(['Key', 'Value', 'SEM'])
            for key in methods:
                value = methods[key]['1.QligFEP_CDK2']
                value_without_nan = np.array(value)[~np.isnan(value)]
                SEM = np.std(value_without_nan) / (np.sqrt(len(value_without_nan)))
                csv_writer.writerow([key, value, SEM])
    
    def plot_data(self):
        y_axis = {}
        x_axis = range(0,self.lambda_sum+1)
        avg = []
        for replicate in self.energies:
            for FEPstage in self.FEPstages:
                if not FEPstage in self.energies[replicate] and not replicate in self.failed:
                    self.failed.append(replicate)
        for replicate in self.failed:
            try:
                del self.energies[replicate]
            except:
                continue
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
        
    def write_re2pdb(self):
        curdir = os.getcwd()
        os.chdir(self.FEP + '/analysis')
        if not os.path.exists('pdbs'):
            os.mkdir('pdbs')

        libfiles = glob.glob('../inputfiles/*.lib')
        re_files = glob.glob('../FEP*/*/*/*.re')
        topology = glob.glob('../inputfiles/*.top')[0]
        
        with open('../inputfiles/qprep.inp') as f:
            protlib = f.readline()

        with open('re2pdb.inp', 'w') as outfile:
            outfile.write('{}'.format(protlib))
            
            for libfile in libfiles:
                outfile.write('rl {}\n'.format(libfile))

            outfile.write('rt {}\n'.format(topology))

            for re_file in re_files:
                pdb_out = re_file.split('/')[-1][:-3]
                repeat = '{:02d}'.format(int(re_file.split('/')[3]))
                fep_stage = re_file.split('/')[1]
                pdb_out = 'pdbs/{}_{}_{}'.format(repeat, fep_stage, pdb_out)
                outfile.write('rx {}\n'.format(re_file))
                outfile.write('wp {}.pdb\n'.format(pdb_out))
                outfile.write('y\n')
            
            outfile.write('mask none\n')
            outfile.write('mask not excluded\n')
            outfile.write('wp pdbs/complexnotexcluded.pdb\n')
            outfile.write('y\n')
            
            outfile.write('q\n')

        cluster_options = getattr(s, self.cluster)
        qprep = cluster_options['QPREP']
        options = ' < re2pdb.inp > re2pdb.out'
        # Somehow Q is very annoying with this < > input style so had to implement
        # another function that just calls os.system instead of using the preferred
        # subprocess module....
        IO.run_command(qprep, options, string = True)
            
        os.chdir(curdir)
        
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='protPREP',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Analyse FEP == ')

    
    parser.add_argument('-F', '--FEP',
                        dest = "FEP",
                        required = True,
                        help = "name of FEP directory (FEP_$)")
    
    parser.add_argument('-pdb', '--PDB',
                        dest = "PDB",
                        required = False,
                        default = False,
                        action = 'store_true',
                        help = "Add this argument if you want .pdb files of the trajectory")
    
    parser.add_argument('-c', '--color',
                        dest = "color",
                        required = False,
                        default = 'blue',
                        choices = ['blue', 'red'],
                        help = "color for the plot")
    
    parser.add_argument('-C', '--cluster',
                        dest = "cluster",
                        required = False,
                        help = "cluster information")
    parser.add_argument('-esc', '--end-state-catastrophe',
                        dest = "esc",
                        required = False,
                        help = "Add this argument in case you have singularities in the final lambda windows resulting in only nan values")
    
    args = parser.parse_args()
    run = Run(FEP = args.FEP,
              color = args.color,
              cluster = args.cluster,
              PDB = args.PDB,
              esc = args.esc
             )
    
    run.create_environment()
    run.read_FEPs()
    run.SEM_calculator('dG.csv', args.FEP)
    run.read_mdlog()
    
    if plot == True:
        run.plot_data()
        
    if args.PDB == True:
        run.write_re2pdb()