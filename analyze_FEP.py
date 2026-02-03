#!/usr/bin/env python

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
    plot = False
except:
    print('cannot import matplotlib, skipping plot generation')
    plot = False

class Run(object):
    """
    Store script arguments and set global variables
    """
    def __init__(self, FEP, temp, cluster, color, PDB, esc, start, *args, **kwargs):
        self.temp = temp
        self.cluster = cluster
        self.esc = esc
        self.FEP = FEP.strip('/')
        self.FromGly = False
        self.start = start
        if self.FEP[4:7] == 'GLY':
            self.FromGly = True
        #self.FromGly = False
        FEPfiles = sorted(glob.glob(self.FEP + '/inputfiles/*.fep'))
        self.FEPstages = [f'FEP{x}' for x in range(1, len(FEPfiles)+1)]
        
        inputs = glob.glob(self.FEP + '/inputfiles/md*.inp')
        self.failed = []
        self.energies = {}

        self.lambda_sum = len(FEPfiles) * (len(inputs)-1)
        
        colors = {'blue':['navy','lightblue'],
                  'red' :['darkred','mistyrose']
                 }
        
        self.color = colors[color]
        
    def create_environment(self):
        """
        Creates analysis/ directory for generated plots and structure files
        """
        self.analysisdir = self.FEP + '/analysis'
        if os.path.isdir(self.analysisdir) != True:
            os.mkdir(self.analysisdir)
    
    def read_FEPs(self):
        """
        Reads the qfep output files (qfep.out) and prints the energies per replicate
        """
        methods_list = ['dG', 'dGf', 'dGr', 'dGbar']
        methods = {'dG'     : {},
                   'dGf'    : {},
                   'dGr'    : {},
                   'dGbar' :  {}
                  }
        results = {}
        out = []
        
        # Loops over the FEPstages, temperatures and replicates
        FEPs = sorted(glob.glob(f'{self.FEP}/FEP*/{self.temp}/*/qfep.out'))
        for filename in FEPs:
            file_parse = filename.split('/')
            FEP = file_parse[1]
            temperature = file_parse[2]
            replicate = file_parse[3]
            
            # In case of end-state catastrophe with single FEPstage, use read_qfep_esc(), otherwise default read_qfep()
            try:
                if self.esc:
                    energies = IO.read_qfep_esc(filename)
                else:
                    energies = IO.read_qfep(filename, self.FromGly, self.start, FEP)
            # Append failed replicates and set their energies to NaN
            except:
                print("Could not retrieve energies for: " + filename)
                energies = [np.nan, np.nan, np.nan, np.nan]
                if not replicate in self.failed:
                    self.failed.append(replicate)
            
            # Store the energies from read_qfep in the methods dictionary
            for i, method in enumerate(methods_list):
                try:
                    methods[method][FEP].append(energies[i])
                except:
                    methods[method][FEP] = [energies[i]]
                    
            # Construct for the energy plot
            if not replicate in self.energies:
                self.energies[replicate] = {}
                
            self.energies[replicate][FEP] = IO.read_qfep_verbose(filename)
        
        # print the energies per method, FEPstage, and replicate, and calculate average total FEP energy and s.e.m.
        for method in methods:
            dG_array = []
            energies = '{0:6.2f} {1:6.2f} {2:6.2f} {3:6.2f} {4:6.2f} {5:6.2f} {6:6.2f} {7:6.2f} {8:6.2f} {9:6.2f}'
            #energies = '{0:6.2f} {1:6.2f} {2:6.2f} {3:6.2f} {4:6.2f}'
            for key in sorted(methods[method]):
                print(f'{method:5} {key}', energies.format(*methods[method][key]))
                dG_array.append(methods[method][key])
            dG_array = [[float(i) for i in fep_stage] for fep_stage in dG_array]
            dG_array = np.array(dG_array)
            dG = f.avg_sem(dG_array)
            results[method]='{:6.2f}{:6.2f}'.format(*dG)

        # Print average total FEP energy and s.e.m.
        for method in methods_list:
            out.append(results[method])
        print(self.FEP, '{} {} {} {}'.format(*out))
        print('crashes  {}'.format(len(self.failed)))
        
    def read_mdlog(self):
        """
        Read the qdyn md.log files to get plotting details
        """
        mapping = {}
        cnt = -1
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
        """
        Plot the energy progression during the full FEP
        """
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
        """
        Write out pdb structures using qdyn .re files as input 
        """
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
        IO.run_command(qprep, options, string = True)
            
        os.chdir(curdir)
        
            
if __name__ == "__main__":
    """
    Parser for user input arguments
    """
    parser = argparse.ArgumentParser(
        prog='analyze_FEP.py',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '===== Analyse FEP output files from Q; Extracts free energy change (ΔG) for a FEP mutation simulation output directory (FEP_{WT}2{MUT}) =====')
    
    # Required command line arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument('-F', '--FEP',
                          dest = "FEP",
                          required = True,
                          help = "Specify the FEP directory (FEP_WT2MUT)")

    required.add_argument('-T', '--temp',
                         dest = "temp",
                         required = True,
                         help = "Specificy for which temperature to analyze")

    required.add_argument('-l', '--start',
                          dest = "start",
                          default = '1',
                          choices = ['1', '0.5', '0'],
                          help = "Specify starting point of FEP simulations; either start (1), middle (0.5) or end (0) lambda point")
    
    # Optional command line arguments
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('-C', '--cluster',
                          dest = "cluster",
                          required = False,
                          default = s.DEFAULT,
                          help = "Use to specify HPC cluster of analysis")

    optional.add_argument('-pdb', '--PDB',
                          dest = "PDB",
                          required = False,
                          default = False,
                          action = 'store_true',
                          help = "Use to generate PDB files of the trajectory")

    optional.add_argument('-c', '--color',
                          dest = "color",
                          required = False,
                          default = 'blue',
                          choices = ['blue', 'red'],
                          help = "Use to specify the color for the energy propagation plot")

    optional.add_argument('-esc', '--end-state-catastrophe',
                          dest = "esc",
                          required = False,
                          help = "Use in case of singularities in the final lambda windows resultin in NaN energy values")
    
    args = parser.parse_args()
    run = Run(FEP = args.FEP,
              temp = args.temp,
              color = args.color,
              cluster = args.cluster,
              PDB = args.PDB,
              esc = args.esc,
              start = args.start)
    
    run.create_environment()
    run.read_FEPs()
    
    if plot == True:
        run.read_mdlog()
        run.plot_data()    
    if args.PDB == True:
        run.write_re2pdb()
