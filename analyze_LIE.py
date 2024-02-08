#!/usr/bin/env python

import os
import argparse
import glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from rdkit import Chem
from rdkit.Chem import Fragments
from rdkit.Chem import Draw

import functions as f
import settings as s
import IO

class Run(object):
    """
    Analyses the LIE calculation: determines the beta for the ligand (sdf flie) 
    and reads the Q energies (log files). Then calcuates the average Van der Waals 
    and electrostatic energies, and the free energy. 
    """
    def __init__(self, LIEdir, cluster, *args, **kwargs):
        self.LIEdir  = LIEdir.strip('/')
        self.cluster = cluster
        self.vdw     = []
        self.el      = []
        self.a       = 0.18
        self.b0      = 0.33
        self.beta    = 0
        self.g       = 0
        self.path    = os.getcwd()

    def create_environment(self):
        self.analysisdir = self.LIEdir + '/analysis'
        # Add overwrite function?
        if os.path.isdir(self.analysisdir) != True:
            os.mkdir(self.analysisdir)

    def get_ligand(self):
        ligand = glob.glob('*sdf')[0]
        self.ligand  = self.path+'/'+ligand
        return self.ligand
        
    def read_LIE(self):
        LIEs = {}
        self.MDs = 0
        self.Reps = sorted(glob.glob(self.LIEdir + '/FEP1/*/*'))
        for rep in self.Reps:
            LIEs[rep] = sorted(glob.glob(rep + '/md_LIE_*.log'))
            if len(LIEs[rep]) > self.MDs:
                self.MDs = len(LIEs[rep])
        
        for x,rep in enumerate(LIEs):
            for i,LIE_MD in enumerate(LIEs[rep]):

                with open(LIE_MD) as infile:
                    vdw_n   = []
                    el_n    = []                
                    for line in infile:
                        line = line.split()
                        if len(line) > 3:
                            if line[0] == 'Q-surr.' and line[1] == '1':
                                vdw_n.append(float(line[4]))
                                el_n.append(float(line[3]))

                    if i == 0: # if ../md_LIE_01.log file, append new list to energy array's
                        self.vdw.append(vdw_n)
                        self.el.append(el_n)
                    else: # if non ../md_LIE_01.log file, join list to current replicate list
                        self.vdw[x] = [*self.vdw[x], *vdw_n]
                        self.el[x]  = [*self.el[x], *el_n]

        self.vdw = [x for x in self.vdw if x != []] # discards empty lists of failed reps
        self.el  = [x for x in self.el if x != []]
    
    def Q_test(self, data):     # Function for calculation of outliers
        Q = {3: 0.941, 4: 0.765, 5: 0.642, 6: 0.560, 7: 0.507, 8: 0.437, 9: 0.437, 10: 0.412}
        Q1 = abs(data[0]-data[1])/(data[-1]-data[0])
        Q2 = abs(data[-1]-data[-2])/(data[-1]-data[0])
        if Q1 > Q[len(data)]:
            data.pop(0)
        if Q2 > Q[len(data)]:
            data.pop(-1)
        return data
        
    def calc_LIE(self):
        avg_vdw = []
        avg_el  = []
        
        for rep in self.vdw:
            avg_vdw.append(np.nanmean(rep)) # list of average vdw per rep
            
#        avg_vdw.sort()     # Remove outliers - Turn this part on if you want to remove outliers (vdw part)
#        while True:
#            a = avg_vdw[:]
#            self.Q_test(avg_vdw)
#            if a == avg_vdw:
#                break
        
        vdw = np.nanmean(avg_vdw) # total average vdw
#        vdw = np.nanmedian(avg_vdw) # total median vdw     # Turn on if you want to use the median instead of the average (vdw)
        vdw_sem = np.nanstd(avg_vdw)/np.sqrt(len(avg_vdw))
        
        for rep in self.el:
            avg_el.append(np.nanmean(rep)) # list of average el per rep
        
#        avg_el.sort()      # Remove outliers - Turn this part on if you want to remove outliers (el part)
#        while True:
#            a = avg_el[:]
#            self.Q_test(avg_el)
#            if a == avg_el:
#                break
        
        el = np.nanmean(avg_el) # total average el
#        el = np.nanmedian(avg_el) # total median el        # Turn on if you want to use the median instead of the average (el)
        el_sem = np.nanstd(avg_el)/np.sqrt(len(avg_el)) 
        
        print('vdw: {:.3f} +/- {:.3f} | el: {:.3f} +/- {:.3f} | beta {:.3f}'.format(vdw, vdw_sem, el, el_sem, self.beta))
        dG_LIE = self.a * vdw + self.beta * el + self.g
        dG_sem = self.a * vdw_sem + self.beta * el_sem
        print('dG: {:.3f} +/- {:.3f}'.format(dG_LIE, dG_sem))
        self.energy = {'vdw': vdw, 'el': el}
        return self.energy

    def write_re2pdb(self):
        curdir = os.getcwd()
        os.chdir(self.LIEdir + '/analysis')
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
                pdb_out = 'pdbs/{}_{}'.format(repeat, pdb_out)
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
    
class Beta(object):
    """
    Calculates the beta value for a ligand.
    """    
    def __init__(self, ligand):
        self.ligand    = ligand
        self.mol       = Chem.rdmolfiles.MolFromMolFile(self.ligand, removeHs=False)
        self.FG_Smarts = {                              # non-zero functional groups
            'AL':   '[#6][OX2H]',                               # Alcohols
            'AD1':  '[NX3;H2][CX3](=[OX1])[#6]',                # Primary Amide
            'AN12': '[NX3;H2,H1;!$(NC=O)]',                     # Primary of Secondary Amine
            'CA':   '[CX3](=O)[OX2H1]'                          # Carboxylic Acid
        }
        self.OG_Smarts = {                             # other (zero) groups
            'AD2':  '[NX3;H1][CX3](=[OX1])[#6]',                # Secondary Amide
            'AD3':  '[NX3;H0][CX3](=[OX1])[#6]',                # Tertiary Amide
            'AN3':  '[NX3;H0;!$(NC=[!#6]);!$(NC#[!#6])]',       # Tertiary Amines
            'KT':   '[#6][CX3](=O)[#6]',                        # Ketone
            'ADH':  '[CX3H1](=O)[#6]',                          # Aldehyde
            'TH':   '[#16X2H]',                                 # Thiol
            'ET':   '[OD2]([#6])[#6]',                          # Ether
            'ES':   '[#6][CX3](=O)[OX2H0][#6]',                 # Ester
            'NL':   '[NX1]#[CX2]',                              # Nitrile
            'NO':   '[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]',  # Nitro
            'SU':   '[#16X2H0]'                                 # Sulfide
        }
        self.FG = {}
        self.OG = {}
        
        self.charge = Chem.rdmolops.GetFormalCharge(self.mol)
        if self.charge < 0:
            self.FG['ANI']  = abs(self.charge)
            self.FG['CAT']  = 0
        elif self.charge > 0:
            self.FG['ANI']  = 0
            self.FG['CAT']  = abs(self.charge)
        else:
            self.FG['ANI']  = 0
            self.FG['CAT']  = 0

    def debug(self, FG):
        X = Chem.MolFromSmarts(self.OG_Smarts[FG])
        print(FG, Chem.MolToSmiles(X))
        match = self.mol.GetSubstructMatches(X)
        print(match, len(match))
        
    def Get_FG(self, X):
        FG = Chem.MolFromSmarts(self.FG_Smarts[X])
        match = self.mol.GetSubstructMatches(FG)
        self.FG[X] = len(match)
        return self.FG[X]
    
    def Get_OG(self, X):
        OG = Chem.MolFromSmarts(self.OG_Smarts[X])
        match = self.mol.GetSubstructMatches(OG)
        self.OG[X] = len(match)
        return self.OG[X]

    def beta(self):
        b0  = 0.43
        w_t = sum(self.OG.values()) + self.FG['CA'] + self.FG['AN12'] + self.FG['AD1'] + self.FG['AL'] + 11*(self.FG['ANI'] + self.FG['CAT'])
        db  = {
            'AL':   [-0.06, 1], 
            'AN12': [-0.04, 1],
            'AD1':  [-0.02, 1],
            'CA':   [-0.03, 1],
            'ANI':  [0.02, 11],
            'CAT':  [0.09, 11]
        }
        w_b = 0
        for i in db:
            w_b += db[i][0]*db[i][1]*self.FG[i]
        try:
            beta = b0 + (w_b/w_t)
        except:
            beta = b0
        return beta
    
class Plot(object):
    """
    Plots the energy profile of the LIE simulation.
    Plot 1 shows the average energy of all replicates combined.
    Plot 2 shows the average energy per replicate.
    """    
    def __init__(self, LIEdir, data, *args, **kwargs):
        self.LIEdir  = LIEdir
        self.vdw     = []
        self.el      = []
        self.avg_vdw = []
        self.avg_el  = []
        self.fl_vdw  = []
        self.fl_el   = []
        self.data    = data

    def fl_en_rep(self, en):
        a = np.cumsum(en, axis=1)
        b = np.indices(a.shape)
        b = b[1] + 1
        c = a / b
        return c

    def fl_en_avg(self, en):
        a = np.cumsum(en)
        b = np.indices(a.shape)
        b = b + 1
        c = a / b
        return c
        
    def calc_LIE(self):
        vdw = self.data[0]
        el  = self.data[1]
        self.MDs = self.data[2]

        for rep in vdw: # adds nan to equalize all list lengths, in order to define arrays
            if len(rep) < (self.MDs*10001):
                x = (self.MDs*10001) - len(rep)
                for i in range(x):
                    rep.append(np.nan)
        for rep in el:
            if len(rep) < (self.MDs*10001):
                x = (self.MDs*10001) - len(rep)
                for i in range(x):
                    rep.append(np.nan)
        
        self.vdw = np.array(vdw)
        self.el = np.array(el)
        
        self.x_axis = range(0, len(self.vdw[0]))
        
        self.avg_vdw = np.nanmean(self.vdw, axis=0)
        self.avg_el  = np.nanmean(self.el, axis=0)
        
        self.avg_fl_vdw = self.fl_en_avg(self.avg_vdw)  # FIX: after a replicate has converged, make sure it's endpoint still affects the floating                                                          average afterwards
        self.avg_fl_el = self.fl_en_avg(self.avg_el)    #  in order that ongoing replicates do not affect the floating average to severely
        
        self.fl_vdw = self.fl_en_rep(self.vdw)
        self.fl_el  = self.fl_en_rep(self.el)
        
        avg_vdw = []
        avg_el  = []
        
        for rep in self.vdw:
            avg_vdw.append(np.nanmean(rep)) # list of average vdw per rep
        vdw = np.nanmean(avg_vdw) # total average vdw
        self.vdw_sem = np.nanstd(avg_vdw)/np.sqrt(len(avg_vdw))
        
        for rep in self.el:
            avg_el.append(np.nanmean(rep)) # list of average el per rep
        el = np.nanmean(avg_el) # total average el
        self.el_sem = np.nanstd(avg_el)/np.sqrt(len(avg_el)) 

    def new_plot(self): # plots the floating mean of all individual reps and the total average floating mean
        mpl.style.use('seaborn')
        for i in range(len(self.fl_vdw)):
            plt.figure('QLIE floating mean all reps', figsize=(20,20))
            plt.plot(self.x_axis, self.fl_vdw[i], color='tab:blue', alpha=0.7, linewidth=1)
            plt.plot(self.x_axis, self.fl_el[i], color='tab:orange', alpha=0.7, linewidth=1)
        plt.errorbar(self.x_axis, self.avg_fl_vdw[0], yerr=self.vdw_sem, color='tab:blue', label='vdw', linewidth=2, errorevery=(self.MDs*1000), zorder=3)
        plt.errorbar(self.x_axis, self.avg_fl_el[0], yerr=self.el_sem, color='tab:orange', label='el', linewidth=2, errorevery=(self.MDs*1000), zorder=3)
        plt.xlabel('time (fs)')
        plt.ylabel('Q-surr (kcal/mol)')
        plt.title('Average Energy Development QLIE and floating mean '+self.LIEdir)
        plt.legend()
        plt.savefig(self.LIEdir+'/analysis/dG_reps.png', dpi=300,transparent=True)
        
    def average_plot(self): # plots the total average energy and it's floating mean
        mpl.style.use('seaborn')
        plt.figure('QLIE average energy and floating mean', figsize=(20,20))
        plt.plot(self.x_axis, self.avg_vdw, color='tab:blue', alpha=0.7, linewidth=1)
        plt.plot(self.x_axis, self.avg_el, color='tab:orange', alpha=0.7, linewidth=1)
        plt.errorbar(self.x_axis, self.avg_fl_vdw[0], yerr=self.vdw_sem, color='tab:blue', label='vdw', linewidth=2, errorevery=(self.MDs*1000), zorder=3)
        plt.errorbar(self.x_axis, self.avg_fl_el[0], yerr=self.el_sem, color='tab:orange', label='el', linewidth=2, errorevery=(self.MDs*1000), zorder=3)
        plt.xlabel('time (fs)')
        plt.ylabel('Q-surr (kcal/mol)')
        plt.title('Energy Development QLIE floating mean all reps '+self.LIEdir)
        plt.legend()
        plt.savefig(self.LIEdir+'/analysis/dG_avg.png', dpi=300,transparent=True)

                
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='analyzeLIE',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Analyse LIE == ')
    
    parser.add_argument('-L', '--LIEdir',
                        dest = "LIEdir",
                        required = True,
                        help = "name of LIE directory")
    
    parser.add_argument('-p', '--plot',
                        dest = "plot",
                        required = False,
                        action = 'store_true',
                        help = "include argument to obtain images of energy development")
    
    parser.add_argument('-pdb', '--PDB',
                        dest = "PDB",
                        required = False,
                        default = False,
                        action = 'store_true',
                        help = "Add this argument if you want .pdb files of the trajectory")
    
    parser.add_argument('-C', '--cluster',
                        dest = "cluster",
                        required = False,
                        help = "cluster information")
    
    args = parser.parse_args()
    run = Run(LIEdir  = args.LIEdir,
              PDB     = args.PDB,
              cluster = args.cluster
             )

    run.create_environment()
    run.get_ligand()
    beta = Beta(ligand = run.ligand)
    for FG in beta.FG_Smarts:
        beta.Get_FG(FG)
    for OG in beta.OG_Smarts:
        beta.Get_OG(OG)
    run.beta = beta.beta()
    run.read_LIE()
    run.calc_LIE()
    
    if args.plot == True:
        plot = Plot(data = [run.vdw, run.el, run.MDs],
                    LIEdir = run.LIEdir)
        plot.calc_LIE()
        plot.new_plot()
        plot.average_plot()
    
    if args.PDB == True:
        run.write_re2pdb()