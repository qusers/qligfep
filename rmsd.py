import os
import argparse
import sys
import glob
import numpy as np

import IO

class rmsd(object):
    def __init__(self, LIEdir, rep, option, *args, **kwargs):
        self.LIEdir   = LIEdir.strip('/')
        self.rep      = rep
        self.option   = option
        self.top      = []
        self.protein  = []
        self.ligand   = []
        self.backbone = []
        self.residues = []
        self.dualtop  = self.LIEdir+'/FEP1/298/'+rep+'/dualtop.top'
        self.eq_re    = self.LIEdir+'/FEP1/298/'+rep+'/eq5.re'
        self.eq_top   = self.LIEdir+'/FEP1/298/'+rep+'/eq5.top'
        self.eq_inp   = self.LIEdir+'/inputfiles/eq5re2top.inp'
        self.eq_out   = self.LIEdir+'/inputfiles/eq5re2top.log'
        self.MD       = glob.glob(self.LIEdir+'/FEP1/298/'+rep+'/md_LIE_*.dcd')[-1]
        self.libs     = glob.glob(self.LIEdir + '/inputfiles/*.lib')
        self.inp      = self.LIEdir+'/inputfiles/rmsd_'+option+'.inp'
        self.out      = self.LIEdir+'/FEP1/298/'+rep+'/rmsd_'+option+'.txt'
        
        # Parses the topology pdb file topology list.
        top = LIEdir+'/inputfiles/top_p.pdb'
        with open(top, 'r') as infile:
            for line in infile:
                self.top.append(IO.pdb_parse_in(line))
                
        libs = glob.glob(self.LIEdir + '/inputfiles/*.lib')
        eq_re    = self.LIEdir+'/FEP1/298/'+rep+'/eq5.re'
        eq_top   = self.LIEdir+'/FEP1/298/'+rep+'/eq5.top'
        eq_commands = ['rl '+libs[0], 'rl '+libs[1], 'rt '+self.dualtop, 'rx '+eq_re, 'wt '+eq_top, 'q']
        with open(self.eq_inp, 'w') as outfile:
            for line in eq_commands:
                outfile.write(line+'\n')
        
    
    # Writes the qcalc input file commands with the protein and ligand atom numbers.
    def General(self):
        for atom in self.top:
            if len(atom[0]) != 6:
                continue
            try:
                IO.AA(atom[4])
                self.protein.append(atom[1])
            except:
                if atom[4] == 'LIG':
                    self.ligand.append(atom[1])
                else:
                    continue
        self.commands = [self.eq_top, '1', str(self.protein[0])+' '+str(self.protein[-1]),
                         str(self.ligand[0])+' '+str(self.ligand[-1]), 'end', 'go', self.MD, '.']
        with open(self.inp, 'w') as outfile:
            for line in self.commands:
                outfile.write(line+'\n')
    
    # Writes the qcalc input file commands with the alpha carbon atom numbers.
    def Backbone(self):
        for atom in self.top:
            if atom[2] == 'CA':
                self.backbone.append(atom[1])
        self.commands = [self.eq_top, '1']
        for CA in self.backbone:
            self.commands.append(str(CA))
        self.commands.append('end')
        self.commands.append('go')
        self.commands.append(self.MD)
        self.commands.append('.')
        with open(self.inp, 'w') as outfile:
            for line in self.commands:
                outfile.write(line+'\n')

    # Writes x qcalc input files commands for sets of 10 residues.
    def Residues(self, option):
        # Option 2 for all the residue atoms.
        if option == 2:    
            resnum = 1
            residue = []
            for atom in self.top:
                if atom[6] == resnum:
                    residue.append(atom[1])
                else:
                    self.residues.append(str(residue[0])+' '+str(residue[-1]))
                    residue.clear()
                    residue.append(atom[1])
                    resnum += 1

                if atom[2] == 'OXT':
                    self.residues.append(str(residue[0])+' '+str(residue[-1]))
                    break
        
        # Option 1 for only the residue side chains beta carbons.
        elif option == 1:
            for atom in self.top:
                if atom[2] == 'CB':
                    self.residues.append(str(atom[1]))
        self.commands = []
        
        # Determines the input commands in sets 0f 10 residues.
        x = -1
        for i,resn in enumerate(self.residues, 1):
            if i//10 != x:
                self.commands.append([self.eq_top])
            self.commands[x].append('1')
            self.commands[x].append(resn)
            self.commands[x].append('end')
            if i//10 != x:
                if i == 1:
                    x += 1
                    continue
                self.commands[x].append('go')
                self.commands[x].append(self.MD)
                self.commands[x].append('.')
                x += 1
        self.commands[x].append('go')
        self.commands[x].append(self.MD)
        self.commands[x].append('.')
        
        for i,j in enumerate(self.commands, 1):
            if i <= 9:
                with open(self.LIEdir+'/inputfiles/rmsd_'+self.option+'_00'+str(i)+'.inp', 'w') as outfile:
                    for line in j:
                        outfile.write(line+'\n')
            else:
                with open(self.LIEdir+'/inputfiles/rmsd_'+self.option+'_0'+str(i)+'.inp', 'w') as outfile:
                    for line in j:
                        outfile.write(line+'\n')
        
    # Writes the qcalc input file commands with the ligand atom numbers.
    def Ligand(self):
        for atom in self.top:
            if atom[4] == 'LIG':
                self.ligand.append(atom[1])
                
        self.commands = [self.eq_top, '1', str(self.ligand[0])+' '+str(self.ligand[-1]),
                         'end', 'go', self.MD, '.']
        with open(self.inp, 'w') as outfile:
            for line in self.commands:
                outfile.write(line+'\n')                
    
    def get_top(self, infile, outfile):
        qprep = '/home/koenekoop/software/q6/bin/qprep'
        inp = ' < '+infile
        out = ' > '+outfile
        IO.run_command(qprep, inp + out, string=True)
        
    
    # Executes qcalc and writes the result to the self.out text file.
    def qcalc(self, infile, outfile):
        qcalc = '/home/koenekoop/software/q6/bin/qcalc'
        inp = ' < '+infile
        out = ' > '+outfile
        IO.run_command(qcalc, inp + out, string=True)
    
    # Calculates the average rmsd value with the output.txt file.
    def rmsd(self):
        frames = []
        with open(self.out) as infile:
            for line in infile:
                if line[0:3] == 'LIE':
                    frames.append(float(line[26:34]))
        rmsd = np.average(frames)
        std  = np.std(frames)
        return ['{:.6f}'.format(rmsd), '{:.4f}'.format(std)]
    
    # Calculates the average rmsd value for each residue with the output.txt file.
    def rmsd_residues(self):
        residues = {}
        with open(self.out) as infile:
            for line in infile:
                if line[0:8] == 'Residues':
                    r = int(line[9:12])
                    resn = [x for x in range(int(line[9:12]),int(line[13:16])+1)]
                    for i in resn:
                        residues[i] = []
                if line[0:3] == 'LIE':
                    a = 26
                    b = 34
                    for i in range(r, r+10):
                        try:
                            residues[i].append(float(line[a:b]))
                        except:
                            try:
                                del residues[i]
                            except:
                                pass
                        a += 9
                        b += 9
        
        for i in residues:
            rmsd = np.average(residues[i])
            std  = np.std(residues[i])
            print(str(i)+':', '{:.6f} +/- {:.4f}'.format(rmsd, std))
            

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='RMSD',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Calculates the RMSD of the option requested == ')
    
    parser.add_argument('-L', '--LIEdir',
                        dest = "LIEdir",
                        required = True,
                        help = "name of LIE directory (FEP_$)")
    
    parser.add_argument('-r', '--rep',
                        dest = "rep",
                        required = True,
                        help = "rep to evaluate")
    
    parser.add_argument('-o', '--option',
                        dest = "option",
                        required = True,
                        choices = ['general', 'backbone', 'residues', 'ligand'],
                        help = "RMSA option to calculate")
    
    args = parser.parse_args()
    run = rmsd(LIEdir = args.LIEdir,
               rep    = args.rep,
               option = args.option)
    run.get_top(run.eq_inp, run.eq_out)
    
    if run.option == 'general':
        run.General()
        run.qcalc(run.inp, run.out)
        print(run.option+' rmsd: '+run.rmsd()[0], '+/-', run.rmsd()[1])

    elif run.option == 'backbone':
        run.Backbone()
        run.qcalc(run.inp, run.out)
        print(run.option+' rmsd: '+run.rmsd()[0], '+/-', run.rmsd()[1])

    elif run.option == 'residues':
        while True:
                option = int(input('Side chain CBs (1) or All residue atoms (2)? '))
                if option == 1 or option == 2:
                    break
                else:
                    print("That's not a valid option!")
        run.Residues(option)
        inputs = sorted(glob.glob(run.LIEdir+'/inputfiles/rmsd_'+run.option+'*.inp'))
        for i,infile in enumerate(inputs):
            if i <= 9:
                run.qcalc(infile, run.LIEdir+'/FEP1/298/'+run.rep+'/rmsd_'+run.option+'_00'+str(i)+'.txt')
            else:
                run.qcalc(infile, run.LIEdir+'/FEP1/298/'+run.rep+'/rmsd_'+run.option+'_0'+str(i)+'.txt')
        outfiles = sorted(glob.glob(run.LIEdir+'/FEP1/298/'+run.rep+'/rmsd_'+run.option+'_0*.txt'))
        with open(run.out, 'w') as outfile:
            for i,file in enumerate(outfiles, 1):
                with open(file) as infile:
                    if (i*10-9) <= 9:
                        outfile.write('Residues 00'+str(i*10-9)+'-0'+str(i*10)+'\n')
                    elif (i*10-9) >= 10 and (i*10) <= 99:
                        outfile.write('Residues 0'+str(i*10-9)+'-0'+str(i*10)+'\n')
                    elif (i*10-9) >= 10 and (i*10) == 100:
                        outfile.write('Residues 0'+str(i*10-9)+'-'+str(i*10)+'\n')
                    else:
                        outfile.write('Residues '+str(i*10-9)+'-'+str(i*10)+'\n')
                    outfile.write('----------------------------- Calculation results ----------------------------\n')
                    outfile.write('file                 frame  1: RMSD(A)  2: RMSD(A)  3: RMSD(A)  4: RMSD(A)  5: RMSD(A)  6: RMSD(A)  7: RMSD(A)  8: RMSD(A)  9: RMSD(A) 10: RMSD(A)\n')
                    for line in infile:
                        if line[0:3] == 'LIE':
                            outfile.write(line)
                    outfile.write('\n')
        for file in outfiles:
            os.remove(file)
        print('rmsd per residue:')
        run.rmsd_residues()

    elif run.option == 'ligand':
        run.Ligand()
        run.qcalc(run.inp, run.out)
        print(run.option+' rmsd: '+run.rmsd()[0], '+/-', run.rmsd()[1])
        
    else:
        print('Error: invalid argument --option')
        
    print('outputfile: '+os.getcwd()+'/'+run.out)