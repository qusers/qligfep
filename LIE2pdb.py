import argparse
import glob
import numpy as np
import os

import functions as f
import settings as s
import IO

class Run(object):
    """
    """
    def __init__(self, LIEdir, *args, **kwargs):
        self.LIEdir = LIEdir
        self.FEP = LIEdir.strip('/')
        self.cluster='CSB'
        self.path   = os.getcwd()
        
        
    def create_environment(self):
        self.analysisdir = self.FEP + '/analysis'
        # Add overwrite function?
        if os.path.isdir(self.analysisdir) != True:
            os.mkdir(self.analysisdir)
            
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
                
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='analyzeLIE',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Analyse LIE == ')

    
    parser.add_argument('-L', '--LIEdir',
                        dest = "LIEdir",
                        required = True,
                        help = "name of LIE directory (FEP_$)")
    
    
    args = parser.parse_args()
    run = Run(LIEdir = args.LIEdir)
    run.create_environment()
    run.write_re2pdb()
