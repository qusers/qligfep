import argparse
import glob
import numpy as np

import functions as f
import settings as s
import IO

class Run(object):
    """
    """
    def __init__(self, LIEdir, *args, **kwargs):
        self.LIEdir = LIEdir
        self.vdw    = []
        self.el     = []
        self.a      = 0.18
        self.b      = 0.33
        self.g      = 0
        
    def read_LIE(self):
        LIEs = sorted(glob.glob(self.LIEdir + '/*/*/*/md_LIE_01.log'))
        for LIE in LIEs:
            with open(LIE) as infile:
                vdw_n   = []
                el_n    = []                
                for line in infile:
                    line = line.split()
                    if len(line) > 3:
                        if line[0] == 'Q-surr.' and line[1] == '1':
                            #print line
                            vdw_n.append(float(line[3]))
                            el_n.append(float(line[4]))
                            
                self.vdw.append(vdw_n)
                self.el.append(el_n)
                
        self.vdw = np.array(self.vdw)
        self.el = np.array(self.el)
        
    def calc_LIE(self):
        avg_vdw = np.mean(self.vdw, axis=0)
        avg_el = np.mean(self.el, axis=0)
        
        vdW_sem = np.nanstd(avg_vdw, ddof =1)/np.sqrt(len(avg_vdw))
        vdW = np.nanmean(avg_vdw)
        
        el_sem = np.nanstd(avg_el, ddof =1)/np.sqrt(len(avg_el))
        el = np.nanmean(avg_el)
        
        print vdW, vdW_sem, el, el_sem 
        dG_LIE = self.a * vdW + self.b * el * self.g
        print dG_LIE
                
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='analyzeLIE',
        version='1.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Analyse LIE == ')

    
    parser.add_argument('-L', '--LIEdir',
                        dest = "LIEdir",
                        required = True,
                        help = "name of LIE directory (FEP_$)")
    
    
    args = parser.parse_args()
    run = Run(LIEdir = args.LIEdir)
    
    run.read_LIE()
    run.calc_LIE()
