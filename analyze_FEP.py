import argparse
import glob
import numpy as np

import functions as f
import settings as s
import IO

class Run(object):
    """
    """
    def __init__(self, FEP, *args, **kwargs):
        self.FEP = FEP
        
        
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
                if FEP == 'FEP2':
                    energies = [0.00,0.00,0.00,0.00,0.00]
                else:
                    print "Could not retrieve energies for: " + filename
                    energies = [np.nan, np.nan, np.nan, np.nan, np.nan]
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

        print self.FEP, '{} {} {} {} {}'.format(*out)
                
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='protPREP',
        version='1.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Analyse FEP == ')

    
    parser.add_argument('-F', '--FEP',
                        dest = "FEP",
                        required = True,
                        help = "name of FEP directory (FEP_$)")
    
    
    args = parser.parse_args()
    run = Run(FEP = args.FEP)
    
    run.read_FEPs()
