import glob
import statistics
# import pandas
import argparse
import os

class Run(object):
    """
    """
    def __init__(self, FEP, *args, **kwargs):
        self.FEP = FEP.strip('/')
        self.dG_BAR = []
        self.sampling = []

    def return_dG_BAR_method(self):
        repfolders = sorted(glob.glob(self.FEP + '/*/'))
        directories = [repfolders[0] + entry for entry in os.listdir(repfolders[0]) if os.path.isdir(repfolders[0] + entry)]
        for folder in directories:    
            qfepfilesrep = sorted(glob.glob(folder + '/qfep*.out'))
            for file in qfepfilesrep:
                # print(file)
                with open(file, "r") as energies:
                    start = 0
                    dgBAR_values = []
                    for line in energies:
                        if "BAR Bennet" in line:
                            start = 1
                        elif start == 1:
                            if line != "\n":
                                dGBAR = line.strip('\n')
                                # print(dGBAR)
                                # no_empty_string = [for numbers in dGBAR]
                                # print(list(filter(None,dGBAR.split(" "))))
                                dGBAR = list(filter(None,dGBAR.split(" ")))[1]
                                # print(dGBAR)
                                dgBAR_values.append(float(dGBAR))
                    # print(dgBAR_va
                    if dgBAR_values != []:
                        self.dG_BAR.append(dgBAR_values)

    def retreive_sampling(self):
        repfolders = sorted(glob.glob(self.FEP + '/*/'))
        directories = [repfolders[0] + entry for entry in os.listdir(repfolders[0]) if os.path.isdir(repfolders[0] + entry)]
        for folder in directories:    
            qfepfilesrep = sorted(glob.glob(folder + '/qfep*.out'))
            for file in qfepfilesrep:
                # print(file)
                with open(file, "r") as energies:
                    start = 0
                    for line in energies:
                        if "BAR Bennet" in line:
                            start = 1
                        elif start == 1:
                            if line != "\n":
                                dGBAR = line.strip('\n')
                                dGBAR = list(filter(None,dGBAR.split(" ")))[0]
                                # print(dGBAR)
                                self.sampling.append(str(dGBAR))
                    break

    def print_results(self):
        # print(self.sampling)
        deltaGs = zip(*self.dG_BAR)
        stdevs = []
        averages = []
        for dGs in reversed(list(deltaGs)):
            dG_strings = [str(dG) for dG in dGs]
            print(', '.join(dG_strings))  
            stdevs.append(str(statistics.stdev(dGs)))
            averages.append(str(statistics.mean(dGs)))
        sampling_stdev = zip(reversed(self.sampling),stdevs,averages)
        for samp_std in (sampling_stdev):
            # print(samp_std)
            print(', '.join(samp_std))   
           
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='retreive_restraints',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = '       == Analyse FEP == ')
    
    parser.add_argument('-F', '--FEP',
                        dest = "FEP",
                        required = True,
                        help = "name of FEP directory (FEP_$)")
    
    args = parser.parse_args()
    run = Run(FEP = args.FEP)

    run.return_dG_BAR_method()
    run.retreive_sampling()
    run.print_results()
    