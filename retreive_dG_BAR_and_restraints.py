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
        self.lib_for_df = {}
        self.dG_BAR = {}

    def return_restraints(self):
        repfolders = sorted(glob.glob(self.FEP + '/*/'))
        directories = [repfolders[0] + entry for entry in os.listdir(repfolders[0]) if os.path.isdir(repfolders[0] + entry)]
        for folder in directories:
            logfilesrep = sorted(glob.glob(folder + '/md*.log'))
            rep_restraints = []
            for file in logfilesrep:
                with open(file) as logfile:
                    for line in logfile:
                        if line[:10] == "restraints":
                            solute_restraint = line.split(" ")[-1].strip('\n')
                            rep_restraints.append(float(solute_restraint))
                if rep_restraints != []:
                    self.lib_for_df[file.split("/")[-2]] = rep_restraints
                    print_statement = ""
                else:
                    self.lib_for_df[file.split("/")[-2]] = 'nan'
                    print_statement = 'Could not retrieve restraints for:' + folder
            if print_statement == "":
                continue
            else:
                print(print_statement)

    def return_dG_BAR_method(self):
        repfolders = sorted(glob.glob(self.FEP + '/*/'))
        directories = [repfolders[0] + entry for entry in os.listdir(repfolders[0]) if os.path.isdir(repfolders[0] + entry)]
        for folder in directories:    
            qfepfilesrep = sorted(glob.glob(folder + '/qfep*.out'))
            for file in qfepfilesrep:
                # print (qfepfilesrep)
                with open(file, "r") as energies:
                    for line in energies:
                        pass
                    dGBAR = line
                    dGBAR = dGBAR.strip('\n').split(" ")[-1]
                    try:
                        dGBAR = float(dGBAR)
                        self.dG_BAR[file.split("/")[-2]] = dGBAR
                    except ValueError:
                        self.dG_BAR[file.split("/")[-2]] = "nan"

    def print_results(self):
        averages_restraints = []
        dG_BAR_values = []
        for rep_key in self.lib_for_df:
            if self.lib_for_df[rep_key] == [] or "nan" in self.lib_for_df[rep_key] or self.dG_BAR[rep_key] == "nan":
                averages_restraints.append('nan')
                dG_BAR_values.append('nan')
            else:
                average = statistics.mean(self.lib_for_df[rep_key])
                averages_restraints.append(str(round(average,3)))
                dG_BAR_values.append(str(self.dG_BAR[rep_key]))
        print(', '.join(dG_BAR_values))      
        print(', '.join(averages_restraints))

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
    run.return_restraints()
    run.print_results()