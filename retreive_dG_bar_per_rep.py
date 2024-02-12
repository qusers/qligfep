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

    def return_restraints(self):
        md_files = sorted(glob.glob(self.FEP + '/md*.log'))
        qfepfile = sorted(glob.glob(self.FEP + '/qfep*.out'))
        solute_restraints = []
        shell_restraints = []
        total_restraints = []
        replicate = self.FEP.split[-1]
        for file in md_files:
            with open(file) as logfile:
                for line in logfile:
                    if line[:10] == "restraints":
                        splitted_list = line.split(" ")
                        no_empty_strings = list(filter(None, splitted_list))
                        solute_restraint = no_empty_strings[-1].strip('\n')
                        shell_restraint = no_empty_strings[-2]
                        total_restraint = no_empty_strings[1]
                        solute_restraints.append(float(solute_restraint))
                        shell_restraints.append(float(shell_restraint))
                        total_restraints.append(float(total_restraint))
        with open(qfepfile[0]) as qfepfile_out:
            for line in energies:
                pass
            dGBAR = line
            dGBAR = dGBAR.strip('\n').split(" ")[-1]
            try:
                dGBAR = float(dGBAR)
            except ValueError:
                dGBAR = "nan"

        average_solute = statistics.mean(solute_restraints)
        average_shell = statistics.mean(shell_restraints)
        average_total = statistics.mean(total_restraints)
        list_for_df = [dG_BAR, average_solute, average_shell, average_total]
        print(average_solute)
        # df = pd.read_csv('input.csv')
        # df[replicate] = list_for_df

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