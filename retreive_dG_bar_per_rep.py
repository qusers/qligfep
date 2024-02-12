import glob
import statistics
import pandas as pd
import argparse
import os
import shutil

class Run(object):
    """
    """
    def __init__(self, FEP, *args, **kwargs):
        self.FEP = FEP.strip('/')

    def return_restraints(self):
        md_files = sorted(glob.glob(self.FEP + '/md*.log'))
        runfolder = str(self.FEP)
        qfepfile = sorted(glob.glob(self.FEP + '/qfep*.out')) 
        solute_restraints = []
        shell_restraints = []
        total_restraints = []
        replicate = runfolder.split("/")[-1]
        FEP_path = []
        path_seperation = "/"
        for path_pieces in runfolder.split("/"):
            FEP_path.append(path_pieces)
            if "FEP" in path_pieces:
                break
        FEP_path_total = path_seperation.join(FEP_path)
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
            for line in qfepfile_out:
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
        list_for_df = [dGBAR, average_solute, average_shell, average_total]

        if replicate == "1":
            shutil.copyfile('/gpfs/scratch1/nodespecific/int6/yricky/template_outputfile.csv', '/gpfs/scratch1/nodespecific/int6/yricky/' + FEP_path_total)
            df = pd.read_csv('/gpfs/scratch1/nodespecific/int6/yricky/' + FEP_path_total + "/template_outputfile.csv")
            df[replicate] = list_for_df
        else:
            df = pd.read_csv('/gpfs/scratch1/nodespecific/int6/yricky/' + FEP_path_total + "/template_outputfile.csv")
            df[replicate] = list_for_df
            df.to_csv(''/gpfs/scratch1/nodespecific/int6/yricky/' + FEP_path_total + "/template_outputfile.csv"', index=False)

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
    run.return_restraints()
