import glob
import statistics
# import pandas
import argparse
import os

slurmfile = glob.glob('slurm*.out')

with open(slurmfile[0]) as slurm:
    for line in slurm:
        if "STOP qprep ended normally" in line:
            continue
        else:
            print(line.strip("\n"))