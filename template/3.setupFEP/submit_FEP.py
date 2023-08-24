#!/usr/bin/env python3
import os

# Get all folders starting with "FEP"
folders = [folder for folder in os.listdir('.') if folder.startswith('FEP') and os.path.isdir(folder)]

# Loop through folders
for folder in folders:
    print(f"Entering folder: {folder}")
    os.chdir(folder)
    print(f"Running FEP_submit.sh in {folder}")
    os.system("sbatch FEP_submit.sh")
    # Return to the original directory
    os.chdir('..')
