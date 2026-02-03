# QligFEP v2.0 and QresFEP v2.0

This collection of python command line functions is designed with the
aim to faciliate a robust and fast setup of FEP calculations for the
software package **Q**. These modules use python 3, python 2 is no 
longer supported, and an old version of the code using python 2
is now only available in the python2 branch.

This pacakge includes at the moment four main modules:  
- QligFEP.py: module to generate ligand FEP calculations using a
dual topology approach, 
see Jespers et al. (https://doi.org/10.1186/s13321-019-0348-5).

- QresFEP.py: module to generate protein FEP calculations using either a
single topology approach, 
see Jespers et al. (https://doi.org/10.1021/acs.jctc.9b00538),
or using a dual topology approach,
see Koenekoop et al. (https://doi.org/10.1038/s42004-025-01771-0).

- QLIE.py: module to generate ligand LIE calculations.

Future versions will include several translation tools for new forcefields 
(at the moment we support opls, charmm, amber and openFF).

A few toplevel scripts are included in the scripts folder to faciliate
high throughput setup. Additionally, a tutorials folder is included, with one tutorial
showing a detailed description of the setup procedure for ligand mutations as published in
Jespers et al. (QligFEP). This tutorial includes the generation
of ligand parameters using OPLS, how to prepare a protein system, and
how to run ligand FEP calculations. These examples are 
based on ligand binding of CDK2 inhibitors.
Another tutorial shows a detailed description of the setup procedure for residue mutations 
as published in Koenekoop et al. (QresFEP). This tutorial includes the generation
of mutant amino acids, how to prepare a protein system, and
how to run residue FEP calculations. These examples are 
based on thermal stability shifts of T4L mutants.

# Installing QligFEP and QresFEP  

- Install a working version of Q, e.g.: https://github.com/esguerra/Q6
- Clone this repository:

        git clone https://github.com/qusers/qligfep.git

- Initialize the repository:

        bash qligfep_init.sh

    This will add qligfep repo executables to your $PATH and update the settings.py script with:
    - Q_PATH:   the path to your Q installation
    - SCHROD:   the path to your Schrödinger installation
    - DEFAULT:  your working HPC cluster (or localhost)
    
    You can modify you setting.py script further manually to add more clusters or cluster options.

## Requirements
**Schrödinger**
- PyMOL
- Protein Preparation Wizard  
- ffld_server

**Software**
- Q
- cgenff
- Python3.XX


**Python packages**
- mdtraj

contact: w.jespers@rug.nl

