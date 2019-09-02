# QligFEP and QresFEP 

This collection of python command line functions is designed with the
aim to faciliate a robust and fast setup of FEP calculations for the
software package **Q**. It includes modules which are described in further
detail below:  

- QligFEP.py: module to generate ligand FEP calculations using a
dual topology approach, see Jespers et al. (https://doi.org/10.1186/s13321-019-0348-5).  

- QresFEP.py: module to generate protein FEP calculations using a
single topology approach, see Jespers et al. (https://doi.org/10.1021/acs.jctc.9b00538). 

- opls2Q and charmm2Q.py: tools to translate from OPLS (ffld_server) 
and CHARMM (cgenff) to Q readable file formats.

A few toplevel scripts are included in the scripts folder to faciliate
high throughput setup. Additionally, a tutorials folder is included
with a detailed description of the setup procedure as published in
Jespers et al. (QresFEP/QligFEP). This tutorial includes the generation
of ligand parameters using OPLS, how to prepare a protein system, and
how to run ligand and protein FEP calculations. These examples are 
based on ligand binding of CDk2 inhibitors.

# Installing QligFEP and QresFEP  

- Install a working version of Q, e.g.:

https://github.com/esguerra/Q6

- Clone this repository:  

    git clone https://github.com/qusers/qligfep.git


In settings.py:  

- Change SCHROD_DIR to the Schrodinger location, if you want to be
able to generate OPLS ligand parameters using ffld_server.  

- Change Q_DIR to the location of the q executables. This can be
particularly useful if you use setupFEP from a local machine on
a mounted directory. (In which case, the executables of the preperation
part and running part of Q are at several places).  

- You can add slurm specific parameters in the CLUSTER INPUTS section,
according to the given example.   

## Requirements  
- ffld_server  
- cgenff  
- Protein Preparation Wizard  
- Python2.7.XX  
- Q  

contact:willem.jespers@icm.uu.se

