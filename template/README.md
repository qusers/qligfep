Before starting:
To make your workflow easier, copy the `template` folder in the repository. This contains everything you will use for the common usecase for QligFEP. 
Locate the .pdb files for your ligands in `1.ligprep` folder and your_protein.pdb in the `2.protprep` folder. 

# 1. Ligand Preparation
We initiate the process with ligands generated in Maestro, typically provided as PDB files. However, to conduct MD simulations, we require specific parameters that are not available initially. Hence, the ligands undergo a crucial step known as ligand preparation.


The folder contains pre-generated 3D-coordinates of the ligands, which are not generated by setupFEP. Therefore, it is our responsibility to ensure that these coordinates are reliable. It is essential to bear in mind that FEP handles small changes in the system effectively, but reliable results are obtained only when the phase space sufficiently overlaps. The ligands in this folder share a similar core region, with variations in small substituents being investigated. For AMBER parameters, please refer to qtools:


### Generate the files using a script from the script folder. 
To generate the necessary files using a the script located in the script folder, follow these steps:
1. Move to 1.ligprep folder from the template. If you're not using the template provided, create a folder to store your ligands and call it `1.ligprep`:
* `mkdir 1.ligprep`
2. Navigate to the newly created directory:
* `cd 1.ligprep`
3. Place all your ligand.pdb files inside this folder. If your files are in your local machine and you want to use a Cluster, you can use the following command:
* `scp /path/to/local/file username@cluster_ip:/path/on/cluster/`
4. Finally, run the following command to execute the script:
* `python $qligfep/scripts/generate_opls.py`.
This script will produce the necessary .lib, .prm, and .pdb files required for the third stage.
The .log files are generated by ffld_server.

Important!
  * Make sure that you've exported qligfep to your .bashrc file (or any other of the common Unix shells) so you can call qligfep from a different folder by using `$qligfep`. 
If you haven't, you can do it following these steps:
    1. `cd` to your home folder
    2. `nano ~/.bashrc`
    3. Add this to the end of the file: `export qligfep=/home/USERNAME/QLIGFEP/qligfep`
    - Keep in mind that `/home/USERNAME/QLIGFEP/qligfep` should be replaced with the PATH to your qligfep folder. If you're not sure of what the path to it is, go to your qligfep folder and run `pwd` to find out. Then, simply replace `/home/USERNAME/QLIGFEP/qligfep` with the returned path from pwd. 
    5. Use `source ~/.bashrc` to apply the changes
  

# 2. Protprep
The second step involves preparing the protein for simulations under spherical boundary conditions (SBC). Please note that this script does not assign protonation states, so a preparation step should be executed beforehand (e.g., using Maestro's Protein Preparation Wizard).

The center of the sphere can be determined based on a residue in the protein or, as in this case, by explicitly providing the center of geometry (COG) coordinates from a ligand. For example, to obtain the COG coordinates used in this example, run:

`python "$qligfep/scripts/COG.py" -i ../1.ligprep/17.pdb`

! Make sure to **add -i before the pdb file**. 
This will return the coordinates [0.535:26.772:8.819], which can be directly put in protprep.py:

`python "$qligfep/protprep.py" -p 1h1s_PrepWiz.pdb -r 22 -c [0.535:26.772:8.819] -w -P CSB`

After this, you should have three new files in 2.protprep: protPREP.log  protein.pdb  water.pdb.

**Note: Prior to executing the above command, ensure that you have navigated to the folder containing the "1h1s_PrepWiz.pdb" file.**

Furthermore, it is essential to bear in mind that the ligands under investigation share a similar core region, with only minor variations in small substituents. During previous phases, these ligands were superimposed, likely using Maestro. As a result, their centers of geometry are nearly identical. Therefore, we need to perform this process for just one ligand, saving time and effort.

# 3.setupFEP

In the next stage we will prepare the input files for the actual simulations. This stage uses $qligfep/QligFEP.py, you can use the -h flag to get more specifics on all the input variables. This folder includes a useful top level script, generate.py, that takes a .txt file as input that specifies the ligand pairs that will be used in this simulation. You can simply run:

`python setup.py`

And this should result in a 1.protein and 2.water folder, containing the two systems necessary to perform the calculation. Simply go into the folder (1.PROTEIN OR 2.WATER) and run

`./FEP_submit.sh`

This script will either submit your jobs to a slurm queue (see the setup section in $setupFEP/README.md for a more detailed description). Alternatively the jobs can be run on your local machine.

A typical procedure is to put the two systems in a top layer folder to easily track running calculations (e.g. 1.TESTRUN) in the case of this example. Of course, these toplayer scripts can be easily adjusted for any other use (e.g. to calculate solvation free energies, you add a water and vacuum leg, instead of protein water.


# 4. Analyze

Run `python analyze_FEP.py -F tutorials/1.QligFEP_CDK2/3.setupFEP/1.protein/FEP_22-17 -C CSB` from the `QLIGFEP/qligfep` folder. 
Otherwise it won't work because functions.py cannot be imported. 

Additionally, you can try using collect_dG.py to analyze several folders at once. See the instructions provided on how to use this script below the page. The output of the analysis will contain dG, dGf, dGr, dGos, and dGbar, representing different ways of calculating dG. For most cases, dGbar is used. There will be 10 values corresponding to 10 replicas (runs), and some values may be "nan," which is normal. Thanks to changes in this forked repository, SEM is calculated only with existing values (nan values are dropped).

To change the number of runs, modify the FEP_submit.sh file located in the path: /home/USER/QLIGFEP/qligfep/tutorials/1.QligFEP_CDK2/3.setupFEP/1.protein/FEP_17-22


## How to Use collect_dG.py

`collect_dG.py` is a python script that allows you to call analyze_FEP.py for different folders (usually 1.protein/FEP_pair and 2.water/FEP_pair) with a single command  and stores the results in different .csv files so it's easier to work with the results later on.

To use collect_dG.py, you need to provide at least 4 arguments when calling it on the command line:
- collect_dG.py
- folder_name (the folder where the .csv files will be stored)
- path_to_folder1
- path_to_folder2

and is called from the qligfep folder like this: `python collect_dG.py folder_name path_to_folder1 path_to_folder2`

#### Example for tutorial.
After going trough steps 1, 2, 3 and from `qligfep`folder write
`python collect_dG.py dG-FEP_17-22_2 tutorials/1.QligFEP_CDK2/3.setupFEP/1.protein/FEP_17-22 tutorials/1.QligFEP_CDK2/3.setupFEP/2.water/FEP_17-22`

Now, there should be 2 .csv files dG.csv and dG-1.csv inside the dG-FEP_17-22 folder. Each .csv file contains the path to the folder from which the dG values were calculated along with the calculation results. In later versions, it is expected to create a single .csv file with all the results (problems with merging the files are being addressed).
