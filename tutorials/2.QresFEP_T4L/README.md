# QresFEP Tutorial

This tutorial is a step-by-step guide to setup, run and analyze QresFEP calculations, 
as reported in Koenekoop et al. (https://www.nature.com/articles/s42004-025-01771-0).

The workflow consist of the following 4 steps:

1. Prepare the protein to run under spherical boundary conditions with Q
2. Generate mutant amino acid structures
3. Generate all the necessary FEP simulation input files
4. Analyse the FEP simulation output files

This folder contains all the required files to follow these steps. We will now proceed with the first step.

IMPORTANT: Make sure to execute the root repository initialization script before starting this tutorial:

    [Navigate to your /home/user/software/qligfep directory and execute...]
    bash qligfep_init.sh

# 1. Prepare the protein

The subject of this tutorial is the bacteriophage T4 lysozyme (T4L), an often used protein for testing and benchmarking studies, as it was also one of the protein systems used to benchmark QresFEP.

Here we will use the T4L structure of PDB entry [2LZM](https://www.rcsb.org/structure/2LZM), a high resolution structure (1.7 Å) which was not only used in our benchmark but also by others.

A PDB file contains information and coordinates for all atoms in the protein structure (ATOM lines), and possibly also of several other moieties (so called co-crystallizing agents), 
e.g. water molecules, ligands, and/or ions (HETATM lines).
Check the raw T4L PDB file *2LZM.pdb* that is directly downloaded from the RCSB Protein Data Bank [here](https://files.rcsb.org/download/2LZM.pdb).
The raw PDB file usually contains formatting which is not compatible with the forcefield formating used by our MD simulation software Q 
(see details on the MD engine Q on the related [publication](https://www.sciencedirect.com/science/article/pii/S1093326399000121?via%3Dihub). For this, we need to preform some PDB file preparation steps, which typically include:

i) Removing redundant parts of the raw PDB file that should not be present in the structure subject to simulation

ii) Adding missing parts to our structure (certain amino acid side chains,
or even whole amino acid loops from low-resolution regions might be missing). 

This can be handled with several solutions - the one we regularly use makes use of Schrödinger's [Protein Preparation Wizard software](https://www.schrodinger.com/life-science/learn/white-papers/protein-preparation-workflow/). 
In addition to the previous, this solution will add hydrogen atoms - including optimized protonation states for ionizable residues based on pKa predictions - and to any existing structural water molecules - all
based on predicted hydrogen bonding networks. The output structure file *2LZM_-_hbond-opt.pdb* from the PrepWizard will still require some reformatting, as for example all histidine residues will be named as HIS,
yet the forcefield we will use, [OPLS-AA/M](https://zarbi.chem.yale.edu/oplsaam.html), formats the different protonation states of histidine as either HID, HIE, or HIP. In some cases, we might also want to remove some
C-terminal (OXT) or N-terminal (H1, H2) atoms. The final "prepared" structure file (in our example *2LZM_prep.pdb*) is now ready to be treated by the fist script of the QresFEP package, `protprep.py`. This will call the `qprep` program within the Q suite to solvate our system and prepare for running under spherical boundary conditions. Type:

    protprep.py -h

to run this command and see the associated help message, and get familiar with its required command-line options. First of all, with the -p PROT argument we have to indicate a protein PDB file, and with the -r SPHERERADIUS argument we indicate
the radius of our intended sphere (defaults to 25 Å). With the -c SPHETECENTER we can indicate where to center the sphere; either on a residue (use RESN:$ where $ is the residue number), on a specific atom (use ATN:$), or by using
explicit coordinates (use X:Y:Z). In case we choose to center on a residue, we also need to indicate the chain of that residue by using the -mc MUTCHAIN argument. Then with the -f FORCEFIELD argument we can indicate which forcefield
we want to use (defaults to OPLSAAM). We also have to option to not include any crystal structure waters in our structure by turning the -w NOWATER argument on (this is strongly discouraged in most cases). Finally, with the --noclean
flag we can prevent intermediate `qprep` files from being deleted, and with the -P PREPLOCATION flag we can specify which Q executables should be used (as indicated in settings.py - will use its indicated default).
Now that we know what all options are for, prepare our T4L structure by running the following command:

    protprep.py -p 2LZM_prep.pdb -r 25 -c RESN:39 -mc A -f OPLSAAM

The following files should have been generated:

- *protein.pdb*   (contains the protein structure)
- *protPREP.log*  (log file of protprep.py, will be used to extract information by QresFEP.py)
- *water.pdb*     (contains the water molecules)

# 2. Generate mutant amino acid structures

Since this tutorial deals with amino acid mutations, we need to generate mutants. QresFEP works by requesting PDB files of only the mutant amino acid(s) - which need to be spacially aligned to be superpositioned on conformation coordinates 
of the wild-type amino acid. This can be achieved by several ways, and feel free to use any molecular visualisation software of preferance, but we strongly recommend to use [PyMOL](https://www.pymol.org/) which works fairly easily as a
command-line tool as well. This utility allows for automation using PyMOL, and by using the `mutagenesis.py` script, we can do exactly that.

    python mutagenesis.py <structure.pdb> <mutations.txt>

The `mutagenesis.py` script takes two arguments: a PDB structure file - in our case *2LZM_prep.pdb*; and a text file with a list of mutations. Take a look at either *mutations_neutral.txt* or *mutations_charged.txt* - there are 10 random mutations
taken from the T4L benchmark studies for charge-maintaining mutations (Chem Comm 2025) and charge-changing mutations (to be published). Choose whichever option you prefer, and run the following command:

    python mutagenesis.py 2LZM_prep.pdb mutations_neutral.txt

The following file should have been generated:

- mutagenesis.pml   (PyMOL file to generate mutant amino acids PDB files)

With the generated PyMOL run file, we can now run it to generate our mutant PDB files. Execute the following command:

    pymol -c mutagenesis.pml

The following files should have been generated:

- *ALA39.pdb*     *ALA44.pdb*       *ALA120.pdb*      *GLY25.pdb*       *GLY55.pdb*      *GLY105.pdb*      *LEU153.pdb*      *SER26.pdb*       *THR3.pdb*        *THR75.pdb*

Now that the mutant PDB files are generated, there is one more issue that needs be fixed. The output PDB file format by PyMOL uses a different atom name ordering. Take for example a look at *ALA39.pdb*, which should looks like this:

    ATOM      1  N   ALA A  39      38.865  31.620  20.999  1.00  0.00           N
    ATOM      2  CA  ALA A  39      39.713  30.684  21.730  1.00  0.00           C
    ATOM      3  C   ALA A  39      39.134  30.378  23.091  1.00  0.00           C
    ATOM      4  O   ALA A  39      39.394  29.424  23.739  1.00  0.00           O
    ATOM      5  CB  ALA A  39      41.123  31.295  21.807  1.00  0.00           C
    ATOM      6  H   ALA A  39      39.275  32.454  20.602  1.00  0.00           H
    ATOM      7  HA  ALA A  39      39.755  29.732  21.169  1.00  0.00           H
    ATOM      8 1HB  ALA A  39      41.547  31.477  20.801  1.00  0.00           H
    ATOM      9 2HB  ALA A  39      41.129  32.265  22.340  1.00  0.00           H
    ATOM     10 3HB  ALA A  39      41.832  30.630  22.334  1.00  0.00           H
    TER
    END

The side-chain hyrdogens are named 1HB, 2HB and 3HB, whereas our forcefield library uses the atom names HB1, HB2, and HB3. You can check other mutant PDB files as well to see more side-chain hydrogen names ordered differently.
To fix this, we just need to run a simple sed command, which should correct the atom name ordering:

    sed -i -E 's/([1-3])(H[A-Z]) / \2\1/g; s/([1-3])(H[A-Z][0-9])/\2\1/g' *.pdb     (this command corrects all side-chain hydrogens)
    sed 's/ HA / HA2/g' -i GLY*.pdb                                                 (this command corrects Glycine HA to HA2)

Take a look again at *ALA39.pdb*, and if executed correctly, it should now look as follows, end we can move on to the next step:

    ATOM      1  N   ALA A  39      38.865  31.620  20.999  1.00  0.00           N
    ATOM      2  CA  ALA A  39      39.713  30.684  21.730  1.00  0.00           C
    ATOM      3  C   ALA A  39      39.134  30.378  23.091  1.00  0.00           C
    ATOM      4  O   ALA A  39      39.394  29.424  23.739  1.00  0.00           O
    ATOM      5  CB  ALA A  39      41.123  31.295  21.807  1.00  0.00           C
    ATOM      6  H   ALA A  39      39.275  32.454  20.602  1.00  0.00           H
    ATOM      7  HA  ALA A  39      39.755  29.732  21.169  1.00  0.00           H
    ATOM      8  HB1 ALA A  39      41.547  31.477  20.801  1.00  0.00           H
    ATOM      9  HB2 ALA A  39      41.129  32.265  22.340  1.00  0.00           H
    ATOM     10  HB3 ALA A  39      41.832  30.630  22.334  1.00  0.00           H
    TER
    END

# 3. Generate FEP simulation input files

Now that we have our prepared protein structure and mutant amino acids generated, we can use QresFEP to generate our FEP simulation input files.

    QresFEP.py -h

Run this command to see the help message for QresFEP, and which command-line options it requires. Here is a brief description:

- First of all, we need to indicate a mutation with the -m MUTATION argument, which is given as e.g. LEU39ALA (the first neutral mutation).
In this case two things are vital: 
 1. protprep was centered on RESN:39 and 
 2. *ALA39.pdb* is present in the folder. 
- Then, the original PDB chain needs to be indicated again with the -mc MUTCHAIN argument. 
- Next, we need to indicate if we
want to perform the simulation in the protein system, or as a tripeptide (which serves as the reference state - to complete a thermodynamic cycle). We do this with the -S SYSTEM argument. 
 1. In case we choose a tripeptide system, we should
also indicate which type of flanking residues we prefer in our tripeptide using the -t TRIPEPTIDE flag; this can either be alanine (option A), glycine (option G), no flanking residues (option X - effectively a monopeptide), or
the natural sequence neighbouring amino acids (option Z). [Beware: if running charge-changing mutations, option Z results in less stable results.]
- By turning on the -d DUAL flag, the dual/hybrid topology FEP protocol will be utilised (check reference [Koenekoop et al., 2025](https://www.nature.com/articles/s42004-025-01771-0); this is the recommended protocol, and also allows all amino acids mutations. Note that, if this flag is not specifically activated, the single topology FEP protocol for annihilation to Alanine will be used, as implemented in the original version of QresFEP reported in [Jespers et al., 2019](https://pubs.acs.org/doi/full/10.1021/acs.jctc.9b00538).
- In the case a cofactor (e.g. ligand) should be included, use the -c COFACTOR argument. 
- With the -f FORCEFIELD argument the forcefield of use is indicated. Note that one should be consistent wit the corresponding option indicated in step 1. (protein Preparation), whicle in this case the parameters of the ligand will be obtained with the corresponding routine associated to that forcefield.
- The -w WINDOWS argument is used to indicate the number of FEP windows per
FEP stage. Note here that if flag -d is indicated (i.e. dual topology protocol), the associated FEP transformation is done in 2 consecutive FEP stages, in which atomic charges are gradually annihilated/created, and Van der Waals parameters transition through a soft-core stage. This is translated in an effective doubling of the value of this flag (i.e. -w 25 means a total 2 x 25 = 50 FEP windows). On the other hand, using the single topology approach involves a number of FEP stages (5 - 9 stages, depending on the residue side-chain size).
- The -s SAMPLING option determines how the FEP windows will be distributed along the transformation, either in evenly distributed steps (-s linear) or different statistical strategies for unevenly distributed steps (sigmoidal,exponential,reverse_exponential). The differences are explained in detail in the Q manual 
- With the -l START argument one determines the lambda window where the simulations start: lambda state 1 (i.e. wt residue), 0 (mut residue) or 0.5 (the 50/50 mix state). The direction(s) of the tranformation is subsequently determined by the starting value: 0 --> 1; 1 --> 0; 0 <-- 0.5 --> 1. 
- The -ts TIMESTEP argument determines the time step, which is either 1 or 2 fs, in which case the value of the SHAKE algorithm dynamically adjusted accordingly.
- The temperature is set with the -T TEMPERATURE argument. 
- The number of simulation replicates is indicated with the -r REPLICATES argument
- Finally the HPC cluster where the simulations are run is optionally determined with the -C CLUSTER flag
- Optionally, the location of the `qprep` executable can be indicated with the -P PREPLOCATION flag (otherwise the `settings.py` will point to the location of the Q executables).

Now run the following command:

    QresFEP.py -m LEU39ALA -mc A -S protein -t A -d -f OPLSAAM -w 25 -s exponential -l 1 -ts 2fs -T 298 -r 10
    

The following files should have been generated:


| File                                              | Explanation |
| ------------------------------------------------- | ----------- |
| FEP_LEU39ALA/                                     | The FEP directory within it simulation input files and submit files|
| FEP_LEU39ALA/inputfiles/                          | Directory containing `qprep` and `qdyn` inputfiles| 
| FEP_LEU39ALA/inputfiles/complexnotexcluded.pdb    | PDB file of atoms within the simulation sphere only - necessary for trajectory visualisation with molecular visualisation software | 
| FEP_LEU39ALA/inputfiles/complex.pdb               | Input pdb file for `qprep` | 
| FEP_LEU39ALA/inputfiles/dualtop.top               | Topology file for `qdyn` generated by `qprep` | 
| FEP_LEU39ALA/inputfiles/eq*.inp                   | equilibration input files for `qdyn`| 
| FEP_LEU39ALA/inputfiles/FEP*fep                   | FEP input files for `qdyn` - describe the transforation strategy| 
| FEP_LEU39ALA/inputfiles/md_{1/2}_*.inp            | production input files for `qdyn`; the first digit {1/2} indicates the FEP stage, the second and third number indicate the FEP window) |
| FEP_LEU39ALA/inputfiles/OPLSAAM.lib               | forcefield library input file for `qprep` |
| FEP_LEU39ALA/inputfiles/OPLSAAM_merged.prm        | forcefield parameter input file for `qprep` |
| FEP_LEU39ALA/inputfiles/qfep.inp                  | input file for `qfep` |
| FEP_LEU39ALA/inputfiles/qprep.inp                 | input file for `qprep` - it was run during `QresFEP.py` to generate the topology file (*dualtop.top*) |
| FEP_LEU39ALA/inputfiles/qprep.out                 | Output file from `qprep` |
| FEP_LEU39ALA/inputfiles/L2A.lib                   | hybrid residue library file for `qprep` |
| FEP_LEU39ALA/inputfiles/runCLUSTER.sh             | SLURM runfile for simulation submission on the HPC CLUSTER |
| FEP_LEU39ALA/inputfiles/top_p.pdb                 | pdb file of the topology (*dualtop.top*) |
| FEP_LEU39ALA/inputfiles/water.pdb                 | input water pdb file for `qprep`|

If up until this part all files have been generated properly, we can be confident that our QligFEP repository was succesfully installed and initiated.

Now it is time to setup all simulations in an automated fashion.
For this we will use the `mutagenesis.sh` script, which requires that:

- All mutant PDB files have being generated (i.e., the PyMOL mutagenesis has been succesfully executed
- A text file indicating the mutations to perform is present, wiht each mutation in the form WTxxMUT (i.e. LEU39ALA, MET120ALA, etc). In our case, *mutations_neutral.txt* indicates the collection of 9 mutations in this tutorial. 
- Note that one can run the *mutations_charged.txt* set of mutations in principle, though those imply charge-changing mutations (like Asp to Ala, etc), which use a new routine to account for the charge-changing situation that will be soon reported.
If you made sure that this is indeed the case, then simply run the script. (You can change *mutations_neutral.txt* to *mutations_charged.txt* if you set up for the charged mutations instead).
Now we can run the script:

    bash mutagenesis.sh

This will take some time to generate all QresFEP inputfiles, but for each loop, you should see the follwing being printed out to your terminal

    LEU39ALA
    STOP qprep ended normally
    Generating dual topology inputfiles
    FEP_LEU39ALA
    STOP qprep ended normally
    STOP qprep ended normally
    Generating dual topology inputfiles
    FEP_LEU39ALA
    STOP qprep ended normally
    PyMOL executed successfully
    STOP qprep ended normally

The following files should have been generated:

| File                      | Explanation |
| ------------------------- | ----------- |
|  protein/FEP_LEU39ALA/    | protein directory for protein runs, within it the FEP directory of the LEU39ALA mutation|
| protein/FEP_* ...         | The protein FEP directories of all other mutations)|
| tripeptide/FEP_LEU39ALA/  | Tripeptide directory for reference tripeptide runs, within it the FEP directory of the LEU39ALA mutation|
| tripeptide/FEP_* ...      | The tripeptide FEP directories of all other mutations|

Now we can start submitting our FEP simulations. We typically run these simulations on an HPC cluster, so make sure that your `settings.py` script is updated with correct setting for your CLUSTER.
You can then submit the simulations. 

We will initially test that we can run the 10 replicate simulations corresponding to one leg of the thermodynamic cycle of one mutation. Place yourselt into one of the FEP directories and simply run the submit command therein:

    cd protein/FEP_LEU39ALA
    bash FEP_submit.sh

If this submits 10 simulations, we are basically done!

In principle, you should do this for all mutations times two (i.e., protein and tripeptide directories therein). For our convenience, we can simply use a batch script that does this iteration for us. We provide one example, but feel free to design such a script in any way that you seem fit for your specific case. For this tutorial we can use the example script `submit.sh`:

    bash submit.sh

You should have 10 (mutations) x 2 (legs) x 10 (replicas) = 200 jobs submitted. Depending on how many cores per job (e.g., 8 or 16) you have requested, as well as how long the equilibration phase you set, these jobs will run for several hours to complete.
Set this in `settings.py` to optimize for your cluster entry and have all the mutations done as quick as possible on your cluster. 

# 4. Analyse FEP simulation output files

When the simulations have succesfully finished, output files and directories should have been generated. Take a look at the provided example output in the directory FEP_example.
Besides copied files from the inputfiles/ directory, it contains the following new simulation output files:

| File                                      | Explanation |
| ----------------------------------------- | ----------- |
| FEP_example/slurm-001.out                 | SLURM output file, with output print of the simulation steps|
| FEP_example/FEP{1/2}/                     | output directories of each FEP stage (FEP1 and FEP2) |
| FEP_example/FEP{1/2}/298/                 | temperature directory - if ran at several temperatures, multiple directories will be listed here |
| FEP_example/FEP{1/2}/298/{1}              | replicate directory - if ran for a series of replicates, multiple directories will be listed here |
| FEP_example/FEP{1/2}/298/{1}/*.log        | `qdyn` output log files for each simulation step (equilibration and FEP windows) |
| FEP_example/FEP{1/2}/298/{1}/*.dcd        | `qdyn` output trajectory files - load into molecular visualisation software to observe |
| FEP_example/FEP{1/2}/298/{1}/*.re         | `qdyn` restart files - saved coordinates and velocities to resume simulation on the next window |
| FEP_example/FEP{1/2}/298/{1}/*.en         | `qdyn` output energy files - contains the energy evolution of the simulation |
| FEP_example/FEP{1/2}/298/{1}qfep.inp      | modifies `qfep` input file to read all generated energy files |
| FEP_example/FEP{1/2}/298/{1}/qfep.out     | `qfep` output file summarizing the energy of the FEP stage transformation |

We will now make use of the analysis tool `analyze_FEP.py` to facilitate processing the simulation output from the generated `qfep` output files.
Run the following command to see the help menu:

    analyze_FEP.py -h

This script requires the FEP directory as input argument to perform the analysis for with the -F FEP flag. 

- Which temperature should be analysed for is set with the -T TEMPERATURE flag.
- The simulation direction is indicated with the -l START flag to determine the correct sign of the free energy difference. 
- Optionally, we can choose to indicate the cluster of analysis with the -C CLUSTER flag. 
- to write out PDB files of the simulation trajectories with the -pdb PDB flag, and 
- wether to plot energies in either red of blue color with -c COLOR.
- Note: End-state-catastrophe is not relevant for multiple stage amino acid FEP. 

You can analyze your own simulation results, or use the `FEP_example` directory provided:

    analyze_FEP -F FEP_examlpe -T 298 -l 1

The following should be printed out to your terminal (scroll the horizontal bar for an explanatinon of each line entry):

    dG    FEP1  65.80  64.16  64.14  63.97  65.09  62.16  63.71  66.45  65.96  63.24     (Zwanzig dG for each replicate in FEP1 stage)
    dG    FEP1  65.80  64.16  64.14  63.97  65.09  62.16  63.71  66.45  65.96  63.24     (Zwanzig dG for each replicate in FEP2 stage)
    dGf   FEP1  65.87  64.12  64.19  64.01  65.13  62.10  63.83  66.45  66.00  63.16     (forward dG for each replicate in FEP1 stage)
    dGf   FEP2   6.06   3.24   4.36   3.36   3.83   5.83   6.80   4.95   5.07   6.92     (forward dG for each replicate in FEP2 stage)
    dGr   FEP1 -65.74 -64.19 -64.09 -63.93 -65.05 -62.23 -63.59 -66.44 -65.93 -63.33     (reverse dG for each replicate in FEP1 stage)
    dGr   FEP2  -5.74  -2.70  -4.19  -2.95  -3.14  -5.19  -6.46  -4.31  -4.60  -6.38     (reverse dG for each replicate in FEP2 stage)
    dGbar FEP1  65.81  64.17  64.13  63.98  65.10  62.15  63.72  66.44  65.98  63.23     (Bennett Acceptamce Ratop dG for each replicate in FEP1 stage)
    dGbar FEP2   5.85   2.94   4.31   3.15   3.46   5.49   6.61   4.62   4.79   6.63     (Bennett Acceptamce Ratop dG for each replicate in FEP2 stage)
    FEP_example  69.27  0.54  69.53  0.54 -69.02  0.53  69.26  0.53                      (Average and Standard Error for dG, dGf, dGr, and dGbar)
    crashes  0                                                                           (failed simulations)

We can also perform this step in an automated manner, either write your own script to do this in a way you seem fit, or use the provided `analyze.sh` script:

    bash analyze.sh

The following files should have been created

| File                                  | Explanation |
| ------------------------------------- | ----------- |
| protein/analyze_protein.txt           | Analyze_FEP summary of each FEP directory listed in the protein directory|
| tripeptide/analyze_tripeptide.txt     | Analyze_FEP summary of each FEP directory listed in the tripeptide directory|

For each mutation, the free energy difference between the protein and the reference tripeptide systems due to a mutation equals the relative folding free energy shift, expressed as ΔΔG. We can calculate this as follows:

ΔΔG_mutation = ΔG_FEP_protein - ΔG_FEP_tripeptide (in kcal/mol)