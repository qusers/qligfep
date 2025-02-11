#!/bin/bash
#
#SBATCH -c NTASKS
#SBATCH --job-name=JOBNAME 
#SBATCH -A ACCOUNT
#SBATCH -p PARTITION
#              d-hh:mm:ss
#SBATCH --time=TIME

## Load modules for qdynp
MODULES

## define qdynp location
QDYN
temperature=298
run=10

seed=1

workdir=/WORKDIR
inputfiles=/INPUTFILES



rundir=runs/$run
mkdir -p $rundir
cd $rundir

cp $inputfiles/md*.inp .
cp $inputfiles/*.top .
cp $inputfiles/eq*.inp .
#sed -i s/SEED_VAR/"$[1 + $[RANDOM % 9999]]"/ eq1.inp
sed -i s/SEED_VAR/"$seed"/ eq1.inp
sed -i s/T_VAR/"$temperature"/ *.inp
sed -i s/FEP_VAR/"$fepfile"/ *.inp

echo $run
#time srun $qdyn eq1.inp > eq1.log
#EQ_FILES

#RUN_FILES
