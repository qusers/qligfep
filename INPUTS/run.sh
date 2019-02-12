#!/bin/bash
#
#SBATCH --nodes=NODES
#SBATCH --ntasks-per-node=NTASKS
#SBATCH -A ACCOUNT 
#              d-hh:mm:ss
#SBATCH --time=TIME

export TMP=/tmp
export TEMP=/tmp
export TMPDIR=/tmp
## Load modules for qdynp
MODULES

## define qdynp location
QDYN
fepfiles=(FEPS)
temperature=298
run=10
finalMDrestart=md_0000_1000.re

workdir=/home/jespers/adenosine/1.A1-A2A_selectivity/A1/5.FEP/holo
inputfiles=/home/jespers/adenosine/1.A1-A2A_selectivity/A1/5.FEP/holo/inputfiles
length=${#fepfiles[@]}
length=$((length-1))
for index in $(seq 0 $length);do
fepfile=${fepfiles[$index]}
fepdir=$workdir/FEP$((index+1))
mkdir -p $fepdir
cd $fepdir
tempdir=$fepdir/$temperature
mkdir -p $tempdir
cd $tempdir

rundir=$tempdir/$run
mkdir -p $rundir
cd $rundir

cp $inputfiles/md*.inp .
cp $inputfiles/*.top .
cp $inputfiles/qfep.inp .
cp $inputfiles/$fepfile .
cp $inputfiles/run_0500-1000.sh .
cp $inputfiles/run_0500-0000.sh .

if [ $index -lt 1 ]; then
cp $inputfiles/eq*.inp .
sed -i s/SEED_VAR/"$[1 + $[RANDOM % 9999]]"/ eq1.inp
else
lastfep=FEP$index
cp $workdir/$lastfep/$temperature/$run/$finalMDrestart $rundir/eq5.re
fi

sed -i s/T_VAR/"$temperature"/ *.inp
sed -i s/FEP_VAR/"$fepfile"/ *.inp
if [ $index -lt 1 ]; then
#time mpirun -np 16 $qdyn eq1.inp > eq1.log
#EQ_FILES
fi
#RUN_FILES
timeout 30s QFEP < qfep.inp > qfep.out
done
