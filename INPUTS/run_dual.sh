#!/bin/bash
#
#SBATCH --nodes=NODES
#SBATCH --ntasks-per-node=NTASKS
#SBATCH -A ACCOUNT 
#SBATCH -p CLUSTER-AMD
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

workdir=/WORKDIR
inputfiles=/INPUTFILES
check=/STOPLIE
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

cp $inputfiles/*.top .
cp $inputfiles/qfep.inp .
cp $inputfiles/$fepfile .
cp $inputfiles/eq*.inp .

sed -i s/SEED_VAR/"$[1 + $[RANDOM % 9999]]"/ eq1.inp
if [ $index -lt 1 ]; then
sed -i s/'0.000 1.000'/'1.000 0.000'/ eq*inp
cp $inputfiles/md*_F.inp .
else
sed -i s/'1.000 0.000'/'0.000 1.000'/ eq*inp
cp $inputfiles/md*_R.inp .
fi  

sed -i s/T_VAR/"$temperature"/ *.inp
sed -i s/FEP_VAR/"$fepfile"/ *.inp
if [ $index -lt 1 ]; then
echo $run
echo $fepfile
#EQ_FILES
#RUN_FILES
timeout 30s QFEP < qfep.inp > qfep.out
else
echo $fepfile
#EQ_FILES
#RUN_FILES_R
timeout 30s QFEP < qfep.inp > qfep.out
fi
done