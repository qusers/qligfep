#!/bin/bash
#
#SBATCH --nodes=NODES
#SBATCH --ntasks-per-node=NTASKS
#SBATCH --job-name=JOBNAME 
#SBATCH -A ACCOUNT
#SBATCH -p PARTITION
#              d-hh:mm:ss
#SBATCH --time=TIME

## Load modules for qdynp
MODULES

## define qdynp location
QDYN
fepfiles=(FEPS)
temperature=298
run=10
finalMDrestart=RESTART

seed=1

workdir=/WORKDIR
inputfiles=/INPUTFILES

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

cp $inputfiles/md_$((index+1))_*.inp .
cp $inputfiles/*.top .
cp $inputfiles/qfep.inp .
cp $inputfiles/$fepfile .

if [ $index -lt 1 ]; then
cp $inputfiles/eq*.inp .
#sed -i s/SEED_VAR/"$[1 + $[RANDOM % 9999]]"/ eq1.inp
sed -i s/SEED_VAR/"$seed"/ eq1.inp
else
lastfep=FEP$index
cp $workdir/$lastfep/$temperature/$run/$finalMDrestart $rundir/eq5.re
fi

sed -i s/T_VAR/"$temperature"/ *.inp
sed -i s/FEP_VAR/"$fepfile"/ *.inp
echo $fepfile
if [ $index -lt 1 ]; then
echo $run
#EQ_FILES
#RUN_1_FILES
else
#RUN_2_FILES
fi
printf "%s\n" $(ls -1 md*en) | tac >> qfep.inp
QFEP < qfep.inp > qfep.out
#rm *.log
#rm *.dcd
#mv md_1_0000_1000.re md_1_0000_1000.re.keep
#mv md_2_0000_1000.re md_2_0000_1000.re.keep
#rm *.re
#mv md_1_0000_1000.re.keep md_1_0000_1000.re
#mv md_2_0000_1000.re.keep md_2_0000_1000.re
#rm *.inp
done
