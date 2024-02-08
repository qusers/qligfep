#!/bin/bash
#
#SBATCH --nodes=NODES
#SBATCH --ntasks-per-node=NTASKS
#SBATCH --job-name=JOBNAME 
#SBATCH -A ACCOUNT
#SBATCH -p PARTITION
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
wt_res="RES_WT"
mut_res="RES_MUT"

seed=1

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

cp $inputfiles/md*.inp .
cp $inputfiles/*.top .
cp $inputfiles/qfep.inp .
cp $inputfiles/$fepfile .

if [ $index -lt 1 ]; then
cp $inputfiles/eq*.inp .
#sed -i s/SEED_VAR/"$[1 + $[RANDOM % 9999]]"/ eq1.inp
sed -i s/SEED_VAR/"$seed"/ eq1.inp
sed -i s/MUT_RES/"$mut_res"/ md_1000_0000.inp
sed -i s/MUT_RES/"$mut_res"/ md_0998_0002.inp
sed -i s/MUT_RES/"$mut_res"/ md_0996_0004.inp
sed -i s/MUT_RES/"$mut_res"/ md_0994_0006.inp
sed -i s/MUT_RES/"$mut_res"/ md_0991_0009.inp
sed -i s/MUT_RES/""/ md_*.inp
sed -i s/WT_RES/""/ md_*.inp
else
lastfep=FEP$index
cp $workdir/$lastfep/$temperature/$run/$finalMDrestart $rundir/eq5.re
sed -i s/WT_RES/"$wt_res"/ md_0000_1000.inp
sed -i s/WT_RES/"$wt_res"/ md_0002_0998.inp
sed -i s/WT_RES/"$wt_res"/ md_0004_0996.inp
sed -i s/WT_RES/"$wt_res"/ md_0006_0994.inp
sed -i s/WT_RES/"$wt_res"/ md_0009_0991.inp
sed -i s/MUT_RES/""/ md_*.inp
sed -i s/WT_RES/""/ md_*.inp
fi

sed -i s/T_VAR/"$temperature"/ *.inp
sed -i s/FEP_VAR/"$fepfile"/ *.inp
if [ $index -lt 1 ]; then
echo $run
#time srun $qdyn eq1.inp > eq1.log
#EQ_FILES
fi
echo $fepfile
#RUN_FILES
timeout 30s QFEP < qfep.inp > qfep.out
#rm *.log
#rm *.dcd
#mv md_0000_1000.re md_0000_1000.re.keep
#rm *.re
#mv md_0000_1000.re.keep md_0000_1000.re
done
