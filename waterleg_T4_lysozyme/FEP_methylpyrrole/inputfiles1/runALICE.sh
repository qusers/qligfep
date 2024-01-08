#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#              d-hh:mm:ss
#SBATCH --time=0-3:00:00

export TMP=/tmp
export TEMP=/tmp
export TMPDIR=/tmp
## Load modules for qdynp
module load OpenMPI/3.1.3-GCC-8.2.0-2.31.1

## define qdynp location
qdyn=/home/jespersw/software/q6/bin/qdynp
fepfiles=(FEP1.fep)
temperature=298
run=10
finalMDrestart=md_0000_1000.re

workdir=/home/jespers/adenosine/1.A1-A2A_selectivity/A1/5.FEP/holo
inputfiles=/home/jespers/adenosine/1.A1-A2A_selectivity/A1/5.FEP/holo/inputfiles
fepfile=FEP1.fep
fepdir=$workdir/FEP1
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

cp $inputfiles/eq*.inp .
sed -i s/SEED_VAR/"$[1 + $[RANDOM % 9999]]"/ eq1.inp

sed -i s/T_VAR/"$temperature"/ *.inp
sed -i s/FEP_VAR/"$fepfile"/ *.inp

#EQ_FILES
time mpirun -np 24 $qdyn eq1.inp > eq1.log
time mpirun -np 24 $qdyn eq2.inp > eq2.log
time mpirun -np 24 $qdyn eq3.inp > eq3.log
time mpirun -np 24 $qdyn eq4.inp > eq4.log
time mpirun -np 24 $qdyn eq5.inp > eq5.log

#RUN_FILES
time mpirun -np 24 $qdyn md_1000_0000.inp > md_1000_0000.log
time mpirun -np 24 $qdyn md_0980_0020.inp > md_0980_0020.log
time mpirun -np 24 $qdyn md_0960_0040.inp > md_0960_0040.log
time mpirun -np 24 $qdyn md_0940_0060.inp > md_0940_0060.log
time mpirun -np 24 $qdyn md_0920_0080.inp > md_0920_0080.log
time mpirun -np 24 $qdyn md_0900_0100.inp > md_0900_0100.log
time mpirun -np 24 $qdyn md_0880_0120.inp > md_0880_0120.log
time mpirun -np 24 $qdyn md_0860_0140.inp > md_0860_0140.log
time mpirun -np 24 $qdyn md_0840_0160.inp > md_0840_0160.log
time mpirun -np 24 $qdyn md_0820_0180.inp > md_0820_0180.log
time mpirun -np 24 $qdyn md_0800_0200.inp > md_0800_0200.log
time mpirun -np 24 $qdyn md_0780_0220.inp > md_0780_0220.log
time mpirun -np 24 $qdyn md_0760_0240.inp > md_0760_0240.log
time mpirun -np 24 $qdyn md_0740_0260.inp > md_0740_0260.log
time mpirun -np 24 $qdyn md_0720_0280.inp > md_0720_0280.log
time mpirun -np 24 $qdyn md_0700_0300.inp > md_0700_0300.log
time mpirun -np 24 $qdyn md_0680_0320.inp > md_0680_0320.log
time mpirun -np 24 $qdyn md_0660_0340.inp > md_0660_0340.log
time mpirun -np 24 $qdyn md_0640_0360.inp > md_0640_0360.log
time mpirun -np 24 $qdyn md_0620_0380.inp > md_0620_0380.log
time mpirun -np 24 $qdyn md_0600_0400.inp > md_0600_0400.log
time mpirun -np 24 $qdyn md_0580_0420.inp > md_0580_0420.log
time mpirun -np 24 $qdyn md_0560_0440.inp > md_0560_0440.log
time mpirun -np 24 $qdyn md_0540_0460.inp > md_0540_0460.log
time mpirun -np 24 $qdyn md_0520_0480.inp > md_0520_0480.log
time mpirun -np 24 $qdyn md_0500_0500.inp > md_0500_0500.log
time mpirun -np 24 $qdyn md_0480_0520.inp > md_0480_0520.log
time mpirun -np 24 $qdyn md_0460_0540.inp > md_0460_0540.log
time mpirun -np 24 $qdyn md_0440_0560.inp > md_0440_0560.log
time mpirun -np 24 $qdyn md_0420_0580.inp > md_0420_0580.log
time mpirun -np 24 $qdyn md_0400_0600.inp > md_0400_0600.log
time mpirun -np 24 $qdyn md_0380_0620.inp > md_0380_0620.log
time mpirun -np 24 $qdyn md_0360_0640.inp > md_0360_0640.log
time mpirun -np 24 $qdyn md_0340_0660.inp > md_0340_0660.log
time mpirun -np 24 $qdyn md_0320_0680.inp > md_0320_0680.log
time mpirun -np 24 $qdyn md_0300_0700.inp > md_0300_0700.log
time mpirun -np 24 $qdyn md_0280_0720.inp > md_0280_0720.log
time mpirun -np 24 $qdyn md_0260_0740.inp > md_0260_0740.log
time mpirun -np 24 $qdyn md_0240_0760.inp > md_0240_0760.log
time mpirun -np 24 $qdyn md_0220_0780.inp > md_0220_0780.log
time mpirun -np 24 $qdyn md_0200_0800.inp > md_0200_0800.log
time mpirun -np 24 $qdyn md_0180_0820.inp > md_0180_0820.log
time mpirun -np 24 $qdyn md_0160_0840.inp > md_0160_0840.log
time mpirun -np 24 $qdyn md_0140_0860.inp > md_0140_0860.log
time mpirun -np 24 $qdyn md_0120_0880.inp > md_0120_0880.log
time mpirun -np 24 $qdyn md_0100_0900.inp > md_0100_0900.log
time mpirun -np 24 $qdyn md_0080_0920.inp > md_0080_0920.log
time mpirun -np 24 $qdyn md_0060_0940.inp > md_0060_0940.log
time mpirun -np 24 $qdyn md_0040_0960.inp > md_0040_0960.log
time mpirun -np 24 $qdyn md_0020_0980.inp > md_0020_0980.log
time mpirun -np 24 $qdyn md_0000_1000.inp > md_0000_1000.log
timeout 30s /home/jespersw/software/q6/bin/qfep < qfep.inp > qfep.out
done
