#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#              d-hh:mm:ss
#SBATCH --time=0-12:00:00

export TMP=/tmp
export TEMP=/tmp
export TMPDIR=/tmp
## Load modules for qdynp
module load 2023 iimpi/2023a

softcore_leg() {
local index=$1
sed -i s/T_VAR/"$temperature"/ *.inp
sed -i s/FEP_VAR/"$fepfile"/ *.inp

#EQ_FILES
time mpirun -np 16 $qdyn eq1.inp > eq1.log
time mpirun -np 16 $qdyn eq2.inp > eq2.log
time mpirun -np 16 $qdyn eq3.inp > eq3.log
time mpirun -np 16 $qdyn eq4.inp > eq4.log
time mpirun -np 16 $qdyn eq5.inp > eq5.log

#RUN_FILES
time mpirun -np 16 $qdyn md_1000_0000.inp > md_1000_0000.log
time mpirun -np 16 $qdyn md_0980_0020.inp > md_0980_0020.log
time mpirun -np 16 $qdyn md_0960_0040.inp > md_0960_0040.log
time mpirun -np 16 $qdyn md_0940_0060.inp > md_0940_0060.log
time mpirun -np 16 $qdyn md_0920_0080.inp > md_0920_0080.log
time mpirun -np 16 $qdyn md_0900_0100.inp > md_0900_0100.log
time mpirun -np 16 $qdyn md_0880_0120.inp > md_0880_0120.log
time mpirun -np 16 $qdyn md_0860_0140.inp > md_0860_0140.log
time mpirun -np 16 $qdyn md_0840_0160.inp > md_0840_0160.log
time mpirun -np 16 $qdyn md_0820_0180.inp > md_0820_0180.log
time mpirun -np 16 $qdyn md_0800_0200.inp > md_0800_0200.log
time mpirun -np 16 $qdyn md_0780_0220.inp > md_0780_0220.log
time mpirun -np 16 $qdyn md_0760_0240.inp > md_0760_0240.log
time mpirun -np 16 $qdyn md_0740_0260.inp > md_0740_0260.log
time mpirun -np 16 $qdyn md_0720_0280.inp > md_0720_0280.log
time mpirun -np 16 $qdyn md_0700_0300.inp > md_0700_0300.log
time mpirun -np 16 $qdyn md_0680_0320.inp > md_0680_0320.log
time mpirun -np 16 $qdyn md_0660_0340.inp > md_0660_0340.log
time mpirun -np 16 $qdyn md_0640_0360.inp > md_0640_0360.log
time mpirun -np 16 $qdyn md_0620_0380.inp > md_0620_0380.log
time mpirun -np 16 $qdyn md_0600_0400.inp > md_0600_0400.log
time mpirun -np 16 $qdyn md_0580_0420.inp > md_0580_0420.log
time mpirun -np 16 $qdyn md_0560_0440.inp > md_0560_0440.log
time mpirun -np 16 $qdyn md_0540_0460.inp > md_0540_0460.log
time mpirun -np 16 $qdyn md_0520_0480.inp > md_0520_0480.log
time mpirun -np 16 $qdyn md_0500_0500.inp > md_0500_0500.log
time mpirun -np 16 $qdyn md_0480_0520.inp > md_0480_0520.log
time mpirun -np 16 $qdyn md_0460_0540.inp > md_0460_0540.log
time mpirun -np 16 $qdyn md_0440_0560.inp > md_0440_0560.log
time mpirun -np 16 $qdyn md_0420_0580.inp > md_0420_0580.log
time mpirun -np 16 $qdyn md_0400_0600.inp > md_0400_0600.log
time mpirun -np 16 $qdyn md_0380_0620.inp > md_0380_0620.log
time mpirun -np 16 $qdyn md_0360_0640.inp > md_0360_0640.log
time mpirun -np 16 $qdyn md_0340_0660.inp > md_0340_0660.log
time mpirun -np 16 $qdyn md_0320_0680.inp > md_0320_0680.log
time mpirun -np 16 $qdyn md_0300_0700.inp > md_0300_0700.log
time mpirun -np 16 $qdyn md_0280_0720.inp > md_0280_0720.log
time mpirun -np 16 $qdyn md_0260_0740.inp > md_0260_0740.log
time mpirun -np 16 $qdyn md_0240_0760.inp > md_0240_0760.log
time mpirun -np 16 $qdyn md_0220_0780.inp > md_0220_0780.log
time mpirun -np 16 $qdyn md_0200_0800.inp > md_0200_0800.log
time mpirun -np 16 $qdyn md_0180_0820.inp > md_0180_0820.log
time mpirun -np 16 $qdyn md_0160_0840.inp > md_0160_0840.log
time mpirun -np 16 $qdyn md_0140_0860.inp > md_0140_0860.log
time mpirun -np 16 $qdyn md_0120_0880.inp > md_0120_0880.log
time mpirun -np 16 $qdyn md_0100_0900.inp > md_0100_0900.log
time mpirun -np 16 $qdyn md_0080_0920.inp > md_0080_0920.log
time mpirun -np 16 $qdyn md_0060_0940.inp > md_0060_0940.log
time mpirun -np 16 $qdyn md_0040_0960.inp > md_0040_0960.log
time mpirun -np 16 $qdyn md_0020_0980.inp > md_0020_0980.log
time mpirun -np 16 $qdyn md_0000_1000.inp > md_0000_1000.log
timeout 30s /home/wjespers/software/Q/bin/qfep < qfep2.inp > qfep2.out
}

charge_leg() {
local index=$1
lastfep=FEP$((index-1))
cp $workdir/$lastfep/$temperature/$run/$finalMDrestart $rundir/eq5.re
sed -i s/T_VAR/"$temperature"/ *.inp
sed -i s/FEP_VAR/"$fepfile"/ *.inp
#RUN_FILES
time mpirun -np 16 $qdyn md_1000_0000.inp > md_1000_0000.log
time mpirun -np 16 $qdyn md_0980_0020.inp > md_0980_0020.log
time mpirun -np 16 $qdyn md_0960_0040.inp > md_0960_0040.log
time mpirun -np 16 $qdyn md_0940_0060.inp > md_0940_0060.log
time mpirun -np 16 $qdyn md_0920_0080.inp > md_0920_0080.log
time mpirun -np 16 $qdyn md_0900_0100.inp > md_0900_0100.log
time mpirun -np 16 $qdyn md_0880_0120.inp > md_0880_0120.log
time mpirun -np 16 $qdyn md_0860_0140.inp > md_0860_0140.log
time mpirun -np 16 $qdyn md_0840_0160.inp > md_0840_0160.log
time mpirun -np 16 $qdyn md_0820_0180.inp > md_0820_0180.log
time mpirun -np 16 $qdyn md_0800_0200.inp > md_0800_0200.log
time mpirun -np 16 $qdyn md_0780_0220.inp > md_0780_0220.log
time mpirun -np 16 $qdyn md_0760_0240.inp > md_0760_0240.log
time mpirun -np 16 $qdyn md_0740_0260.inp > md_0740_0260.log
time mpirun -np 16 $qdyn md_0720_0280.inp > md_0720_0280.log
time mpirun -np 16 $qdyn md_0700_0300.inp > md_0700_0300.log
time mpirun -np 16 $qdyn md_0680_0320.inp > md_0680_0320.log
time mpirun -np 16 $qdyn md_0660_0340.inp > md_0660_0340.log
time mpirun -np 16 $qdyn md_0640_0360.inp > md_0640_0360.log
time mpirun -np 16 $qdyn md_0620_0380.inp > md_0620_0380.log
time mpirun -np 16 $qdyn md_0600_0400.inp > md_0600_0400.log
time mpirun -np 16 $qdyn md_0580_0420.inp > md_0580_0420.log
time mpirun -np 16 $qdyn md_0560_0440.inp > md_0560_0440.log
time mpirun -np 16 $qdyn md_0540_0460.inp > md_0540_0460.log
time mpirun -np 16 $qdyn md_0520_0480.inp > md_0520_0480.log
time mpirun -np 16 $qdyn md_0500_0500.inp > md_0500_0500.log
time mpirun -np 16 $qdyn md_0480_0520.inp > md_0480_0520.log
time mpirun -np 16 $qdyn md_0460_0540.inp > md_0460_0540.log
time mpirun -np 16 $qdyn md_0440_0560.inp > md_0440_0560.log
time mpirun -np 16 $qdyn md_0420_0580.inp > md_0420_0580.log
time mpirun -np 16 $qdyn md_0400_0600.inp > md_0400_0600.log
time mpirun -np 16 $qdyn md_0380_0620.inp > md_0380_0620.log
time mpirun -np 16 $qdyn md_0360_0640.inp > md_0360_0640.log
time mpirun -np 16 $qdyn md_0340_0660.inp > md_0340_0660.log
time mpirun -np 16 $qdyn md_0320_0680.inp > md_0320_0680.log
time mpirun -np 16 $qdyn md_0300_0700.inp > md_0300_0700.log
time mpirun -np 16 $qdyn md_0280_0720.inp > md_0280_0720.log
time mpirun -np 16 $qdyn md_0260_0740.inp > md_0260_0740.log
time mpirun -np 16 $qdyn md_0240_0760.inp > md_0240_0760.log
time mpirun -np 16 $qdyn md_0220_0780.inp > md_0220_0780.log
time mpirun -np 16 $qdyn md_0200_0800.inp > md_0200_0800.log
time mpirun -np 16 $qdyn md_0180_0820.inp > md_0180_0820.log
time mpirun -np 16 $qdyn md_0160_0840.inp > md_0160_0840.log
time mpirun -np 16 $qdyn md_0140_0860.inp > md_0140_0860.log
time mpirun -np 16 $qdyn md_0120_0880.inp > md_0120_0880.log
time mpirun -np 16 $qdyn md_0100_0900.inp > md_0100_0900.log
time mpirun -np 16 $qdyn md_0080_0920.inp > md_0080_0920.log
time mpirun -np 16 $qdyn md_0060_0940.inp > md_0060_0940.log
time mpirun -np 16 $qdyn md_0040_0960.inp > md_0040_0960.log
time mpirun -np 16 $qdyn md_0020_0980.inp > md_0020_0980.log
time mpirun -np 16 $qdyn md_0000_1000.inp > md_0000_1000.log
timeout 30s /home/wjespers/software/Q/bin/qfep < qfep3.inp > qfep3.out
}

## define qdynp location
qdyn=/home/wjespers/software/Q/bin/qdynp
fepfiles=("FEP1.fep" "FEP2.fep" "FEP3.fep")
temperature=298
run=10
finalMDrestart=md_0000_1000.re

workdir=/gpfs/scratch1/nodespecific/int6/yricky/softcore_with_long_endpoint_sampling/butylbenzene_flag_0_seq5_soft5_15_radius_eq_15
inputfiles=/gpfs/scratch1/nodespecific/int6/yricky/softcore_with_long_endpoint_sampling/butylbenzene_flag_0_seq5_soft5_15_radius_eq_15/inputfiles

for ((index=2; index<=${#fepfiles[@]}; index++)); do
# do this for every FEP file, use index in the list.

fepfile=FEP$index.fep
fepstep=FEP$index
fepdir=$workdir/FEP$index
mkdir -p $fepdir
cd $fepdir
tempdir=$fepdir/$temperature
mkdir -p $tempdir
cd $tempdir
mdfilesdir=$inputfiles/mdfiles$index
rundir=$tempdir/$run
mkdir -p $rundir
cd $rundir

cp $mdfilesdir/md*.inp .
cp $inputfiles/*.top .
cp $inputfiles/qfep$index.inp .
cp $inputfiles/$fepfile .

# Run the appropriate job function based on the filename
if [ "$fepstep" == "FEP1" ]; then
vdW_leg $index
elif [ "$fepstep" == "FEP2" ]; then
softcore_leg $index
elif [ "$fepstep" == "FEP3" ]; then
charge_leg $index
fi
done