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
qdyn=/data1/projects/pi-gerard/Q/bin/qdynp
fepfiles=(FEP3.fep)
temperature=298
run=10
finalMDrestart=md_0000_1000.re

workdir=/home/jespers/adenosine/1.A1-A2A_selectivity/A1/5.FEP/holo
inputfiles=/home/jespers/adenosine/1.A1-A2A_selectivity/A1/5.FEP/holo/inputfiles
fepfile=FEP3.fep
fepdir=$workdir/FEP3
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

lastfep=FEP2
cp $workdir/$lastfep/$temperature/$run/$finalMDrestart $rundir/eq5.re

sed -i s/T_VAR/"$temperature"/ *.inp
sed -i s/FEP_VAR/"$fepfile"/ *.inp

#EQ_FILES

#RUN_FILES
time mpirun -np 24 $qdyn md_1000_0000.inp > md_1000_0000.log
time mpirun -np 24 $qdyn md_0999_0001.inp > md_0999_0001.log
time mpirun -np 24 $qdyn md_0997_0003.inp > md_0997_0003.log
time mpirun -np 24 $qdyn md_0996_0004.inp > md_0996_0004.log
time mpirun -np 24 $qdyn md_0994_0006.inp > md_0994_0006.log
time mpirun -np 24 $qdyn md_0993_0007.inp > md_0993_0007.log
time mpirun -np 24 $qdyn md_0991_0009.inp > md_0991_0009.log
time mpirun -np 24 $qdyn md_0989_0011.inp > md_0989_0011.log
time mpirun -np 24 $qdyn md_0987_0013.inp > md_0987_0013.log
time mpirun -np 24 $qdyn md_0985_0015.inp > md_0985_0015.log
time mpirun -np 24 $qdyn md_0982_0018.inp > md_0982_0018.log
time mpirun -np 24 $qdyn md_0980_0020.inp > md_0980_0020.log
time mpirun -np 24 $qdyn md_0977_0023.inp > md_0977_0023.log
time mpirun -np 24 $qdyn md_0975_0025.inp > md_0975_0025.log
time mpirun -np 24 $qdyn md_0971_0029.inp > md_0971_0029.log
time mpirun -np 24 $qdyn md_0968_0032.inp > md_0968_0032.log
time mpirun -np 24 $qdyn md_0964_0036.inp > md_0964_0036.log
time mpirun -np 24 $qdyn md_0960_0040.inp > md_0960_0040.log
time mpirun -np 24 $qdyn md_0956_0044.inp > md_0956_0044.log
time mpirun -np 24 $qdyn md_0951_0049.inp > md_0951_0049.log
time mpirun -np 24 $qdyn md_0946_0054.inp > md_0946_0054.log
time mpirun -np 24 $qdyn md_0940_0060.inp > md_0940_0060.log
time mpirun -np 24 $qdyn md_0933_0067.inp > md_0933_0067.log
time mpirun -np 24 $qdyn md_0926_0074.inp > md_0926_0074.log
time mpirun -np 24 $qdyn md_0917_0083.inp > md_0917_0083.log
time mpirun -np 24 $qdyn md_0907_0093.inp > md_0907_0093.log
time mpirun -np 24 $qdyn md_0896_0104.inp > md_0896_0104.log
time mpirun -np 24 $qdyn md_0883_0117.inp > md_0883_0117.log
time mpirun -np 24 $qdyn md_0867_0133.inp > md_0867_0133.log
time mpirun -np 24 $qdyn md_0847_0153.inp > md_0847_0153.log
time mpirun -np 24 $qdyn md_0824_0176.inp > md_0824_0176.log
time mpirun -np 24 $qdyn md_0793_0207.inp > md_0793_0207.log
time mpirun -np 24 $qdyn md_0754_0246.inp > md_0754_0246.log
time mpirun -np 24 $qdyn md_0700_0300.inp > md_0700_0300.log
time mpirun -np 24 $qdyn md_0622_0378.inp > md_0622_0378.log
time mpirun -np 24 $qdyn md_0500_0500.inp > md_0500_0500.log
time mpirun -np 24 $qdyn md_0378_0622.inp > md_0378_0622.log
time mpirun -np 24 $qdyn md_0300_0700.inp > md_0300_0700.log
time mpirun -np 24 $qdyn md_0246_0754.inp > md_0246_0754.log
time mpirun -np 24 $qdyn md_0207_0793.inp > md_0207_0793.log
time mpirun -np 24 $qdyn md_0176_0824.inp > md_0176_0824.log
time mpirun -np 24 $qdyn md_0153_0847.inp > md_0153_0847.log
time mpirun -np 24 $qdyn md_0133_0867.inp > md_0133_0867.log
time mpirun -np 24 $qdyn md_0117_0883.inp > md_0117_0883.log
time mpirun -np 24 $qdyn md_0104_0896.inp > md_0104_0896.log
time mpirun -np 24 $qdyn md_0093_0907.inp > md_0093_0907.log
time mpirun -np 24 $qdyn md_0083_0917.inp > md_0083_0917.log
time mpirun -np 24 $qdyn md_0074_0926.inp > md_0074_0926.log
time mpirun -np 24 $qdyn md_0067_0933.inp > md_0067_0933.log
time mpirun -np 24 $qdyn md_0060_0940.inp > md_0060_0940.log
time mpirun -np 24 $qdyn md_0054_0946.inp > md_0054_0946.log
time mpirun -np 24 $qdyn md_0049_0951.inp > md_0049_0951.log
time mpirun -np 24 $qdyn md_0044_0956.inp > md_0044_0956.log
time mpirun -np 24 $qdyn md_0040_0960.inp > md_0040_0960.log
time mpirun -np 24 $qdyn md_0036_0964.inp > md_0036_0964.log
time mpirun -np 24 $qdyn md_0032_0968.inp > md_0032_0968.log
time mpirun -np 24 $qdyn md_0029_0971.inp > md_0029_0971.log
time mpirun -np 24 $qdyn md_0025_0975.inp > md_0025_0975.log
time mpirun -np 24 $qdyn md_0023_0977.inp > md_0023_0977.log
time mpirun -np 24 $qdyn md_0020_0980.inp > md_0020_0980.log
time mpirun -np 24 $qdyn md_0018_0982.inp > md_0018_0982.log
time mpirun -np 24 $qdyn md_0015_0985.inp > md_0015_0985.log
time mpirun -np 24 $qdyn md_0013_0987.inp > md_0013_0987.log
time mpirun -np 24 $qdyn md_0011_0989.inp > md_0011_0989.log
time mpirun -np 24 $qdyn md_0009_0991.inp > md_0009_0991.log
time mpirun -np 24 $qdyn md_0007_0993.inp > md_0007_0993.log
time mpirun -np 24 $qdyn md_0006_0994.inp > md_0006_0994.log
time mpirun -np 24 $qdyn md_0004_0996.inp > md_0004_0996.log
time mpirun -np 24 $qdyn md_0003_0997.inp > md_0003_0997.log
time mpirun -np 24 $qdyn md_0001_0999.inp > md_0001_0999.log
time mpirun -np 24 $qdyn md_0000_1000.inp > md_0000_1000.log
timeout 30s /data1/projects/pi-gerard/Q/bin/qfep < qfep.inp > qfep.out
done
