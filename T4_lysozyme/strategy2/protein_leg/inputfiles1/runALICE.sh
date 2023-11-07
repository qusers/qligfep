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
fepfiles=(FEP1.fep)
temperature=298
run=10
finalMDrestart=md_0000_1000.re

workdir=/data1/s2904160/software_qligfep/both_endpoints_sampling/benzene_endpoint_sampling_no_softcore/protein_leg
inputfiles=/data1/s2904160/software_qligfep/both_endpoints_sampling/benzene_endpoint_sampling_no_softcore/protein_leg/inputfiles1
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

cp $inputfiles/eq*.inp .
sed -i s/SEED_VAR/"$[1 + $[RANDOM % 9999]]"/ eq1.inp

sed -i s/T_VAR/"$temperature"/ *.inp
sed -i s/FEP_VAR/"$fepfile"/ *.inp

#time mpirun -np 16 $qdyn eq1.inp > eq1.log
#EQ_FILES
time mpirun -np 24 $qdyn eq1.inp > eq1.log
time mpirun -np 24 $qdyn eq2.inp > eq2.log
time mpirun -np 24 $qdyn eq3.inp > eq3.log
time mpirun -np 24 $qdyn eq4.inp > eq4.log
time mpirun -np 24 $qdyn eq5.inp > eq5.log

#RUN_FILES
time mpirun -np 24 $qdyn md_0500_0500.inp > md_0500_0500.log

time mpirun -np 12 $qdyn md_0514_0485.inp > md_0514_0485.log &
time mpirun -np 12 $qdyn md_0485_0514.inp > md_0485_0514.log

time mpirun -np 12 $qdyn md_0529_0470.inp > md_0529_0470.log &
time mpirun -np 12 $qdyn md_0470_0529.inp > md_0470_0529.log

time mpirun -np 12 $qdyn md_0544_0456.inp > md_0544_0456.log &
time mpirun -np 12 $qdyn md_0456_0544.inp > md_0456_0544.log

time mpirun -np 12 $qdyn md_0559_0441.inp > md_0559_0441.log &
time mpirun -np 12 $qdyn md_0441_0559.inp > md_0441_0559.log

time mpirun -np 12 $qdyn md_0573_0426.inp > md_0573_0426.log &
time mpirun -np 12 $qdyn md_0426_0573.inp > md_0426_0573.log

time mpirun -np 12 $qdyn md_0588_0412.inp > md_0588_0412.log &
time mpirun -np 12 $qdyn md_0412_0588.inp > md_0412_0588.log

time mpirun -np 12 $qdyn md_0603_0397.inp > md_0603_0397.log &
time mpirun -np 12 $qdyn md_0397_0603.inp > md_0397_0603.log

time mpirun -np 12 $qdyn md_0617_0382.inp > md_0617_0382.log &
time mpirun -np 12 $qdyn md_0382_0617.inp > md_0382_0617.log

time mpirun -np 12 $qdyn md_0632_0367.inp > md_0632_0367.log &
time mpirun -np 12 $qdyn md_0367_0632.inp > md_0367_0632.log

time mpirun -np 12 $qdyn md_0647_0353.inp > md_0647_0353.log &
time mpirun -np 12 $qdyn md_0353_0647.inp > md_0353_0647.log

time mpirun -np 12 $qdyn md_0662_0338.inp > md_0662_0338.log &
time mpirun -np 12 $qdyn md_0338_0662.inp > md_0338_0662.log

time mpirun -np 12 $qdyn md_0676_0323.inp > md_0676_0323.log &
time mpirun -np 12 $qdyn md_0323_0676.inp > md_0323_0676.log

time mpirun -np 12 $qdyn md_0691_0309.inp > md_0691_0309.log &
time mpirun -np 12 $qdyn md_0309_0691.inp > md_0309_0691.log

time mpirun -np 12 $qdyn md_0706_0294.inp > md_0706_0294.log &
time mpirun -np 12 $qdyn md_0294_0706.inp > md_0294_0706.log

time mpirun -np 12 $qdyn md_0720_0279.inp > md_0720_0279.log &
time mpirun -np 12 $qdyn md_0279_0720.inp > md_0279_0720.log

time mpirun -np 12 $qdyn md_0735_0265.inp > md_0735_0265.log &
time mpirun -np 12 $qdyn md_0265_0735.inp > md_0265_0735.log

time mpirun -np 12 $qdyn md_0750_0250.inp > md_0750_0250.log &
time mpirun -np 12 $qdyn md_0250_0750.inp > md_0250_0750.log

time mpirun -np 12 $qdyn md_0762_0237.inp > md_0762_0237.log &
time mpirun -np 12 $qdyn md_0237_0762.inp > md_0237_0762.log

time mpirun -np 12 $qdyn md_0775_0225.inp > md_0775_0225.log &
time mpirun -np 12 $qdyn md_0225_0775.inp > md_0225_0775.log

time mpirun -np 12 $qdyn md_0787_0212.inp > md_0787_0212.log &
time mpirun -np 12 $qdyn md_0212_0787.inp > md_0212_0787.log

time mpirun -np 12 $qdyn md_0800_0200.inp > md_0800_0200.log &
time mpirun -np 12 $qdyn md_0200_0800.inp > md_0200_0800.log

time mpirun -np 12 $qdyn md_0812_0187.inp > md_0812_0187.log &
time mpirun -np 12 $qdyn md_0187_0812.inp > md_0187_0812.log

time mpirun -np 12 $qdyn md_0825_0175.inp > md_0825_0175.log &
time mpirun -np 12 $qdyn md_0175_0825.inp > md_0175_0825.log

time mpirun -np 12 $qdyn md_0837_0162.inp > md_0837_0162.log &
time mpirun -np 12 $qdyn md_0162_0837.inp > md_0162_0837.log

time mpirun -np 12 $qdyn md_0850_0150.inp > md_0850_0150.log &
time mpirun -np 12 $qdyn md_0150_0850.inp > md_0150_0850.log

time mpirun -np 12 $qdyn md_0862_0137.inp > md_0862_0137.log &
time mpirun -np 12 $qdyn md_0137_0862.inp > md_0137_0862.log

time mpirun -np 12 $qdyn md_0875_0125.inp > md_0875_0125.log &
time mpirun -np 12 $qdyn md_0125_0875.inp > md_0125_0875.log

time mpirun -np 12 $qdyn md_0887_0112.inp > md_0887_0112.log &
time mpirun -np 12 $qdyn md_0112_0887.inp > md_0112_0887.log

time mpirun -np 12 $qdyn md_0900_0100.inp > md_0900_0100.log &
time mpirun -np 12 $qdyn md_0100_0900.inp > md_0100_0900.log

time mpirun -np 12 $qdyn md_0912_0087.inp > md_0912_0087.log &
time mpirun -np 12 $qdyn md_0087_0912.inp > md_0087_0912.log

time mpirun -np 12 $qdyn md_0925_0075.inp > md_0925_0075.log &
time mpirun -np 12 $qdyn md_0075_0925.inp > md_0075_0925.log

time mpirun -np 12 $qdyn md_0938_0062.inp > md_0938_0062.log &
time mpirun -np 12 $qdyn md_0062_0938.inp > md_0062_0938.log

time mpirun -np 12 $qdyn md_0947_0053.inp > md_0947_0053.log &
time mpirun -np 12 $qdyn md_0053_0947.inp > md_0053_0947.log

time mpirun -np 12 $qdyn md_0955_0045.inp > md_0955_0045.log &
time mpirun -np 12 $qdyn md_0045_0955.inp > md_0045_0955.log

time mpirun -np 12 $qdyn md_0962_0038.inp > md_0962_0038.log &
time mpirun -np 12 $qdyn md_0038_0962.inp > md_0038_0962.log

time mpirun -np 12 $qdyn md_0967_0033.inp > md_0967_0033.log &
time mpirun -np 12 $qdyn md_0033_0967.inp > md_0033_0967.log

time mpirun -np 12 $qdyn md_0969_0031.inp > md_0969_0031.log &
time mpirun -np 12 $qdyn md_0031_0969.inp > md_0031_0969.log

time mpirun -np 12 $qdyn md_0971_0029.inp > md_0971_0029.log &
time mpirun -np 12 $qdyn md_0029_0971.inp > md_0029_0971.log

time mpirun -np 12 $qdyn md_0974_0026.inp > md_0974_0026.log &
time mpirun -np 12 $qdyn md_0026_0974.inp > md_0026_0974.log

time mpirun -np 12 $qdyn md_0976_0024.inp > md_0976_0024.log &
time mpirun -np 12 $qdyn md_0024_0976.inp > md_0024_0976.log

time mpirun -np 12 $qdyn md_0978_0022.inp > md_0978_0022.log &
time mpirun -np 12 $qdyn md_0022_0978.inp > md_0022_0978.log

time mpirun -np 12 $qdyn md_0979_0021.inp > md_0979_0021.log &
time mpirun -np 12 $qdyn md_0021_0979.inp > md_0021_0979.log

time mpirun -np 12 $qdyn md_0981_0019.inp > md_0981_0019.log &
time mpirun -np 12 $qdyn md_0019_0981.inp > md_0019_0981.log

time mpirun -np 12 $qdyn md_0983_0017.inp > md_0983_0017.log &
time mpirun -np 12 $qdyn md_0017_0983.inp > md_0017_0983.log

time mpirun -np 12 $qdyn md_0985_0015.inp > md_0985_0015.log &
time mpirun -np 12 $qdyn md_0015_0985.inp > md_0015_0985.log

time mpirun -np 12 $qdyn md_0986_0014.inp > md_0986_0014.log &
time mpirun -np 12 $qdyn md_0014_0986.inp > md_0014_0986.log

time mpirun -np 12 $qdyn md_0988_0012.inp > md_0988_0012.log &
time mpirun -np 12 $qdyn md_0012_0988.inp > md_0012_0988.log

time mpirun -np 12 $qdyn md_0989_0011.inp > md_0989_0011.log &
time mpirun -np 12 $qdyn md_0011_0989.inp > md_0011_0989.log

time mpirun -np 12 $qdyn md_0990_0010.inp > md_0990_0010.log &
time mpirun -np 12 $qdyn md_0010_0990.inp > md_0010_0990.log

time mpirun -np 12 $qdyn md_0991_0009.inp > md_0991_0009.log &
time mpirun -np 12 $qdyn md_0009_0991.inp > md_0009_0991.log

time mpirun -np 12 $qdyn md_0992_0008.inp > md_0992_0008.log &
time mpirun -np 12 $qdyn md_0008_0992.inp > md_0008_0992.log

time mpirun -np 12 $qdyn md_0993_0007.inp > md_0993_0007.log &
time mpirun -np 12 $qdyn md_0007_0993.inp > md_0007_0993.log

time mpirun -np 12 $qdyn md_0994_0006.inp > md_0994_0006.log &
time mpirun -np 12 $qdyn md_0006_0994.inp > md_0006_0994.log

time mpirun -np 12 $qdyn md_0995_0005.inp > md_0995_0005.log &
time mpirun -np 12 $qdyn md_0005_0995.inp > md_0005_0995.log

time mpirun -np 12 $qdyn md_0996_0004.inp > md_0996_0004.log &
time mpirun -np 12 $qdyn md_0004_0996.inp > md_0004_0996.log

time mpirun -np 12 $qdyn md_0997_0003.inp > md_0997_0003.log &
time mpirun -np 12 $qdyn md_0003_0997.inp > md_0003_0997.log

time mpirun -np 12 $qdyn md_0998_0002.inp > md_0998_0002.log &
time mpirun -np 12 $qdyn md_0002_0998.inp > md_0002_0998.log

time mpirun -np 12 $qdyn md_0999_0001.inp > md_0999_0001.log &
time mpirun -np 12 $qdyn md_0001_0999.inp > md_0001_0999.log

time mpirun -np 12 $qdyn md_1000_0000.inp > md_1000_0000.log &
time mpirun -np 12 $qdyn md_0000_1000.inp > md_0000_1000.log

timeout 30s /data1/projects/pi-gerard/Q/bin/qfep < qfep.inp > qfep.out
done