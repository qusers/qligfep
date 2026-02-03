#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --job-name=QresFEP 
#SBATCH -A naiss2024-3-13
#SBATCH -p tetralith
#              d-hh:mm:ss
#SBATCH --time=0-05:00:00

## Load modules for qdynp



## define qdynp location
qdyn=/proj/alexandria/users/x_lucko/software/q6/bin/qdynp
fepfiles=(FEP1.fep FEP2.fep)
temperature=298
run=10
finalMDrestart=md_1_0000_1000.re

seed=1

workdir=/proj/alexandria/users/x_lucko/hybrid_fep/1YPC/protein/FEP_ARG62ALA
inputfiles=/proj/alexandria/users/x_lucko/hybrid_fep/1YPC/protein/FEP_ARG62ALA/inputfiles

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
time srun $qdyn eq1.inp > eq1.log
wait
time srun $qdyn eq2.inp > eq2.log
wait
time srun $qdyn eq3.inp > eq3.log
wait
time srun $qdyn eq4.inp > eq4.log
wait
time srun $qdyn eq5.inp > eq5.log
wait
#RUN_1_FILES
echo md_1_1000_0000
time srun $qdyn md_1_1000_0000.inp > md_1_1000_0000.log
wait
echo md_1_0997_0003
time srun $qdyn md_1_0997_0003.inp > md_1_0997_0003.log
wait
echo md_1_0993_0007
time srun $qdyn md_1_0993_0007.inp > md_1_0993_0007.log
wait
echo md_1_0989_0011
time srun $qdyn md_1_0989_0011.inp > md_1_0989_0011.log
wait
echo md_1_0985_0015
time srun $qdyn md_1_0985_0015.inp > md_1_0985_0015.log
wait
echo md_1_0981_0019
time srun $qdyn md_1_0981_0019.inp > md_1_0981_0019.log
wait
echo md_1_0977_0023
time srun $qdyn md_1_0977_0023.inp > md_1_0977_0023.log
wait
echo md_1_0972_0028
time srun $qdyn md_1_0972_0028.inp > md_1_0972_0028.log
wait
echo md_1_0967_0033
time srun $qdyn md_1_0967_0033.inp > md_1_0967_0033.log
wait
echo md_1_0961_0039
time srun $qdyn md_1_0961_0039.inp > md_1_0961_0039.log
wait
echo md_1_0956_0044
time srun $qdyn md_1_0956_0044.inp > md_1_0956_0044.log
wait
echo md_1_0950_0050
time srun $qdyn md_1_0950_0050.inp > md_1_0950_0050.log
wait
echo md_1_0943_0057
time srun $qdyn md_1_0943_0057.inp > md_1_0943_0057.log
wait
echo md_1_0936_0064
time srun $qdyn md_1_0936_0064.inp > md_1_0936_0064.log
wait
echo md_1_0929_0071
time srun $qdyn md_1_0929_0071.inp > md_1_0929_0071.log
wait
echo md_1_0921_0079
time srun $qdyn md_1_0921_0079.inp > md_1_0921_0079.log
wait
echo md_1_0913_0087
time srun $qdyn md_1_0913_0087.inp > md_1_0913_0087.log
wait
echo md_1_0904_0096
time srun $qdyn md_1_0904_0096.inp > md_1_0904_0096.log
wait
echo md_1_0895_0105
time srun $qdyn md_1_0895_0105.inp > md_1_0895_0105.log
wait
echo md_1_0885_0115
time srun $qdyn md_1_0885_0115.inp > md_1_0885_0115.log
wait
echo md_1_0874_0126
time srun $qdyn md_1_0874_0126.inp > md_1_0874_0126.log
wait
echo md_1_0863_0137
time srun $qdyn md_1_0863_0137.inp > md_1_0863_0137.log
wait
echo md_1_0851_0149
time srun $qdyn md_1_0851_0149.inp > md_1_0851_0149.log
wait
echo md_1_0838_0162
time srun $qdyn md_1_0838_0162.inp > md_1_0838_0162.log
wait
echo md_1_0825_0175
time srun $qdyn md_1_0825_0175.inp > md_1_0825_0175.log
wait
echo md_1_0810_0190
time srun $qdyn md_1_0810_0190.inp > md_1_0810_0190.log
wait
echo md_1_0795_0205
time srun $qdyn md_1_0795_0205.inp > md_1_0795_0205.log
wait
echo md_1_0779_0221
time srun $qdyn md_1_0779_0221.inp > md_1_0779_0221.log
wait
echo md_1_0761_0239
time srun $qdyn md_1_0761_0239.inp > md_1_0761_0239.log
wait
echo md_1_0743_0257
time srun $qdyn md_1_0743_0257.inp > md_1_0743_0257.log
wait
echo md_1_0724_0276
time srun $qdyn md_1_0724_0276.inp > md_1_0724_0276.log
wait
echo md_1_0703_0297
time srun $qdyn md_1_0703_0297.inp > md_1_0703_0297.log
wait
echo md_1_0681_0319
time srun $qdyn md_1_0681_0319.inp > md_1_0681_0319.log
wait
echo md_1_0657_0343
time srun $qdyn md_1_0657_0343.inp > md_1_0657_0343.log
wait
echo md_1_0632_0368
time srun $qdyn md_1_0632_0368.inp > md_1_0632_0368.log
wait
echo md_1_0606_0394
time srun $qdyn md_1_0606_0394.inp > md_1_0606_0394.log
wait
echo md_1_0578_0422
time srun $qdyn md_1_0578_0422.inp > md_1_0578_0422.log
wait
echo md_1_0548_0452
time srun $qdyn md_1_0548_0452.inp > md_1_0548_0452.log
wait
echo md_1_0516_0484
time srun $qdyn md_1_0516_0484.inp > md_1_0516_0484.log
wait
echo md_1_0482_0518
time srun $qdyn md_1_0482_0518.inp > md_1_0482_0518.log
wait
echo md_1_0446_0554
time srun $qdyn md_1_0446_0554.inp > md_1_0446_0554.log
wait
echo md_1_0408_0592
time srun $qdyn md_1_0408_0592.inp > md_1_0408_0592.log
wait
echo md_1_0367_0633
time srun $qdyn md_1_0367_0633.inp > md_1_0367_0633.log
wait
echo md_1_0324_0676
time srun $qdyn md_1_0324_0676.inp > md_1_0324_0676.log
wait
echo md_1_0278_0722
time srun $qdyn md_1_0278_0722.inp > md_1_0278_0722.log
wait
echo md_1_0229_0771
time srun $qdyn md_1_0229_0771.inp > md_1_0229_0771.log
wait
echo md_1_0177_0823
time srun $qdyn md_1_0177_0823.inp > md_1_0177_0823.log
wait
echo md_1_0121_0879
time srun $qdyn md_1_0121_0879.inp > md_1_0121_0879.log
wait
echo md_1_0062_0938
time srun $qdyn md_1_0062_0938.inp > md_1_0062_0938.log
wait
echo md_1_0000_1000
time srun $qdyn md_1_0000_1000.inp > md_1_0000_1000.log
wait
else
#RUN_2_FILES
echo md_2_1000_0000
time srun $qdyn md_2_1000_0000.inp > md_2_1000_0000.log
wait
echo md_2_0938_0062
time srun $qdyn md_2_0938_0062.inp > md_2_0938_0062.log
wait
echo md_2_0879_0121
time srun $qdyn md_2_0879_0121.inp > md_2_0879_0121.log
wait
echo md_2_0823_0177
time srun $qdyn md_2_0823_0177.inp > md_2_0823_0177.log
wait
echo md_2_0771_0229
time srun $qdyn md_2_0771_0229.inp > md_2_0771_0229.log
wait
echo md_2_0722_0278
time srun $qdyn md_2_0722_0278.inp > md_2_0722_0278.log
wait
echo md_2_0676_0324
time srun $qdyn md_2_0676_0324.inp > md_2_0676_0324.log
wait
echo md_2_0633_0367
time srun $qdyn md_2_0633_0367.inp > md_2_0633_0367.log
wait
echo md_2_0592_0408
time srun $qdyn md_2_0592_0408.inp > md_2_0592_0408.log
wait
echo md_2_0554_0446
time srun $qdyn md_2_0554_0446.inp > md_2_0554_0446.log
wait
echo md_2_0518_0482
time srun $qdyn md_2_0518_0482.inp > md_2_0518_0482.log
wait
echo md_2_0484_0516
time srun $qdyn md_2_0484_0516.inp > md_2_0484_0516.log
wait
echo md_2_0452_0548
time srun $qdyn md_2_0452_0548.inp > md_2_0452_0548.log
wait
echo md_2_0422_0578
time srun $qdyn md_2_0422_0578.inp > md_2_0422_0578.log
wait
echo md_2_0394_0606
time srun $qdyn md_2_0394_0606.inp > md_2_0394_0606.log
wait
echo md_2_0368_0632
time srun $qdyn md_2_0368_0632.inp > md_2_0368_0632.log
wait
echo md_2_0343_0657
time srun $qdyn md_2_0343_0657.inp > md_2_0343_0657.log
wait
echo md_2_0319_0681
time srun $qdyn md_2_0319_0681.inp > md_2_0319_0681.log
wait
echo md_2_0297_0703
time srun $qdyn md_2_0297_0703.inp > md_2_0297_0703.log
wait
echo md_2_0276_0724
time srun $qdyn md_2_0276_0724.inp > md_2_0276_0724.log
wait
echo md_2_0257_0743
time srun $qdyn md_2_0257_0743.inp > md_2_0257_0743.log
wait
echo md_2_0239_0761
time srun $qdyn md_2_0239_0761.inp > md_2_0239_0761.log
wait
echo md_2_0221_0779
time srun $qdyn md_2_0221_0779.inp > md_2_0221_0779.log
wait
echo md_2_0205_0795
time srun $qdyn md_2_0205_0795.inp > md_2_0205_0795.log
wait
echo md_2_0190_0810
time srun $qdyn md_2_0190_0810.inp > md_2_0190_0810.log
wait
echo md_2_0175_0825
time srun $qdyn md_2_0175_0825.inp > md_2_0175_0825.log
wait
echo md_2_0162_0838
time srun $qdyn md_2_0162_0838.inp > md_2_0162_0838.log
wait
echo md_2_0149_0851
time srun $qdyn md_2_0149_0851.inp > md_2_0149_0851.log
wait
echo md_2_0137_0863
time srun $qdyn md_2_0137_0863.inp > md_2_0137_0863.log
wait
echo md_2_0126_0874
time srun $qdyn md_2_0126_0874.inp > md_2_0126_0874.log
wait
echo md_2_0115_0885
time srun $qdyn md_2_0115_0885.inp > md_2_0115_0885.log
wait
echo md_2_0105_0895
time srun $qdyn md_2_0105_0895.inp > md_2_0105_0895.log
wait
echo md_2_0096_0904
time srun $qdyn md_2_0096_0904.inp > md_2_0096_0904.log
wait
echo md_2_0087_0913
time srun $qdyn md_2_0087_0913.inp > md_2_0087_0913.log
wait
echo md_2_0079_0921
time srun $qdyn md_2_0079_0921.inp > md_2_0079_0921.log
wait
echo md_2_0071_0929
time srun $qdyn md_2_0071_0929.inp > md_2_0071_0929.log
wait
echo md_2_0064_0936
time srun $qdyn md_2_0064_0936.inp > md_2_0064_0936.log
wait
echo md_2_0057_0943
time srun $qdyn md_2_0057_0943.inp > md_2_0057_0943.log
wait
echo md_2_0050_0950
time srun $qdyn md_2_0050_0950.inp > md_2_0050_0950.log
wait
echo md_2_0044_0956
time srun $qdyn md_2_0044_0956.inp > md_2_0044_0956.log
wait
echo md_2_0039_0961
time srun $qdyn md_2_0039_0961.inp > md_2_0039_0961.log
wait
echo md_2_0033_0967
time srun $qdyn md_2_0033_0967.inp > md_2_0033_0967.log
wait
echo md_2_0028_0972
time srun $qdyn md_2_0028_0972.inp > md_2_0028_0972.log
wait
echo md_2_0023_0977
time srun $qdyn md_2_0023_0977.inp > md_2_0023_0977.log
wait
echo md_2_0019_0981
time srun $qdyn md_2_0019_0981.inp > md_2_0019_0981.log
wait
echo md_2_0015_0985
time srun $qdyn md_2_0015_0985.inp > md_2_0015_0985.log
wait
echo md_2_0011_0989
time srun $qdyn md_2_0011_0989.inp > md_2_0011_0989.log
wait
echo md_2_0007_0993
time srun $qdyn md_2_0007_0993.inp > md_2_0007_0993.log
wait
echo md_2_0003_0997
time srun $qdyn md_2_0003_0997.inp > md_2_0003_0997.log
wait
echo md_2_0000_1000
time srun $qdyn md_2_0000_1000.inp > md_2_0000_1000.log
wait
fi
printf "%s\n" $(ls -1 md*en) | tac >> qfep.inp
/proj/alexandria/users/x_lucko/software/q6/bin/qfep < qfep.inp > qfep.out
#rm *.log
#rm *.dcd
#mv md_1_0000_1000.re md_1_0000_1000.re.keep
#rm *.re
#mv md_1_0000_1000.re.keep md_1_0000_1000.re
done
