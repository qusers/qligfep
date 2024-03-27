#! /bin/bash

temperatures=(298)
runs=10
restartfile=md_0000_1000.re
workdir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
inputfiles=$workdir/inputfiles
submitfile=$inputfiles/runSNELLIUS.sh

sed -i s/finalMDrestart=.*/finalMDrestart="$restartfile"/g $submitfile
sed -i s#workdir=.*#workdir="$workdir"#g $submitfile
sed -i s#inputfiles=.*#inputfiles="$inputfiles"#g $submitfile
for temp in ${temperatures[*]};do
sed -i s/temperature=.*/temperature="$temp"/g $submitfile
for i in $(seq 1 $runs);do
sed -i s/run=.*/run="$i"/g $submitfile
sbatch $submitfile
done
done
