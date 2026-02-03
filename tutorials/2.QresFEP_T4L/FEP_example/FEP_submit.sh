#! /bin/bash

temperatures=(298)
runs=10
#restartfile=md_0000_1000.re
workdir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
inputfiles=$workdir/inputfiles
submitfile=$inputfiles/runTETRA.sh

# Random seeds for lysozyme F104A mutation protocol optimization
seeds=(971 2856 8253 1405 5342 7486 5374 1123 8397 8884)

#sed -i s/finalMDrestart=.*/finalMDrestart="$restartfile"/g $submitfile
sed -i s#workdir=.*#workdir="$workdir"#g $submitfile
sed -i s#inputfiles=.*#inputfiles="$inputfiles"#g $submitfile
for temp in ${temperatures[*]};do
sed -i s/temperature=.*/temperature="$temp"/g $submitfile
for i in $(seq 1 $runs);do
sed -i s/run=.*/run="$i"/g $submitfile
sed -i s/seed=1/seed="${seeds[$i-1]}"/g $submitfile
sbatch $submitfile
sed -i s/seed="${seeds[$i-1]}"/seed=1/g $submitfile
done
done
