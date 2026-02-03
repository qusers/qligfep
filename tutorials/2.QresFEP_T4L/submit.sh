#!/bin/bash

dir=$(pwd)
for fep in protein/FEP_*/
do
        echo $fep
        cd $fep
        bash FEP_submit.sh
        cd $dir
done
for fep in tripeptide/FEP_*/
do
         echo $fep
        cd $fep
        bash FEP_submit.sh
        cd $dir
done
