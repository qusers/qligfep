#!/bin/bash

dir=$(pwd)
[ -f protein/analyze_protein.txt ] && rm protein/analyze_protein.txt
for fep in protein/FEP_*/
do
        echo "${fep#protein/}"
	cd protein
        analyze_FEP.py -F "${fep#protein/}" -T 298 -l 1 >> analyze_protein.txt
	cd $dir
done

[ -f tripeptide/analyze_tripeptide.txt ] && rm tripeptide/analyze_tripeptide.txt
for fep in tripeptide/FEP_*/
do
        echo "${fep#tripeptide/}"
        cd tripeptide
        analyze_FEP.py -F "${fep#tripeptide/}" -T 298 -l 1 >> analyze_tripeptide.txt
        cd $dir
done
