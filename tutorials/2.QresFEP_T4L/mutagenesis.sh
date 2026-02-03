#!/bin/bash

# Correct PyMOL output formatting of generated PDB files
sed 's/ HA / HA2/g' -i GLY*.pdb
sed -i -E 's/([1-3])(H[A-Z]) / \2\1/g; s/([1-3])(H[A-Z][0-9])/\2\1/g' *.pdb

mkdir -p protein
mkdir -p tripeptide

while read mutation; do
  echo $mutation
  pos=$(echo $mutation | tr -dc '0-9')
  protprep.py -p 2LZM_prep.pdb -r 25 -c RESN:$pos -f OPLSAAM -mc A
  wait

  QresFEP.py -m $mutation -mc A -S protein -t A -d -f OPLSAAM -w 25 -s exponential -l 1 -ts 2fs -T 298 -r 10
  wait
  mv FEP* protein/

  QresFEP.py -m $mutation -mc A -S water -t A -d -f OPLSAAM -w 25 -s exponential -l 1 -ts 2fs -T 298 -r 10
  wait
  mv FEP* tripeptide/
done < mutations_neutral.txt
