#!/bin/bash
software=$PWD/../../software
for XXXX in */ ; do
cd $XXXX
    echo "#!/bin/bash -l
    #source /home/apps/gromacs-4.6.7/bin/GMXRC.bash
    #source /home/apps/bin/apps.sh
    module load gromacs/4.6.7
    python2 $software/qligfep/scripts/pymemdyn.py
    " > temp.sh
    chmod +x temp.sh
    sbatch  -c 8 -t 47:59:00 -J pymemdyn temp.sh
    cd ../
done