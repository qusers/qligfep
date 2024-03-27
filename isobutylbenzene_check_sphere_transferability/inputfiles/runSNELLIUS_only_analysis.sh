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

## define qdynp location
qdyn=/home/wjespers/software/Q/bin/qdynp
fepfiles=("FEP1.fep" "FEP2.fep" "FEP3.fep")
temperature=298
run=10
finalMDrestart=md_0000_1000.re

workdir=/gpfs/scratch1/nodespecific/int6/yricky/softcore_with_long_endpoint_sampling/butylbenzene_flag_0_seq5_soft5_15_radius_eq_15
inputfiles=/gpfs/scratch1/nodespecific/int6/yricky/softcore_with_long_endpoint_sampling/butylbenzene_flag_0_seq5_soft5_15_radius_eq_15/inputfiles

cd /gpfs/scratch1/nodespecific/int6/yricky
python retreive_dG_BAR_and_restraints.py -F two_restrained_types_FEP23_25_sphere/benzene/FEP1
python retreive_dG_BAR_and_restraints.py -F two_restrained_types_FEP23_25_sphere/benzene/FEP2
python retreive_dG_BAR_and_restraints.py -F two_restrained_types_FEP23_25_sphere/benzene/FEP3
python retreive_shellrestr.py -F two_restrained_types_FEP23_25_sphere/benzene/FEP1
python retreive_shellrestr.py -F two_restrained_types_FEP23_25_sphere/benzene/FEP2
python retreive_shellrestr.py -F two_restrained_types_FEP23_25_sphere/benzene/FEP3
python retreive_total_restr.py -F two_restrained_types_FEP23_25_sphere/benzene/FEP1
python retreive_total_restr.py -F two_restrained_types_FEP23_25_sphere/benzene/FEP2
python retreive_total_restr.py -F two_restrained_types_FEP23_25_sphere/benzene/FEP3
~/anaconda3/envs/QligFEP/bin/python check_angle.py -F two_restrained_types_FEP23_25_sphere/benzene/FEP1 -I inputfiles_softcore_check/benzene -FF OPLS2005 -convert
~/anaconda3/envs/QligFEP/bin/python check_angle.py -F two_restrained_types_FEP23_25_sphere/benzene/FEP2 -I inputfiles_softcore_check/benzene -FF OPLS2005 -convert
~/anaconda3/envs/QligFEP/bin/python check_angle.py -F two_restrained_types_FEP23_25_sphere/benzene/FEP3 -I inputfiles_softcore_check/benzene -FF OPLS2005 -convert