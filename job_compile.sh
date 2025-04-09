#!/bin/sh
#BATCH --job-name=dwld03_ERA5_MYNN_MYNN
#SBATCH --partition=savio3
#SBATCH --qos=aiolos_savio3_normal 
#SBATCH --account=co_aiolos 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:59:59
#SBATCH --mail-type=all
#SBATCH --mail-user=adamwise@berkeley.edu

mpirun ./compile -j 1 em_les > compile.log

