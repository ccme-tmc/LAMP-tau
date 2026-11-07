#!/bin/bash
#SBATCH -o stdout
#SBATCH --partition=YOUR_PARTITION
#SBATCH -J Dyf_Dyac
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --mail-type=end
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH --time=100:00:00

module purge
module load openmolcas/v23.02-dyf
export MOLCAS_MEM=50000
export OMP_NUM_THREADS=14
export MOLCAS_KEEP_WORKDIR=NO
#export Project=check

date
pymolcas eq_exp.in -f
date

