#!/bin/bash
#SBATCH -o stdout
#SBATCH --partition=YOUR_PARTITION
#SBATCH -J Dyf_DyCp
#SBATCH --get-user-env
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=56
#SBATCH --mail-type=end
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH --time=100:00:00

bin="${VASP_BIN:-vasp_std}"

module purge
module load vasp/6.2.0

date
mpirun -np 112 $bin
date
