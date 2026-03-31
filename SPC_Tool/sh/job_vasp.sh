#!/bin/bash
#SBATCH -o stdout
#SBATCH --partition=cp4
#SBATCH -J Dyf_DyCp
#SBATCH --get-user-env
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=56
#SBATCH --mail-type=end
#SBATCH --mail-user=sunhuaiyang@pku.edu.cn
#SBATCH --time=100:00:00

bin="/fs2/home/jiangh/softwares/vasp/6.2.0/vasp_std"

module purge
module load vasp/6.2.0

date
mpirun -np 112 $bin
date
