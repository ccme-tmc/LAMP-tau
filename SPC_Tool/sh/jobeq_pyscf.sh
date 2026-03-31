#!/bin/bash
#SBATCH -o stdout
#SBATCH --partition=cp4
#SBATCH -J Dyf_Coeq
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --mail-type=end
#SBATCH --mail-user=sunhuaiyang@pku.edu.cn
#SBATCH --time=100:00:00

#bin="/THL6/home/jiangh/program/vasp.5.4.4/bin/vasp_std"
module purge
module load AiYH_code/0607-dyf
#module load loginnode
#export MOLCAS_MEM=50000
#export OMP_NUM_THREADS=56
#export PYTHONPATH=/fs2/home/jiangh/program/AiYH_Code_20230607/

source activate pyscf

shopt -s extglob

export TMPDIR=${PWD}
export TMP=${TMPDIR}
export TEMP=${TMPDIR}

date > pyout
python dmet_cassiso.py >> pyout
date >> pyout

