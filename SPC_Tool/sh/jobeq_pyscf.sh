#!/bin/bash
#SBATCH -o stdout
#SBATCH --partition=YOUR_PARTITION
#SBATCH -J Dyf_Coeq
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --mail-type=end
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH --time=100:00:00

module purge
# module load AiYH_code/0607-dyf  # install LAMP_emb and set PYTHONPATH instead
#module load loginnode
#export MOLCAS_MEM=50000
#export OMP_NUM_THREADS=56
#export PYTHONPATH=/path/to/LAMP_emb

# source activate pyscf  # activate the environment containing PySCF and LAMP_emb

shopt -s extglob

export TMPDIR=${PWD}
export TMP=${TMPDIR}
export TEMP=${TMPDIR}

date > pyout
python dmet_cassiso.py >> pyout
date >> pyout

