#!/bin/bash
#SBATCH -o stdout
#SBATCH --partition=cp4
#SBATCH -J Dyf_Dyac
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --mail-type=end
#SBATCH --mail-user=sunhuaiyang@pku.edu.cn
#SBATCH --time=100:00:00

#bin="/THL6/home/jiangh/program/vasp.5.4.4/bin/vasp_std"
module purge
module load openmolcas/v23.02-dyf
export MOLCAS_MEM=50000
export OMP_NUM_THREADS=14
export MOLCAS_KEEP_WORKDIR=NO
#export Project=check

date
name=mag
imod=${PWD#*perb2/}
cd ../
for ((ipt=1;ipt<=4;ipt++));do
{
	cd $imod/$ipt/
    #echo $PWD
	tempFile="${name}_${imod}_${ipt}.log"
	if [ ! -f "$tempFile" ];then
		pymolcas ${name}_${imod}_${ipt}.in -f
	fi
	cd ../../
}&
done
wait

date

