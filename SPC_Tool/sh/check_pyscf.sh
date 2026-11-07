#!/bin/bash
#SBATCH -o stdout
#SBATCH --partition=YOUR_PARTITION
#SBATCH -J Dyf_Dyd2
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=end
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH --time=100:00:00

#module purge
#module load openmolcas/v22.06-dyf
export MOLCAS_MEM=50000
export OMP_NUM_THREADS=7
export MOLCAS_KEEP_WORKDIR=NO
#export Project=check


cont=0
num_freq=279
name="DyCp"
shopt -s extglob

for ((imod=1;imod<=$num_freq;imod++));do
{
	cd $imod/
	tempFile="1/${name}_rohf.chk"
	if [ ! -f "$tempFile" ];then
		echo $PWD
		#cp ../eq/*_rohf.chk ./6/
		#cp ../eq/*_rohf.chk ./5/
		#qsub job.sh  
	else
		echo "Have submited mode $imod"
	fi
	cd ../
}
done


