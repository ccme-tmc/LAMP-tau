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

date
name=mag
cont=0
for ((imod=278;imod>=1;imod--));do
{
	cp job.sh $imod/
	cd $imod/
	icount=0
	for ((ipt=1;ipt<=4;ipt++));do
		cd $ipt/
		tempFile="${name}_${imod}_${ipt}.log"
		if [ ! -f "$tempFile" ];then
		icount=$[icount+1]
		fi
		cd ../
	done
	if [ $icount -gt 0 ];then
		echo $PWD
		echo $icount
		yhbatch job.sh
    	cont=$[cont+1]
	else
		echo "Have submitted mode $imod"
	fi    
	cd ../
	
  if [ $cont -eq 40 ];then
		break
  fi
}
done
echo $cont

date

