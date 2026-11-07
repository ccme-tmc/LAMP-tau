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
num_freq=216
num_pt=10
name=mag
shopt -s extglob
for ((imod=1;imod<=$num_freq;imod++));do
{
	cd $imod/
	for ((ipt=1;ipt<=$num_pt;ipt++));do
	cd $ipt/
	tempFile="${name}_${imod}_${ipt}.log"
	if [ ! -f "$tempFile" ];then
		echo "Not calculated mode $imod point ${ipt}"
		#yhbatch job.sh
	else
		checktext=$(tail -n 20 $tempFile | grep "Happy landing!")
		if [ "$checktext" != "" ];then
			echo "complete mode $imod point ${ipt}"
				echo $PWD
				tempFile2="${name}_${imod}_${ipt}.log"
				mv $tempFile2 temp.tep				
				sed -n '1,500p' temp.tep > temp1.tep
				sed -n '/MULTIPLET 1/,$p' temp.tep > temp2.tep
				cat temp1.tep temp2.tep > $tempFile2
				rm -rf !(*.xyz|*.in|*.log)
			echo "Clear mode $imod"
		else
			echo "not complete mode $imod point ${ipt}"
			echo $PWD
				#rm -rf !(*.xyz|*.in)
		fi
	fi    
	cd ../
	done
	cd ../
}
done


date

