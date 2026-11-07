#!/bin/bash
#SBATCH -o stdout
#SBATCH --partition=YOUR_PARTITION
#SBATCH -J Dyf_Co2
#SBATCH --get-user-env
##SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mail-type=end
#SBATCH --mail-user=YOUR_EMAIL
#SBATCH --time=100:00:00

module purge
# module load AiYH_code/0607-dyf  # install LAMP_emb and set PYTHONPATH instead

source activate pyscf
export TMPDIR=${PWD}
export TMP=${TMPDIR}
export TEMP=${TMPDIR}
shopt -s extglob
num_pt=10
exFile="SOC.txtnone"

cp ../eq/*_rohf.chk ./6/
cd 6/    		
	echo $PWD
	if [ ! -f "$exFile" ];then
		date > pyout
		python dmet_cassiso.py >> pyout
		date >> pyout
		rm *_dmet_chk.h5 
	else
		echo "Have finished pt"
	fi
	
cd ../


for ((ipt=7;ipt<=$num_pt;ipt++));do
{
	cd $ipt/
		let tipt=ipt-1
		if [ ! -f "$exFile" ];then
			cp ../$tipt/*_rohf.chk ./
			echo $PWD
			date > pyout
			python dmet_cassiso.py >> pyout
			date >> pyout
			rm *_dmet_chk.h5 
		else
			echo "Have finished pt"
		fi
	cd ../
}
done

cp ../eq/*_rohf.chk ./5/
cd 5/    		
	echo $PWD
	if [ ! -f "$exFile" ];then
		date > pyout
		python dmet_cassiso.py >> pyout
		date >> pyout
		rm *_dmet_chk.h5 
	else
		echo "Have finished pt"
	fi
cd ../

for ((ipt=4;ipt>=1;ipt--));do
{
	cd $ipt/
		let tipt=ipt+1
		if [ ! -f "$exFile" ];then
			cp ../$tipt/*_rohf.chk ./
			echo $PWD
			date > pyout
			python dmet_cassiso.py >> pyout
			date >> pyout
			rm *_dmet_chk.h5 
		else
			echo "Have finished pt"
		fi
	cd ../
}
done
