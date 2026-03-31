#!/bin/bash
#SBATCH -o stdout
#SBATCH --partition=cps4
#SBATCH -J Dyf_Co2
#SBATCH --get-user-env
##SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mail-type=end
#SBATCH --mail-user=sunhuaiyang@pku.edu.cn
#SBATCH --time=100:00:00

#bin="/THL6/home/jiangh/program/vasp.5.4.4/bin/vasp_std"
module purge
module load AiYH_code/0607-dyf

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
