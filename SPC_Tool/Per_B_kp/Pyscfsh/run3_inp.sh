#!/bin/bash

name="DyAC"
num_freq=201
num_pt=10
i_charge=0
prefix="perb"
shopt -s extglob
RTP=$(realpath ../../)
UPP=$PWD

cp $RTP/sh/yhjob_pyscf.sh perb/subjob.sh
cp $RTP/sh/check_pyscf.sh perb/check.sh
cp $RTP/sh/clean_pyscf.sh perb/clean.sh
chmod +x perb/subjob.sh

cd $prefix/eq
python $UPP/insert_xyz.py --tmp $RTP/inpfile/dmet_tmp.py --xyz eq.xyz -o dmet_cassiso.py -c $i_charge -t $name 
chmod +x dmet_cassiso.py
cp $RTP/sh/jobeq_pyscf.sh ./jobeq.sh
cd ../../
echo "---------------------------------eq----------------------------------------"

3for ((imod=22;imod<=$num_freq;imod++));do
for imod in $(seq 1 $num_freq);do 
	cd $prefix/$imod
		cp $RTP/sh/job_pyscf.sh ./pyjob.sh
	cd ../../
	for ipt in $(seq 1 $num_pt); do
		cd $prefix/$imod/$ipt
			python $UPP/insert_xyz.py --tmp $RTP/inpfile/dmet_tmp.py --xyz mol_${imod}_${ipt}.xyz -o dmet_cassiso.py -c $i_charge -t $name
			#rm !(*.py|*.xyz)
			chmod +x dmet_cassiso.py
		cd ../
		#rm !(*.sh)
		#rm -r ompi*
		cd ../../
		echo "--------------modeNo. $imod ------------pointNo. $ipt ---------------------"
	done
done
