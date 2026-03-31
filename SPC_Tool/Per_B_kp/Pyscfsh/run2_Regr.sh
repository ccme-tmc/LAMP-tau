#!/bin/bash

num_freq=216
num_pt=10
prefix="perb"

cd $prefix/eq
python ../../Replace_group.py -i "eq.xyz"  -o "eq2.xyz" --inpfile "../../../inpfile"
cd ../../
echo "---------------------------------eq----------------------------------------"

for imod in $(seq 1 $num_freq);do
	for ipt in $(seq 1 $num_pt); do
		cd $prefix/$imod/$ipt
			python ../../../Replace_group.py -i "mol_${imod}_${ipt}.xyz"  -o "mol2_${imod}_${ipt}.xyz" --inpfile "../../../../inpfile"
		cd ../../../
		echo "--------------modeNo. $imod ------------pointNo. $ipt ---------------------"
	done
done
