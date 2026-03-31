bin="/fs2/home/jiangh/softwares/vasp/6.2.0/vasp_std"

module purge
module load vasp/6.2.0
module load vaspkit/1.4.0
module load ase/3.22.1-gcc8.5.0 

CName=Dy

RTP=$(realpath ../../)
OFP=$(realpath ../)

cd ../VASP/opt/

checktext=$(grep "reached required accuracy" stdout)
if [ "$checktext" != "" ];then
	echo "complete opt"
	cd ../../
	#mkdir tmpfile
	cd VASP
	mkdir freq
	cd freq

	cp ../opt/CONTCAR $RTP/inpfile/opt_CONTCAR

	python $OFP/vaspsh/Select_group.py --cry $RTP/inpfile/opt_CONTCAR --ct $CName --out $OFP/tmpfile/ -f
	cp $OFP/tmpfile/freq_POSCAR $RTP/inpfile/
	cp $OFP/tmpfile/group.xyz $RTP/inpfile/eq_opt.xyz
	cp $OFP/tmpfile/dpxs $RTP/inpfile/
	cp $OFP/tmpfile/Atom_list $RTP/inpfile/
	cp $RTP/inpfile/freq_POSCAR ./POSCAR 
	cp $OFP/VASP/opt/POTCAR ./POTCAR
	cp $OFP/VASP/opt/KPOINTS ./KPOINTS
	cp $RTP/inpfile/freq_INCAR ./INCAR
	cp $RTP/sh/job_vasp.sh ./

	date
	sbatch job_vasp.sh
	date
else
	echo "not complete opt"
	cp POSCAR POSCAR_old
	cp CONTCAR POSCAR
	date
	sbatch job_vasp.sh
	date
fi



