bin="/fs2/home/jiangh/softwares/vasp/6.2.0/vasp_std"

module purge
module load vasp/6.2.0
module load vaspkit/1.4.0

RTP=$(realpath ../../)
OFP=$(realpath ../)
cd ../
mkdir tmpfile
mkdir VASP
cd VASP
mkdir opt
cd opt

cp $RTP/inpfile/POSCAR ./
cp $RTP/inpfile/opt_INCAR ./INCAR
cp $RTP/sh/job_vasp.sh ./
vaspkit -task 103
vaspkit -task 102 -kpr 0.03

date
sbatch job_vasp.sh
date
