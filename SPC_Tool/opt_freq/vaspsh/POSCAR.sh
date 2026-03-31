bin="/fs2/home/jiangh/softwares/vasp/6.2.0/vasp_std"

module purge

module load ase/3.22.1-gcc8.5.0 

mkdir ../tmpfile

RTP=$(realpath ../../)
OFP=$(realpath ../)

python gen_POSCAR.py --cry $RTP/inpfile/.cif --out $OFP/tmpfile/
cp $OFP/tmpfile/POSCAR $RTP/inpfile/POSCAR



