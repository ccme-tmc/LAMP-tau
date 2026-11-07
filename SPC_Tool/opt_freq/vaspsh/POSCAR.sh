bin="${VASP_BIN:-vasp_std}"

module purge

module load ase/3.22.1-gcc8.5.0 

mkdir -p ../tmpfile

RTP=$(realpath ../../)
OFP=$(realpath ../)

python gen_POSCAR.py --cry $RTP/inpfile/.cif --out $OFP/tmpfile/
cp $OFP/tmpfile/POSCAR $RTP/inpfile/POSCAR



