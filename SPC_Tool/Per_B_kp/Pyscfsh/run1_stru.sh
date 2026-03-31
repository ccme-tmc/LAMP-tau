#!/usr/bin/bash
module load ase/3.22.1-gcc8.5.0 


python3 gen_structures_xyz.py --pos ../../inpfile/eq_opt.xyz -np 10 -mdp 0.02 
#python3 gen_structures_cell.py --out inpfile/freq.log --pos ../inpfile/opt_CONTCAR -np 4 -mdp 0.2