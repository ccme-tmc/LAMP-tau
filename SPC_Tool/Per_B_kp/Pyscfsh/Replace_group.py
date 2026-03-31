#!/home/dengyf/.conda/envs/my_env/bin/python
#/usr/bin/python3
# coding: utf-8

# In[14]:


import numpy as np
from numpy import polyfit, poly1d
import argparse
from pathlib import Path
import shutil 
import os 
import sys
import re
#import yaml
#from gen_structures import mol_freq

from ase import Atoms
from ase.io import read
from ase.visualize import view
from ase.build import fcc111, add_adsorbate
from ase.io import write



parser = argparse.ArgumentParser(description="""
    To replace the groups.""")
parser.add_argument("--inp", "-i", help="The template Molcas input file. defaulting to mag_tmp.in", default="mol.xyz")
parser.add_argument("--out", "-o", help="The final input file name, defaulting to mag.in", default="mol.xyz")
parser.add_argument("--inpfile", required=True, help="The foldname of inpfile.")
args = parser.parse_args()


dpxs = np.loadtxt(args.inpfile+"/dpxs")
origcell = read(args.inpfile+"/POSCAR")
cell = origcell.get_cell()

orig = read(args.inp)
list1 = np.loadtxt(args.inpfile+"/Atom_list",dtype = int)
li = list1.split()
Reserve_Atnum = np.asarray(li,dtype=int)
num = Reserve_Atnum[0]
positions = tuple(orig.get_positions()[num])+dpxs[num,0]*cell[0,...]+dpxs[num,1]*cell[1,...]+dpxs[num,2]*cell[2,...]
new = Atoms(orig.get_chemical_symbols()[Reserve_Atnum[0]],[positions])
#view(orig)
for inum in range(len(Reserve_Atnum)-1):
    num = Reserve_Atnum[inum+1]
    new.append(orig.get_chemical_symbols()[num])
    new.positions[inum+1] = tuple(orig.get_positions()[num])+dpxs[num,0]*cell[0,...]+dpxs[num,1]*cell[1,...]+dpxs[num,2]*cell[2,...]
#view(new)
"""
Modify_Atnum = [13,17,18, #SiH3
                14,21,22,
                19,32,43, #Et
                24,50,29]
Modify_Symbol = [1,1,1,1,1,1,1,1,1,1,1,1]
Modify_Length = [1.47,1.47,1.47,1.47,1.47,1.47,1.1,1.1,1.1,1.1,1.1,1.1]
Modify_anchor = [9,9,9,11,11,11,10,20,20,12,23,23]

for inum in range(len(Reserve_Atnum)-1):
    num = Reserve_Atnum[inum+1]
    new.append(orig.get_chemical_symbols()[num])
    new.positions[inum+1] = tuple(orig.get_positions()[num])

for inum in range(len(Modify_Atnum)):
    num = Modify_Atnum[inum]
    new.append(symbols[Modify_Symbol[inum]])
    new.positions[inum+len(Reserve_Atnum)] = tuple(orig.get_positions()[num])
    num_anc = Reserve_Atnum.index(Modify_anchor[inum])
    new.set_distance(num_anc,inum+len(Reserve_Atnum),Modify_Length[inum],fix=0)
"""
write(args.out,new)





# In[ ]:





# In[ ]:




