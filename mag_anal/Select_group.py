# %%
import numpy as np
from numpy import polyfit, poly1d
import argparse
from pathlib import Path
import shutil 
import os 
import sys
import re

import ase as ase
from ase import Atoms
from ase import Atom
from ase.io import read
from ase.visualize import view
from ase.build import fcc111, add_adsorbate
from ase.io import write
from ase.io import vasp

symbols = {
    1: "H",    2: "He",    
    3: "Li",    4: "Be",    5: "B",    6: "C",    7: "N",    8: "O",    9: "F",    10: "Ne",    
    11: "Na",    12: "Mg",    13: "Al",    14: "Si",    15: "P",    16: "S",    17: "Cl",    18: "Ar",    
    19: "K",    20: "Ca",    21: "Sc",    22: "Ti",    23: "V",    24: "Cr",    25: "Mn",    26: "Fe",    27: "Co",    28: "Ni",    29: "Cu",    30: "Zn",    31: "Ga",    32: "Ge",    33: "As",    34: "Se",    35: "Br",    36: "Kr",
    37: "Rb",    38: "Sr",    39: "Y",    40: "Zr",    41: "Nb",    42: "Mo",    43: "Tc",    44: "Ru",    45: "Rh",    46: "Pd",    47: "Ag",    48: "Cd",    49: "In",    50: "Sn",    51: "Sb",    52: "Te",    53: "I",    54: "Xe",    
    55: "Cs",   56: "Ba",    
    57: "La",    58: "Ce",    59: "Pr",    60: "Nd",    61: "Pm",    62: "Sm",    63: "Eu",    64: "Gd",    65: "Tb",    66: "Dy",    67: "Ho",    68: "Er",    69: "Tm",    70: "Yb",    
    71: "Lu",    72: "Hf",    73: "Ta",    74: "W",    75: "Re",    76: "Os",    77: "Ir",    78: "Pt",    79: "Au",    80: "Hg",    81: "Tl",    82: "Pb",    83: "Bi",    84: "Po",    
    85: "At",    86: "Rn",
}


element_dict = {
    "H": 1,    "He": 2,    
    "Li": 3,    "Be": 4,    "B": 5,    "C": 6,    "N": 7,    "O": 8,    "F": 9,    "Ne": 10,    
    "Na": 11,    "Mg": 12,    "Al": 13,    "Si": 14,    "P": 15,    "S": 16,    "Cl": 17,    "Ar": 18,    
    "K": 19,    "Ca": 20,    "Sc": 21,    "Ti": 22,    "V": 23,    "Cr": 24,    "Mn": 25,    "Fe": 26,    "Co": 27,    "Ni": 28,    "Cu": 29,    "Zn": 30,    "Ga": 31,    "Ge": 32,    "As": 33,    "Se": 34,    "Br": 35,    "Kr": 36,
    "Rb": 37,    "Sr": 38,    "Y": 39,    "Zr": 40,    "Nb": 41,    "Mo": 42,    "Tc": 43,    "Ru": 44,    "Rh": 45,    "Pd": 46,    "Ag": 47,    "Cd": 48,    "In": 49,    "Sn": 50,    "Sb": 51,    "Te": 52,    "I": 53,    "Xe": 54,    
    "Cs": 55,    "Ba": 56,    
    "La": 57,    "Ce": 58,    "Pr": 59,    "Nd": 60,    "Pm": 61,    "Sm": 62,    "Eu": 63,    "Gd": 64,    "Tb": 65,    "Dy": 66,    "Ho": 67,    "Er": 68,    "Tm": 69,    "Yb": 70,    
    "Lu": 71,    "Hf": 72,    "Ta": 73,    "W": 74,    "Re": 75,    "Os": 76,    "Ir": 77,    "Pt": 78,    "Au": 79,    "Hg": 80,    "Tl": 81,    "Pb": 82,    "Bi": 83,   
}


def Expand_group(orig,Select_Atnum,length,dpxso=None):
    """
    Use length as radius for each atom in the atom set Select_Atnum, select the atoms within the radius into the atom set, and output the new atom set Select_Atnum. 

    Args:
        orig (atoms): the target molecule.
        Select_Atnum (np): the original atom set.
        length (float): the radius of search
        dpxso (np): It store the imformation of the translation of the atoms. If input it,the function regard the molecule as a cell.

    Returns:
        Select_Atnum (np): the final atom set.
        dpxs
    """
    if dpxso is None:
        for j in Select_Atnum:
            for i in range(0,orig.get_number_of_atoms()):
                l0 = orig.get_positions()[j]
                l1 = orig.get_positions()[i]
                dl = np.sqrt(np.sum((l0-l1)**2))
                #print(i,orig.get_chemical_symbols()[i],l1,dl)
                if dl < length:
                    Select_Atnum = np.append(Select_Atnum,i)
                    #print(dl,length,Select_Atnum)
        Select_Atnum = np.unique(Select_Atnum)
        Select_Atnum = np.sort(Select_Atnum)
        print(Select_Atnum)
        return Select_Atnum
    else:
        dpxs = np.zeros([orig.get_number_of_atoms(),3],dtype = int)
        for j in Select_Atnum:
            for i in range(0,orig.get_number_of_atoms()):
                cell = orig.get_cell()
                l0 = orig.get_positions()[j]
                l1 = orig.get_positions()[i]
                l1,dl,dpx = Select_l(l1,l0,cell)
                #dl = np.sqrt(np.sum((l0-l1)**2))
                if dl < length:
                    Select_Atnum = np.append(Select_Atnum,i)
                    dpxs[i,...] = dpx+dpxso[j,...]
                #if i == 1:
                    #print(dpxs[i,...])
                #dl = np.sqrt((l0[0]-l1[0])**2+(l0[1]-l1[1])**2+(l0[2]-l1[2])**2)
                #print(dl)

        Select_Atnum = np.unique(Select_Atnum)
        Select_Atnum = np.sort(Select_Atnum)
        print(Select_Atnum)
        return Select_Atnum,dpxs

def Select_l(l1,l0,cell):
    """
    Find the nearest distance of two position in the period condition

    Args:
        l0: the position not changed
        l1: the position will changed by the periodic lattice
        cell: the lattice vector

    Returns:
        lx: the nearsest position of l0
        dl: the nearest distance
        dpx: the change by the lattice vector
    """
    lx = l1
    dl = np.sqrt(np.sum((l0-l1)**2))
    dpx = np.zeros([3],dtype = int)
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                l = l1+i*cell[0,...]+j*cell[1,...]+k*cell[2,...]
                dl2 = np.sqrt(np.sum((l0-l)**2))
                if dl2 < dl:
                    lx = l
                    dl = dl2
                    dpx = [i,j,k]
    return lx,dl,dpx

def check_occupy(orig):
    """
    Delete the atoms whose occupancy is lower than 0.5

    Args:
        orig (atoms): the target molecule.

    Returns:
        new (atoms): the final molecule.
    """
    n_cell = int(len(orig)/len(orig.info['occupancy']))
    Select_Atnum = np.zeros([(len(orig))],dtype = int)
    for i in range(len(orig)):
        Select_Atnum[i] = int(i)
    num = Select_Atnum[0]


    positions = tuple(orig.get_positions()[num])
    new = Atoms(orig.get_chemical_symbols()[Select_Atnum[0]],[positions],cell = orig.get_cell())#add the first atom in new
    
    
    j = 1
    for i in range(len(Select_Atnum)-1):
        num = Select_Atnum[i+1]        
        numdict = (i+1)//n_cell
        dit = orig.info['occupancy']["{x}".format(x=numdict)]
        if list(dit.values())[0]>0.5:
            #print(i,num,list(dit.values())[0])
            new.append(orig.get_chemical_symbols()[num])
            new.positions[j] = tuple(orig.get_positions()[num])
            j = j+1
    #print(xray.get_number_of_atoms())
    #view(xray)
    return new

def choose_group(cryfile,center_atom_name,outfold,r1,r2):
    """
    Select the set of atoms from the molecule using this way:
    (1)Choose a center atom
    (2)use Expand_group function with length=r1 to select the Select_Atnum once. 
    (3)use Expand_group function with length=r2 to select the Select_Atnum untill the amount of atoms unchanged or 20 times.
    (4)return the new atoms

    Args:
        orig (atoms): the target molecule.
        center_atom_name (str): the name of the center atom, such as 'Dy'.
        outfold (str): the path of the output fold
        r1 (float)
        r2 (float)

    Returns:
        new (atoms): the final molecule.
    """
    orig = read(cryfile)
    #view(orig)
    cell = orig.get_cell()
    initial_num = np.where(orig.numbers == element_dict[center_atom_name])[0][0]
    Select_Atnum = np.asarray([initial_num])
    cell = orig.get_cell()

    dpxs = np.zeros([orig.get_number_of_atoms(),3],dtype = int)
    #以r1为半径进行第一次选择
    Select_Atnum,dpxs = Expand_group(orig,Select_Atnum,r1,dpxs)
    print("Start select group")

    #以r2为半径不断选择，指导原子集数量不变
    for i in range(0,20):
        num1 = len(Select_Atnum)
        Select_Atnum,dpxs = Expand_group(orig,Select_Atnum,r2,dpxs)
        num2 = len(Select_Atnum)
        if num1 == num2:
            break

    #输出目标原子序号Atom_list，晶格移动编号dpxs
    np.savetxt(outfold+"Atom_list",Select_Atnum)
    np.savetxt(outfold+"dpxs",dpxs)
    print("output dpxs")

    #"""
    """
    dpxs = np.loadtxt("dpxs",dtype = int)
    #print(dpxs)
    list1 = "3   7  11  15  19  23  27  31  35  39  43  47  51  55  59  63  75  87\
        91  95  99 103 111 115 119 123 127 135 139 143 147 159 175 207 211 215\
        219 223 235 239 243 247 251 255 259 263 275 279 283 287 291 319 355"
    li = list1.split()
    Select_Atnum = np.asarray(li,dtype=int)
    #print(Select_Atnum2)
    """
    #按照目标原子编号生成新的结构


    num = Select_Atnum[0]
    #读取初始位置
    positions = tuple(orig.get_positions()[num])+dpxs[num,0]*cell[0,...]+dpxs[num,1]*cell[1,...]+dpxs[num,2]*cell[2,...]
    new = Atoms(orig.get_chemical_symbols()[Select_Atnum[0]],[positions])

    #按照目标原子编号依次在new中添加原子
    for i in range(len(Select_Atnum)-1):
        num = Select_Atnum[i+1]
        new.append(orig.get_chemical_symbols()[num])
        new.positions[i+1] = tuple(orig.get_positions()[num])+dpxs[num,0]*cell[0,...]+dpxs[num,1]*cell[1,...]+dpxs[num,2]*cell[2,...]
    print(new.get_number_of_atoms())

    #输出结构
    #view(new)
    return orig,new,Select_Atnum

def choose_atoms(orig,center_atom_name,choose_atom_name):
    """
    Select the set of atoms including one center_atom and all choose_atom_name atoms.


    Args:
        orig (atoms): the target molecule.
        center_atom_name (str): the name of the center atom, such as 'Dy'.
        choose_atom_name (str): the name of the choosed atom, such as 'C'.

    Returns:
        new (atoms): the final molecule.
    """
    initial_num = np.where(orig.numbers == element_dict[center_atom_name])[0][0]
    atom_distance_list = np.zeros([orig.get_number_of_atoms()],dtype = float)

    for i in range(0,orig.get_number_of_atoms()):
        if orig.get_chemical_symbols()[i] == choose_atom_name:
            #print(orig.get_chemical_symbols()[i])
            l0 = orig.get_positions()[0]
            l1 = orig.get_positions()[i]
            dl = np.sqrt(np.sum((l0-l1)**2))
            atom_distance_list[i] = dl

    little_list = np.argsort(atom_distance_list)
    positions = tuple(orig.get_positions()[initial_num])
    new = Atoms(orig.get_chemical_symbols()[initial_num],[positions])
    j = 0
    for i in range(little_list.shape[0]):
        if j >= 10:
            break
        if atom_distance_list[little_list[i]]>0:
            num = little_list[i]
            new.append(orig.get_chemical_symbols()[num])
            new.positions[j+1] = tuple(orig.get_positions()[num])
            j = j+1

    return new

def addH_to_Cp(orig,r1,rH):
    """
    Add H atoms to the atoms with Dy and tow Cp rings.

    Args:
        orig (atoms): the molecule has one Dy atom and tow Cprings.
        r1 (float): the value, which is larger than the side length of theCprings, but should be lower than the distance of Dy-C. 1.5A is a good choice.
        rH (float): the value of the distance between the C and H, 0.95A is a good choice.

    Returns:
        orig (atoms): the final molecule.
    """
    Cp_Atnum1 = np.zeros([1],dtype = int)
    Cp_Atnum2 = np.zeros([1],dtype = int)

    for i in range(0,orig.get_number_of_atoms()):
        if orig.get_chemical_symbols()[i] == "C":
            Cp_Atnum1[0] = i
            break

    for i in range(0,20):
        num1 = len(Cp_Atnum1)
        Cp_Atnum1 = Expand_group(orig,Cp_Atnum1,r1)
        num2 = len(Cp_Atnum1)
        if num1 == num2:
            break


    for i in range(0,orig.get_number_of_atoms()):
        if orig.get_chemical_symbols()[i] == "C":
            if not np.isin(i, Cp_Atnum1):
                print(i)
                Cp_Atnum2[0] = i
                break

    for i in range(0,20):
        num1 = len(Cp_Atnum2)
        Cp_Atnum2 = Expand_group(orig,Cp_Atnum2,r1)
        num2 = len(Cp_Atnum2)
        if num1 == num2:
            break

    center1 = np.zeros([3],dtype = float)
    for j in Cp_Atnum1:
        l0 = orig.get_positions()[j]
        center1 = center1+l0/Cp_Atnum1.shape
    for j in Cp_Atnum1:
        vector = orig.get_positions()[j]-center1
        vector = vector/np.sqrt(np.sum(vector**2))
        posH = orig.get_positions()[j]+vector*rH
        orig.append(Atom('H',posH))

    center2 = np.zeros([3],dtype = float)
    for j in Cp_Atnum2:
        l0 = orig.get_positions()[j]
        center2 = center2+l0/Cp_Atnum2.shape
    for j in Cp_Atnum2:
        vector = orig.get_positions()[j]-center2
        vector = vector/np.sqrt(np.sum(vector**2))
        posH = orig.get_positions()[j]+vector*rH
        orig.append(Atom('H',posH))
    return orig

def Calculate_mode_list(Atnum):
    num = np.shape(Atnum)[0]
    list = np.zeros([3*num],dtype = int)
    for i in range(0,num):
        for j in range(0,3):
            list[3*i+j] = 3*Atnum[i]+1+j
    return list
    




if __name__ == "__main__":
    """parser = argparse.ArgumentParser(description=""
        A script to select groups from Crystal structure."")
    parser.add_argument("--cry", dest="cryfile",type=str, required=True, help="The filename of Crystal structure.")
    parser.add_argument("--ct", dest="center", type=str, default="None", help="Element of the center of Stru.")
    parser.add_argument("--num",  dest="number", type=int, help="The number of the center atom of Stru.")
    parser.add_argument("--out", dest="outfold",type=str, help="The fold where the out file output")
    parser.add_argument("--freq", "-f", action="store_true", help="output freq-POSCAR.")
    args = parser.parse_args()"""
    #选择目标原子集

    cryfile = "./eq_opt.xyz"
    outfold = "./"
    groupfold = outfold
    center_atom_name = "Dy"
    choose_atom_name = "C"
    r1 = 1.5
    rH = 0.95

    #读取晶体数据
    
    orig = read(cryfile)
    view(orig)

    cell = orig.get_cell()
    
    #确定初始原子集，默认为单个args.center命名的原子，一般为Co或Dy
    #initial_num = np.where(orig.numbers == element_dict['Er'])[0][1]

    #print(initial_num)
    
    #得出晶格参数
    #cell = orig.get_cell()
    #vasp.write_vasp(outfold+"POSCAR",orig,sort=True)
    
    #确定初始原子集，默认为单个args.center命名的原子，一般为Co或Dy
    initial_num = np.where(orig.numbers == element_dict[center_atom_name])[0][0]
    Select_Atnum = np.asarray([initial_num])
    
    
    dpxs = np.zeros([orig.get_number_of_atoms(),3],dtype = int)
    #print(dpxs)
    group,Select_Atnum = choose_group(orig,center_atom_name,outfold,2.9,0.1)
    print(Calculate_mode_list(Select_Atnum))
    view(group)
    #new = choose_atoms(group,center_atom_name,choose_atom_name)


    # with open(outfold+"rtx6_ligs_atoms.txt","w") as f:
    #     f.write("Atom           xx         yy         zz        q     Zeff  group   enum        J\n")
    #     f.write("--------------------------------------------------------------------------------\n")
    #     for i in range(0,new.get_number_of_atoms()):
    #         f.write(new.get_chemical_symbols()[i]+" ")
    #         f.write(" ".join(str(s) for s in new.get_positions()[i]-new.get_positions()[0])+" ")
    #         if new.get_chemical_symbols()[i] == 'Dy':
    #             f.write("+3.000    +3.00      0      9      7.5")
    #         else:
    #             f.write("-0.300    +3.25      1      2      1.0")
    #         f.write("\n")

    # #view(new)
    # write(outfold+"eqx6fs.xyz",new)
    

