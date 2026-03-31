#!/home/dengyf/.conda/envs/my_env/bin/python
#/THL6/home/jiangh/members/lisx/install_lisx/anaconda3/bin/python
# -*- coding:utf-8 -*-

import numpy as np
from numpy import polyfit, poly1d
import argparse
from pathlib import Path
import shutil 
import os 
import sys
import re
#import yaml
#from matplotlib import pyplot as plt

from ase import Atoms
from ase.io import read
from ase.visualize import view
from ase.build import fcc111, add_adsorbate
from ase.io import write
from ase.io import vasp


##################################### constants ####################################
# conversion factor from m (reduced mass, amu) and w (frequency, cm^-1) to the amplitude of a harmonic oscillator x ( Angstrom)
hbar = 1.05457266
amu = 1.66053886
c = 2.9979
x_fac = np.sqrt( hbar / (amu * 2.0 * np.pi * c) ) * 10 ** 1.5
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

atMass = {
    'H': 1.008,    'He': 4.003,    
    'Li': 6.941,    'Be': 9.012,    'B': 10.81,     'C': 12.01,    'N': 14.01,    'O': 16.00,     'F': 19.00,     'Ne': 20.18,     
    'Na': 22.99,    'Mg': 24.31,    'Al': 26.98,    'Si': 28.09,    'P': 30.97,    'S': 32.07,    'Cl': 35.45,    'Ar': 39.95,    
    'K': 39.10, 'Ca': 40.08,    'Sc': 44.96,    'Ti': 47.87,    'V': 50.94,    'Cr': 52.00,    'Mn': 54.94,    'Fe': 55.85,    'Ni': 58.69,    'Co': 58.93,    'Cu': 63.55,    'Zn': 65.38,    'Ga': 69.72,    'Ge': 72.63,    'As': 74.92,    'Se': 78.96,    'Br': 79.90,    'Kr': 83.80,
    'Rb': 85.47,    'Sr': 87.62,    'Y': 88.91,    'Zr': 91.22,    'Nb': 92.91,    'Mo': 95.94,    'Tc': 98.00,    'Ru': 101.07,    'Rh': 102.91,    'Pd': 106.42,    'Ag': 107.87,    'Cd': 112.41,    'In': 114.82,    'Sn': 118.71,    'Sb': 121.76,    'I': 126.90,    'Te': 127.60,    'Xe': 131.29,
    'Cs': 132.91,    'Ba': 137.33,
    'La': 138.91,    'Ce': 140.12,    'Pr': 140.91,    'Nd': 144.24,    'Pm': 145.00,    'Sm': 150.36,    'Eu': 151.96,    'Gd': 157.25,    'Tb': 158.93,    'Dy': 162.50,    'Ho': 164.93,    'Er': 167.26,    'Tm': 168.93,    'Yb': 173.05,
    'Lu': 174.97,    'Hf': 178.49,    'Ta': 180.95,    'W': 183.84,    'Re': 186.21,    'Os': 190.23,    'Ir': 192.22,    'Pt': 195.08,    'Au': 196.97,    'Hg': 200.59,    'Tl': 204.38,    'Pb': 207.2,    'Bi': 208.98,
}

def count_string(string, filename):
    """
        return the number of appearances of string in the file filename
    """
    cmd = "n=`grep " + '"' + string + '" ' + filename + " | wc -l` && echo $n > tmp.dat"
    os.system(cmd)
    with open("tmp.dat") as f:
        n = int(f.readline().split()[0])
    os.system("rm tmp.dat")
    return n 

class mol_freq(object):
    """
        the object corresponds to freq calculation of one molecule 
        in Gaussian.
    """ 
    def __init__(self, outfile,posfile):
        self.outfile = outfile
        self.posfile = posfile
    
    def read_geometry(self):
        """
            read in the equilibrium geometry from the associated log file
        """
        self.orig = read(self.posfile)
        self.crds = np.array(self.orig.positions)
    
    def read_modes(self):
        """
            read in normal modes information from the output file
            including the frequencies (cm-1), reduced masses (amu), force constants (mDyne/A) 
            and atom displacement vectors (unitless). 
        """
        n_freq_analysis = count_string("cm-1", self.outfile)
        freqs = np.empty([n_freq_analysis], dtype = float)
        vecs = np.empty([n_freq_analysis,self.orig.get_number_of_atoms(),3], dtype = float)
        print("""
There are %4d frequency analysis in the log file, and we always use the first part.
When this number is greater than 1, the first part contains more digit numbers for the 
displacement vectors. 
""" % n_freq_analysis) 
        if n_freq_analysis == 1:
            print("\nThe vectors have only 2 digits, more accurate results could be printed in log by specifying freq=hpmodes in Gaussian input file!\n")      
        with open(self.outfile, "r") as f:
            for i in range(0,n_freq_analysis):
                line = f.readline()
                freqs[i] = float(line.split()[7])
                line = f.readline()
                for j in range(0,self.orig.get_number_of_atoms()):
                    vecs[i,j,...] = f.readline().split()[3:]
                line = f.readline()
        self.freqs, self.vecs = freqs, vecs 
        umasses = np.empty([n_freq_analysis], dtype = float)
        for i in range(0,n_freq_analysis):
            mass = 0
            for j in range(0,self.orig.get_number_of_atoms()):
                mass += (self.vecs[i,j,...] ** 2).sum()/atMass[symbols[self.orig.get_atomic_numbers()[j]]]
            umasses[i] = 1/mass
        self.masses = umasses

    def perturb(self, npts, mdp):
        """
            to perturb the equilibrium structures along each normal modes. 
            The maximum displacements are the amplitude of the zero-point vibration.
            The number of points are controled by npts. 
        """
        if npts % 2 != 0:
            print("The number of points along a given mode should be an odd number!")
            return
        # check if the output directory already exists
        if Path("./perb").exists():
            print("./perb directory already exist! We now delete it...")
            shutil.rmtree("./perb") 
        Path("perb").mkdir(parents=True)
        scalings = list(np.linspace(-mdp, mdp, npts + 1)) 
        scalings.pop( npts // 2) # remove the middle 0 
        scalings = np.array(scalings)
        Path("./perb/eq").mkdir(parents=True) 
        outfile = "./perb/eq/eq.xyz"
        new_crd = self.crds
        with open(outfile, "w") as f:
            # print out coordinates 
            print(self.orig.get_number_of_atoms(), file=f)
            print("eq", file=f) 
            for i in range(self.orig.get_number_of_atoms()):
                print(" %6s %12.6f %12.6f %12.6f" % (symbols[self.orig.get_atomic_numbers()[i]], new_crd[i,0], new_crd[i,1], new_crd[i, 2]), file=f)
        for imod in range(self.freqs.size):
            Path("./perb/" + str(imod+1)).mkdir(parents=True) 
            # we are dealing with the imod-th mode 
            x = x_fac / np.sqrt(self.masses[imod] * self.freqs[imod])
            norm = np.sqrt( ( self.vecs[imod, :, :] ** 2).sum() ) # in case the vectors are not strictly normalized 
            for ipt in range(npts):
                # the ipt-th displacement
                prefix = "./perb/" + str(imod+1) + "/" + str(ipt+1) 
                Path(prefix).mkdir(parents=True)  
                new_crd = self.crds + self.vecs[imod, :, :]/norm  * x * scalings[ipt]
                outfile = prefix + "/" + "mol_{x}_{y}.xyz".format(x=imod+1,y=ipt+1)
                with open(outfile, "w") as f:
                    # print out coordinates 
                    print(self.orig.get_number_of_atoms(), file=f)
                    print("ModeNo. %5d  StructureNo. %4d scaling= %4.1f%%" % (imod+1, ipt + 1, scalings[ipt] * 100), file=f) 
                    for i in range(self.orig.get_number_of_atoms()):
                        print(" %6s %12.6f %12.6f %12.6f" % (symbols[self.orig.get_atomic_numbers()[i]], new_crd[i,0], new_crd[i,1], new_crd[i, 2]), file=f)
        # "echo out some information..."
        print("All %5d modes are perturbed, the results stored in ./perb/" % self.freqs.size) 
        print("Structures are generated at ", end=" ")
        percents = ["%4.1f%%  " % x for x in scalings * 100] 
        print(" ".join(percents), " of the amplitude of the zero-point motion.")  
        
# the actual work is done here. 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
        A script to generate structures perturbed along normal modes 
        provided by a freq calculation done in Gaussian.""")
    parser.add_argument("--out", dest="outfile", required=True, help="The filename of vasp outcar.")
    parser.add_argument("--pos", dest="posfile", required=True, help="The filename of poscar.")
    parser.add_argument("--npts", "-np", dest="npts", type=int, default=10, help="Number of points along each modes. Default is 10.")
    parser.add_argument("--maxdp", "-mdp", dest="mdp", type=float, default=1, help="The max displacement of each modes. Default is 1.")
    parser.add_argument("--verbose", "-v", action="store_true", help="print out more information of the structures for debugging purpose.")
    args = parser.parse_args()
    
    mol = mol_freq(args.outfile,args.posfile) 
    mol.read_geometry()
    mol.read_modes()
    if args.verbose:
        print("atom numbers:\n", mol.orig.get_number_of_atoms(), "\n")
        print("atom equilibrium coordinates (A):\n", mol.crds, "\n")  
        print("freqs (cm-1):\n", mol.freqs, "\n")
        print("displacement vecs:\n", mol.vecs, "\n") 
        # check if the displacemnet vectors are normalized 
        for i in range(mol.freqs.size):
            print("The mode No. %4d vector norm is: %10.5f" % ((i+1), np.sqrt(( mol.vecs[i, :, :] * mol.vecs[i, :, :] ).sum())) )
    # perturbation starts here.
    mol.perturb(args.npts,args.mdp)  
