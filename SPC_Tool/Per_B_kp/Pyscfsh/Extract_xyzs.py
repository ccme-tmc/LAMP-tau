#!/usr/bin/env python3
# -*- coding:utf-8 -*-

import numpy as np
from numpy import polyfit, poly1d
import argparse
from pathlib import Path
import shutil 
import os 
import sys
import re
#from Analyze_CFPs import Eigenstates
#from matplotlib import pyplot as plt



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

def extract_CFPs(filename):
    """
    extract  CFPs from the file filename, return three list of B_kq, correspond to k=2,4,6)  
    """
    pattern = re.compile(r"[+-]?\d+\.\d+|[+-]?\d+")   #determine the re pattern of number needed to be found
    with open(filename,"r") as f:
        line = f.readline()
        line = f.readline()
        line = f.readline()
        B=float(line)
    return B

def get_Modes_list(nmds,Eigen,mol,select_sigma):
    Modes_list = np.empty((0,),dtype = int)
    for j in range(0,nmds):
        for f in range(0,Eigen.n_states-2,2):
            for i in range(f+2,Eigen.n_states,2):
                omega = abs(Eigen.Vaules[f]-Eigen.Vaules[i])
                if abs(mol.freqs[j]-omega) < select_sigma:
                    Modes_list = np.append(Modes_list,j)            
    Modes_list = np.unique(Modes_list)
    return Modes_list


# the actual work is done here. 
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="""
        A script to extract the perturbation CFPs from a single_aniso calculation done in MOLCAS.
        The path of outfile should be "perb/x/y/mag_x_y.log", 
        and the path of eq outfile should be "name/eq.out" where name is the main fold path you should input; 
        "x" means the number of vb mode, from 1 to max;"y" means the number of points of each mode, from 1 to max. """)
    parser.add_argument("--name", dest="name", required=True, help="The name of your system.")
    parser.add_argument("--eqlog", "-el", dest="eqlog", required=True, help="The path of your eq structure mag file.")
    parser.add_argument("--nmds", "-nd", dest="nmds", type=int, required=True, help="Number of vibration modes.")
    parser.add_argument("--npts", "-np", dest="npts", type=int, default=10, help="Number of points along each modes. Default is 10.")
    parser.add_argument("--verbose", "-v", action="store_true", help="print out more information of the structures for debugging purpose.")
    args = parser.parse_args()
    
    """
    parser = argparse.ArgumentParser(description=""
        A script to extract the perturbation CFPs from a single_aniso calculation done in MOLCAS.."")
    parser.add_argument("--yml", dest="yml", required=True, help="The path of your .yml parameter list input files.")
    args = parser.parse_args()
    """
    """
    Input_filename = "Modes.yml"
    with open(Input_filename,"r") as file:
        finput = file.read()
    pinput = yaml.safe_load(finput)
    
    CFPs_foldname = pinput['CFPs_foldname']
    mag_filename = pinput['mag_filename']
    freq_filename = pinput['freq_filename']
    nmds = pinput['nmds']
    npts = pinput['npts']
    Correction_slope = float(pinput['Correction_slope'])    #Correction_slope of freq
    Corrected_intercept = float(pinput['Corrected_intercept'])    #Corrected_intercept of freq
    If_select_sigma = bool(pinput['If_select_sigma'])
    """
    """
    if If_select_sigma:
        print("your modelist depends on select_sigma")
        select_sigma = pinput['select_sigma']
        mol = mol_freq(freq_filename) 
        mol.read_geometry()
        mol.read_modes()
        mol.freqs = Correction_slope*mol.freqs+Corrected_intercept
        Eigen = Eigenstates(mag_filename)
        Eigen.read_energies()
        Eigen.read_wavefunction()
        Modes_list = get_Modes_list(nmds,Eigen,mol,select_sigma)
    else:
        print("your modelist is all")
        Modes_list = np.linspace(0,nmds-1,num = nmds,dtype = int)
    
    print(Modes_list)
    """

    B0=extract_CFPs(args.eqlog)   #read the eq CFPs
    with open("CFPs.out","w") as f:
        f.write("The Perturbation CFPs\n")
        f.write("Total modes {x}\n".format(x=args.nmds))
        f.write("Total points {y}\n".format(y=args.npts))
        for imod in range(0,args.nmds):
            for ipt in range(0,args.npts):
                print("Working on the Mode {x} Point {y}".format(x=imod+1,y=ipt+1))
                f.write("Mode {x} Point {y}\n".format(x=imod+1,y=ipt+1))
                mag_filename = "perb/{x}/{y}/{name}_mag.txt".format(name=args.name,x=imod+1,y=ipt+1)
                try:
                    f2 =open(mag_filename,"r")
                    f2.close()
                except FileNotFoundError:
                    print("no file "+mag_filename)
                    continue
                B=extract_CFPs(mag_filename) #extract each CFPs
                #B = B-B0
                f.write(str(B)+"\n")
            
    






