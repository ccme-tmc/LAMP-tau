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
from gen_structures_cell import mol_freq
#from Analyze_CFPs import Eigenstates
#from matplotlib import pyplot as plt



def extract_eq_CFPs(filename):
    """
    extract  CFPs from the file filename, return three list of B_kq, correspond to k=2,4,6)  
    """
    pattern = re.compile(r"[+-]?\d+\.\d+|[+-]?\d+")   #determine the re pattern of number needed to be found
    B2 = np.empty([5], dtype = float)   #B_kq,where k = 2
    B4 = np.empty([9], dtype = float)   #B_kq,where k = 4
    B6 = np.empty([13], dtype = float)  #B_kq,where k = 6
    with open(filename,"r", encoding='utf-8') as f:
        while True:
            line = f.readline()
            if not line:
                break
            if len(line.split()) == 4 and " ".join(line.split()) == "Rotated Crystal Field Parameters:": #find the correct position
                for i in range(0,5):
                    line = f.readline()
                for i in range(0,5):
                    line = f.readline()
                    CFPs = line.split('|')[3]
                    B2[i] = float(CFPs)

                for i in range(0,9):
                    line = f.readline()
                    CFPs = line.split('|')[3]
                    B4[i] = float(CFPs)

                for i in range(0,13):
                    line = f.readline()
                    CFPs = line.split('|')[3]
                    B6[i] = float(CFPs)

                break
    return B2,B4,B6

def extract_CFPs(filename,A,B,C):
    """
    extract  CFPs from the file filename, return three list of B_kq, correspond to k=2,4,6)  
    """
    pattern = re.compile(r"[+-]?\d+\.\d+|[+-]?\d+")   #determine the re pattern of number needed to be found
    B2 =np.array(A, dtype = float)  #B_kq,where k = 2
    B4 =np.array(B, dtype = float)  #B_kq,where k = 4
    B6 =np.array(C, dtype = float)  #B_kq,where k = 6
    with open(filename,"r", encoding='utf-8') as f:
        while True:
            line = f.readline()
            if not line:
                break
            if len(line.split()) == 4 and " ".join(line.split()) == "Rotated Crystal Field Parameters:": #find the correct position
                for i in range(0,5):
                    line = f.readline()
                for i in range(0,5):
                    line = f.readline()
                    CFPs = line.split('|')[3]
                    B2[i] = float(CFPs)

                for i in range(0,9):
                    line = f.readline()
                    CFPs = line.split('|')[3]
                    B4[i] = float(CFPs)

                for i in range(0,13):
                    line = f.readline()
                    CFPs = line.split('|')[3]
                    B6[i] = float(CFPs)

                break
    return B2,B4,B6

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
        Extract the perturbation CFPs from the per-displacement calculations.
        The per-point output is read from <prefix>/<mode>/<point>/<stdout>.""")
    parser.add_argument("--eqlog", "-el", dest="eqlog", required=True, help="The path of your eq structure mag file.")
    parser.add_argument("--nmds", "-nd", dest="nmds", type=int, required=True, help="Number of vibration modes.")
    parser.add_argument("--npts", "-np", dest="npts", type=int, default=10, help="Number of points along each modes. Default is 10.")
    parser.add_argument("--prefix", dest="prefix", default="perb", help="Directory containing the per-displacement calculations (default: perb).")
    parser.add_argument("--stdout", dest="stdout", default="stdout1", help="Name of the per-point output file (default: stdout1).")
    parser.add_argument("--out", dest="out", default="CFPs.out", help="Output CFP file (default: CFPs.out).")
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

    B2_0,B4_0,B6_0=extract_eq_CFPs(args.eqlog)   #read the eq CFPs
    #print(B2_0,B4_0,B6_0)
    B2t = np.zeros([5], dtype = float)   #B_kq,where k = 2
    B4t = np.zeros([9], dtype = float)   #B_kq,where k = 4
    B6t = np.zeros([13], dtype = float)  #B_kq,where k = 6
    with open(args.out,"w") as f:
        f.write("The Perturbation CFPs\n")
        f.write("\n")
        f.write("Total modes {x}\n".format(x=args.nmds))
        f.write("Total points {y}\n".format(y=args.npts))
        for imod in range(0,args.nmds):
            for ipt in range(0,args.npts):
                print("Working on the Mode {x} Point {y}".format(x=imod+1,y=ipt+1))
                f.write("Mode {x} Point {y}\n".format(x=imod+1,y=ipt+1))
                mag_filename = os.path.join(
                    args.prefix, str(imod + 1), str(ipt + 1), args.stdout
                )
                try:
                    f2 =open(mag_filename,"r")
                    f2.close()
                except FileNotFoundError:
                    print("no file "+mag_filename)
                    f.write(" ".join(str(s) for s in B2t)+"\n")
                    f.write(" ".join(str(s) for s in B4t)+"\n")
                    f.write(" ".join(str(s) for s in B6t)+"\n")
                    continue
                B2,B4,B6=extract_CFPs(mag_filename,B2_0,B4_0,B6_0) #extract each CFPs
                B2 = B2-B2_0
                B4 = B4-B4_0
                B6 = B6-B6_0                 #get the perturbation CFPs
                f.write(" ".join(str(s) for s in B2)+"\n")
                f.write(" ".join(str(s) for s in B4)+"\n")
                f.write(" ".join(str(s) for s in B6)+"\n")
            
    






