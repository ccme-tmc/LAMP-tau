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
    
def extract_Dtensor(filename):
    """
    extract  Dtensor from the file filename, return three list of B_kq, correspond to k=2,4,6)  
    """
    pattern = re.compile(r"[+-]?\d+\.\d+|[+-]?\d+")   #determine the re pattern of number needed to be found
    Dm = np.empty([3,3], dtype = float)   #D-matrix
    Dt = np.empty([3], dtype = float)   #D-matrix
    Dvec = np.empty([3,3], dtype = float)   #D-matrix
    with open(filename,"r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            if len(line.split()) == 7 and line.split()[0] == "D-matrix": #find the correct position
                for i in range(0,3):
                    line = f.readline()
                    #print(line)
                    for j in range(0,3):
                        Dm[i,j] = float(pattern.findall(line.split("j")[j])[0])
                        #print(pattern.findall(line.split("j")[0])[0])
                line = f.readline()
                line = f.readline()
                for i in range(0,3):
                    Dt[i] = float(pattern.findall(line.split()[i])[0])
                #print(Dt)
                line = f.readline()
                for i in range(0,3):
                    line = f.readline()
                    #print(line)
                    for j in range(0,3):
                        Dvec[i,j] = float(pattern.findall(line.split("j")[j])[0])
                        #print(pattern.findall(line.split("j")[0])[0])
    return Dm,Dt,Dvec
    
def extract_Energy(filename):
    """
    extract  Dtensor from the file filename, return three list of B_kq, correspond to k=2,4,6)  
    """
    pattern = re.compile(r"[+-]?\d+\.\d+|[+-]?\d+")   #determine the re pattern of number needed to be found
    E = np.zeros([5], dtype = float)   #D-matrix
    if_find = True
    with open(filename,"r") as f:
        while True:
            line = f.readline()
            if not line:
                print("no find")
                if_find = False
                return E,if_find
                break
            if len(line.split()) == 4 and line.split()[0] == "Selected": #find the correct position
                line = f.readline()
                #print(line.split())
                if len(line.split()) == 6:
                    E[0] = float(pattern.findall(line.split()[1])[0])
                    E[1] = float(pattern.findall(line.split()[3])[0])
                    E[2] = float(pattern.findall(line.split()[5])[0])
                    line = f.readline()
                    #print(line.split())
                    E[3] = float(pattern.findall(line.split()[1])[0])
                    E[4] = float(pattern.findall(line.split()[3])[0])
                    return E,if_find
                
                if len(line.split()) == 5:
                    E[0] = float(pattern.findall(line.split()[1])[0])
                    E[1] = float(pattern.findall(line.split()[3])[0])
                    line = f.readline()
                    #print(line.split())
                    E[2] = float(pattern.findall(line.split()[0])[0])
                    E[3] = float(pattern.findall(line.split()[2])[0])
                    line = f.readline()
                    E[4] = float(pattern.findall(line.split()[0])[0])
                    return E,if_find


                    
                if len(line.split()) == 4:
                  
                    E[0] = float(pattern.findall(line.split()[0])[0])*10**float(pattern.findall(line.split()[0])[1])
                    E[1] = float(pattern.findall(line.split()[2])[0])*10**float(pattern.findall(line.split()[2])[1])
                    line = f.readline()
                    E[2] = float(pattern.findall(line.split()[0])[0])*10**float(pattern.findall(line.split()[0])[1])
                    E[3] = float(pattern.findall(line.split()[2])[0])*10**float(pattern.findall(line.split()[2])[1])
                    line = f.readline()
                    E[4] = float(pattern.findall(line.split()[0])[0])*10**float(pattern.findall(line.split()[0])[1])
                    return E,if_find


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
    parser.add_argument("--CFPs", "-c", action="store_true", help="extract the CFPs from Molcas")
    parser.add_argument("--Dtensor", "-d", action="store_true", help="extract the Dtensor from Pyscf")
    parser.add_argument("--Energy", "-e", action="store_true", help="extract the Dtensor from Pyscf")
    args = parser.parse_args()

    if args.Dtensor:
        Dm0,Dt0,Dvec0=extract_Dtensor(args.eqlog)   #read the eq CFPs
        with open("Dtensor.out","w") as f:
            f.write("The Perturbation Dtensors\n")
            f.write("Total modes {x}\n".format(x=args.nmds))
            f.write("Total points {y}\n".format(y=args.npts))
            for imod in range(0,args.nmds):
                for ipt in range(0,args.npts):
                    print("Working on the Mode {x} Point {y}".format(x=imod+1,y=ipt+1))
                    f.write("Mode {x} Point {y}\n".format(x=imod+1,y=ipt+1))
                    mag_filename = "perb/{x}/{y}/stdout2".format(name=args.name,x=imod+1,y=ipt+1)
                    try:
                        f2 =open(mag_filename,"r")
                        f2.close()
                    except FileNotFoundError:
                        f.write("None\n")
                        print("no file "+mag_filename)
                        continue
                    Dm1,Dt1,Dvec1=extract_Dtensor(mag_filename) #extract each CFPs
                    Dt = Dt1-Dt0
                    f.write(' '.join(map(str,Dt))+"\n")
                    
        with open("Dmatrix.out","w") as f:
            f.write("The Perturbation Dtensors\n")
            f.write("Total modes {x}\n".format(x=args.nmds))
            f.write("Total points {y}\n".format(y=args.npts))
            for imod in range(0,args.nmds):
                for ipt in range(0,args.npts):
                    print("Working on the Mode {x} Point {y}".format(x=imod+1,y=ipt+1))
                    f.write("Mode {x} Point {y}\n".format(x=imod+1,y=ipt+1))
                    mag_filename = "perb/{x}/{y}/stdout2".format(name=args.name,x=imod+1,y=ipt+1)
                    try:
                        f2 =open(mag_filename,"r")
                        f2.close()
                    except FileNotFoundError:
                        f.write("None\n")
                        print("no file "+mag_filename)
                        continue
                    Dm1,Dt1,Dvec1=extract_Dtensor(mag_filename) #extract each CFPs
                    Dm = Dm1-Dm0
                    f.write(' '.join(map(str,Dm[0,...]))+"\n")
                    f.write(' '.join(map(str,Dm[1,...]))+"\n")
                    f.write(' '.join(map(str,Dm[2,...]))+"\n")
                    
    if args.Energy:
        E0,iftest=extract_Energy(args.eqlog)   #read the eq CFPs
        print(E0)
        with open("Energy.out","w") as f:
            f.write("The Perturbation Energys\n")
            f.write("Total modes {x}\n".format(x=args.nmds))
            f.write("Total points {y}\n".format(y=args.npts))
            for imod in range(0,args.nmds):
                for ipt in range(0,args.npts):
                    print("Working on the Mode {x} Point {y}".format(x=imod+1,y=ipt+1))
                    f.write("Mode {x} Point {y}\n".format(x=imod+1,y=ipt+1))
                    mag_filename = "perb/{x}/{y}/stdout2".format(name=args.name,x=imod+1,y=ipt+1)
                    try:
                        f2 =open(mag_filename,"r")
                        f2.close()
                    except FileNotFoundError:
                        f.write("None\n")
                        print("no file "+mag_filename)
                        continue
                    E1,iftest=extract_Energy(mag_filename) #extract each CFPs
                    if not iftest:
                        continue
                    E = E1-E0
                    f.write(' '.join(map(str,E))+"\n")
                    
    if args.CFPs:
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
                        f.write("None\n")
                        print("no file "+mag_filename)
                        continue
                    B=extract_CFPs(mag_filename) #extract each CFPs
                    #B = B-B0
                    f.write(' '.join(map(str,B))+"\n")
            
    






