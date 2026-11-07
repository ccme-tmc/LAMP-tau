#!/usr/bin/env python3
#/usr/bin/python3
# -*- coding:utf-8 -*-

import argparse
import os 

def fileInsertion(inFile, outFile, repList):
    """
        purpose:
            To create molcas input files, we need to insert files into an existing file at some 
            predefined places marked by strings.
        parameters:
            inFile: name of the input file to be changed.
            outFile: name of the output file
            repList: list of tuples
                e.g. [("str1", "file1.dat"), ("str2", "file2.dat")]
                replaces "str1" line and "str2" line in inFile by the content of file1.dat and file2.dat
                respectively.
            Note:
                match string must occupy a separate line
    """
    # read data
    with open(inFile, "r") as f:
        inpList = f.readlines()
    # read replacement information
    nrep = len(repList)
    strList = []
    fileList = []
    for i in range(nrep):
        strList.append(repList[i][0]) 
        fileList.append(repList[i][1])
    with open(outFile, "w") as f:
        for i in range(len(inpList)):
            if inpList[i].strip() in strList:
                fileIndex = strList.index(inpList[i].strip())
                f.writelines(fileList[fileIndex])
            else:
                f.write(inpList[i])
    return

parser = argparse.ArgumentParser(description="""
    To replace place holders in the template molcas input file with real coordinates.""")
parser.add_argument("--tmp", help="The template Molcas input file. defaulting to mag_tmp.in", default="mag_tmp.in")
parser.add_argument("--xyz",  help="The structure file in xyz format. defaulting to mol.xyz", default="mol.xyz")
parser.add_argument("--out", "-o", help="The final input file name, defaulting to mag.in", default="mag.in")
parser.add_argument("--charge", "-c", help="The Charge of the system", default=0)
parser.add_argument("--title", "-t", help="The Title of Calculation", default=0)
args = parser.parse_args()

repname1 = 'm_name = "{x}"\n'.format(x = args.xyz)
repname2 = 'm_ch = {x}\n'.format(x = args.charge)
repname3 = 'title = "{x}"\n'.format(x = args.title)

rep_list = [("Molecular_Name",repname1),("Molecular_Charge",repname2),("Calculate_Title",repname3)]

fileInsertion(args.tmp, args.out, rep_list)

