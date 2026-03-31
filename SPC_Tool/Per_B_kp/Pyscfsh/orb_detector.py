#!/home/dengyf/.conda/envs/my_env/bin/python
#/usr/bin/python3
# -*- coding:utf-8 -*-

import numpy as np
import re
import argparse

def read_molden(fname):
    orbs = []
    with open(fname,'r') as sav:
        while True:
            line = sav.readline()
            if re.search(r'Sym',line):
                sav.readline()
                sav.readline()
                sav.readline()

                orb_coeff = []
                for i in range(101):
                    coeff = sav.readline().split()
                    orb_coeff.append(float(coeff[1]))
                orbs.append(orb_coeff)
            if not line:
                break

    return np.asarray(orbs)

def Sort_EV(E,V):
    """
    sort the 1-d rate_list E, and get the corresponding eigenvector V
    """
    Elength = np.shape(E)[0]
    while Elength > 0:
        for i in range(0,Elength-1):
            if E[i] < E[i+1]:
                E[i], E[i+1] = E[i+1], E[i]
                V[i], V[i+1] = V[i+1], V[i]
        Elength -= 1
    return E,V

parser = argparse.ArgumentParser(description="""
    To replace place holders in the template molcas input file with real coordinates.""")
parser.add_argument("--title","-t", help="The title", default="DMET")
parser.add_argument("--act","-a",  help="The active electrons numbers", default="0")
parser.add_argument("--orbs", "-o", help="The active orbitals", default="0")
#parser.add_argument("--charge", "-c", help="The Charge of the system", default=0)
args = parser.parse_args()

n_orbital = args.orbs
n_acte = args.act
title = args.title
molden = "{x}_imp_rohf_orbs.molden".format(x=title)
outfile = "{x}_cas_info".format(x=title)
orbs = read_molden(molden)

ans = []
ansnumber = []
for i in range(len(orbs)):
    if np.sum(orbs[i][18:18+5]**2) > 0.1 and np.sum(orbs[i][18:18+5]**2) < 0.8:
        ans.append(str(i+1))
        ansnumber.append(np.sum(orbs[i][18:18+5]**2))
if int(np.shape(ans)[0]) < 5:
    print("The number of ans < n_orbital, please check it")
elif int(np.shape(ans)[0]) >= 5:
    Sort_EV(ansnumber,ans)
    new = ans[0:5]
    new = np.sort(new)
    print(' '.join(new))
    if int(new[4]) < 180:
        with open(outfile,"w") as file:    
            file.write("{x} {y}\n".format(x=n_orbital,y=n_acte))
            file.write(' '.join(new))
            file.write("\n")
