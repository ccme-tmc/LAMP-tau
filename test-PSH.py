# %%
import numpy as np
import pandas as pd
import math
from scipy.sparse.linalg import eigs, eigsh
from scipy.linalg.lapack import zgels, zgelss
from prettytable import PrettyTable

from sympy.physics.quantum.cg import CG
from sympy import S

pd.set_option('display.width',1000)

def setevens_mat_old(k,q,j1,m1,j2,m2):
    if q==0:
        return spherical_mat(k,q,j1,m1,j2,m2)
    if q>0:
        mat1 = spherical_mat(k,q,j1,m1,j2,m2)
        mat2 = spherical_mat(k,-1*q,j1,m1,j2,m2)
        matelem = (-1**q*mat1+mat2)/np.sqrt(2)
        return matelem
    if q<0:
        mat1 = spherical_mat(k,q,j1,m1,j2,m2)
        mat2 = spherical_mat(k,-1*q,j1,m1,j2,m2)
        matelem = (-1**(q+1)*mat1+mat2)/np.sqrt(2)
        matelem = -1*matelem.imag+1j*matelem.real
        return matelem
    
def spherical_mat(k,q,j1,m1,j2,m2):
    cg = CG( S(2*j2)/2, S(2*m2)/2, S(2*k)/2, S(2*q)/2, S(2*j1)/2, S(2*m1)/2, ).doit()
    redelem = spintens(k,j1)
    return complex(cg)*redelem/(np.sqrt(2*j1+1))

def sto_redelem(k,j):
    matelem=math.factorial(int(2*j+k+1))
    matelem=matelem/math.factorial(2*k)
    matelem=matelem/math.factorial(int(2*j-k))
    matelem=np.sqrt(matelem)
    matelem=(-1)**(k)*math.factorial(k)
    matelem=matelem*(-1)**(k)*np.sqrt(math.factorial(2*k))
    matelem=matelem/(2)**(k)
    matelem=matelem/math.factorial(k)
    return matelem
    
def setevens_mat(k,q,j1,m1,j2,m2):
    """
    generate the element <j1,m1|T_kq|j2,m2> of Tesseral Tensor Operators(TTO) matrix 
    
    Args:

    Returns:

    """
    if q==0:
        cg = CG( S(2*j2)/2, S(2*m2)/2, S(2*k)/2, S(2*q)/2, S(2*j1)/2, S(2*m1)/2).doit()
        redelem = spintens(k,j1)
        return complex(cg)*redelem/(np.sqrt(2*j1+1))
    if q>0:
        if abs(int(q/2)-q/2) < 0.1:
            cg1 = CG( S(2*j2)/2, S(2*m2)/2, S(2*k)/2, S(2*q)/2, S(2*j1)/2, S(2*m1)/2).doit()
            cg2 = CG( S(2*j2)/2, S(2*m2)/2, S(2*k)/2, S(-2*q)/2, S(2*j1)/2, S(2*m1)/2).doit()
            cg = complex(cg1)+complex(cg2)
            redelem = spintens(k,j1)
            return complex(cg)*redelem/(np.sqrt(2*j1+1)*np.sqrt(2))
        else:
            cg1 = CG( S(2*j2)/2, S(2*m2)/2, S(2*k)/2, S(2*q)/2, S(2*j1)/2, S(2*m1)/2).doit()
            cg2 = CG( S(2*j2)/2, S(2*m2)/2, S(2*k)/2, S(-2*q)/2, S(2*j1)/2, S(2*m1)/2).doit()
            cg = -complex(cg1)+complex(cg2)
            redelem = spintens(k,j1)
            return complex(cg)*redelem/(np.sqrt(2*j1+1)*np.sqrt(2))
    if q<0:
        if abs(int(q/2)-q/2) < 0.1:
            cg1 = CG( S(2*j2)/2, S(2*m2)/2, S(2*k)/2, S(2*q)/2, S(2*j1)/2, S(2*m1)/2).doit()
            cg2 = CG( S(2*j2)/2, S(2*m2)/2, S(2*k)/2, S(-2*q)/2, S(2*j1)/2, S(2*m1)/2).doit()
            cg = (complex(cg1)-complex(cg2))*1j
            redelem = spintens(k,j1)
            return complex(cg)*redelem/(np.sqrt(2*j1+1)*np.sqrt(2))
        else:
            cg1 = CG( S(2*j2)/2, S(2*m2)/2, S(2*k)/2, S(2*q)/2, S(2*j1)/2, S(2*m1)/2).doit()
            cg2 = CG( S(2*j2)/2, S(2*m2)/2, S(2*k)/2, S(-2*q)/2, S(2*j1)/2, S(2*m1)/2).doit()
            cg = (complex(cg1)+complex(cg2))*1j
            redelem = spintens(k,j1)
            return complex(cg)*redelem/(np.sqrt(2*j1+1)*np.sqrt(2))
    
def spintens(K,J):
    """
    calculate the value of <J||T_kq||J>
    
    Args:

    Returns:

    """
    B=1
    for i in range(-K,K+1):
        B=B*(2*J+1+i)
    B=B/(math.factorial(2*K))

    B=np.sqrt(B)
    B=B*math.factorial(K)
    ak = np.sqrt(2**(-K))

    #ak = np.sqrt(math.factorial(2*K))/math.factorial(2*K)

    #ak = math.factorial(2*K)*math.factorial(int(2*S-K))/math.factorial(int(2*S+K+1))
    #ak = np.sqrt(ak)/math.factorial(K)

    # ak = math.factorial(int(S+K/2))/math.factorial(int(S-K/2))/math.factorial(int(K/2))/math.factorial(int(K/2))
    # ak = ak*math.factorial(int(2*S-K))/math.factorial(int(2*S+K+1))
    # ak = ak*(2*S+1)*np.sqrt(math.factorial(2*K))
    return B*ak

def Out_Orca(Ma,ncol,filename):
    """
    Output the data of Matrix Ma as the format of Orca-Molforge.

    Args:
        Ma (np): the data you want to output.
        ncol (int): the number of the column in the output file, in Orca is 6.
        filename (str): the name of the output filename.

    Returns:
        bool (bool): success finished.
    """
    tcol = Ma.shape[1]
    ntab = tcol//ncol
    nlast = tcol%ncol
    with open(filename,"w") as file:
        print("part", file=file)
        for ta in range(ntab):
            cols = 0+6*ta
            cole = 6+6*ta
            Matemp = Ma[...,cols:cole]
            df = pd.DataFrame(Matemp,columns=range(cols,cole))
            print(df.applymap(lambda x: '%.6f'%x), file=file)
        cols = tcol-nlast
        cole = tcol
        Matemp = Ma[...,cols:cole]
        df = pd.DataFrame(Matemp,columns=range(cols,cole))
        print(df.applymap(lambda x: '%.6f'%x), file=file)
    return True

def read_Orca(Mosize,ncol,filename):
    """
    read the data from Orca-Molforge.

    Args:
        Mosize (int): the size of the molecule orbitals.
        ncol (int): the number of the column in the input file, in Orca is 6.
        filename (str): the name of the input filename.

    Returns:
        Ma (np): the data you read.
    """
    Ma = np.zeros([Mosize,Mosize],dtype=float)
    nrow = Mosize
    tcol = Mosize
    ntab = tcol//ncol
    nlast = tcol%ncol
    with open(filename,"r") as file:
        line = file.readline()
        for ta in range(ntab):
            line = file.readline()
            cols = 0+6*ta
            cole = 6+6*ta
            for nl in range(nrow):
                line = file.readline()
                #print(line.split())
                Ma[nl,cols:cole] = line.split()[1:]
        cols = tcol-nlast
        cole = tcol
        line = file.readline()
        for nl in range(nrow):
                line = file.readline()
                #print(line.split())
                Ma[nl,cols:cole] = line.split()[1:]
    return Ma
 

def gen_ESO_matrix_molfrorge(ESO_foldname,Nstate,lmax):
    """
    generate the Tesseral Tensor Operators(TTO) matrix like Molforge
    
    Args:
        ESO_foldname (str): the fold to save the matrix of Tensor Operators.
        Nstate (int): the number of the state of |JM>.
        lmax (int): the highest order of the Tensor Operators.

    Returns:
        ESO (np): the matrix of Tensor Operators.
    """
    Odim = int(lmax*(lmax+3)/2)
    ESO = np.empty([Nstate,Nstate,Odim], dtype = complex)
    J = float((Nstate-1)/2)
    filenamer = ESO_foldname+"/"+str(Nstate)+"_real.eso"
    filenamei = ESO_foldname+"/"+str(Nstate)+"_imag.eso"
    with open(filenamer,"w") as filer:
        with open(filenamei,"w") as filei:
            i = 0
            for k in range(2,lmax+2,2):
                for q in range(-k,k+1):
                    print("k={x},q={y}".format(x=k,y=q),file=filer)
                    print("k={x},q={y}".format(x=k,y=q),file=filei)
                    for m1 in range(Nstate):
                        Mj1 = m1-J
                        for m2 in range(Nstate):
                            Mj2 = m2-J
                            ESO[m1,m2,i] = setevens_mat(k,q,J,Mj1,J,Mj2)
                    for j in range(0,Nstate):
                        filer.write(" ".join(str(s) for s in ESO.real[j,...,i])+"\n")
                        filei.write(" ".join(str(s) for s in ESO.imag[j,...,i])+"\n")
                    i = i+1
    return ESO

def gen_ESO_matrix_easyspin(ESO_foldname,Nstate,lmax):
    """
    generate the extend Steven's operators(ESO) matrix from easyspin in matlab
        
    Args:
        ESO_foldname (str): the fold to save the matrix of Steven's operators.
        Nstate (int): the number of the state of |JM>.
        lmax (int): the highest order of the Steven's operators.

    Returns:
        ESO2 (np): the matrix of Steven's operators.
    """
    Odim = int(lmax*(lmax+3)/2)
    import matlab.engine
    ESO = np.empty([Nstate,Nstate,Odim], dtype = complex)
    ESO2 = np.empty([Nstate,Nstate,Odim], dtype = complex)
    kq = np.empty([2], dtype = float)
    J = float((Nstate-1)/2)
    eng = matlab.engine.start_matlab()
    i = 0
    for k in range(2,lmax+2,2):
        for q in range(-k,k+1):
            kq[0] = k
            kq[1] = q
            ESO[...,i] = eng.stev(float(J),kq)
            i = i+1
    eng.quit()
    for i in range(Nstate):
        for j in range(Nstate):
            ESO2[i,j] = ESO[Nstate-i-1,Nstate-j-1]
    filenamer = ESO_foldname+"/"+str(Nstate)+"_real.eso"
    filenamei = ESO_foldname+"/"+str(Nstate)+"_imag.eso"
    with open(filenamer,"w") as filer:
        with open(filenamei,"w") as filei:
            i = 0
            for k in range(2,lmax+2,2):
                for q in range(-k,k+1):
                    print("k={x},q={y}".format(x=k,y=q),file=filer)
                    print("k={x},q={y}".format(x=k,y=q),file=filei)
                    for j in range(0,Nstate):
                        filer.write(" ".join(str(s) for s in ESO2.real[j,...,i])+"\n")
                        filei.write(" ".join(str(s) for s in ESO2.imag[j,...,i])+"\n")
                    i = i+1
    eng.quit()
    for i in range(Nstate):
        for j in range(Nstate):
            ESO2[i,j] = ESO[Nstate-i-1,Nstate-j-1]
    return ESO2

def read_ESO_matrix(ESO_foldname,Nstate,lmax):
    """
    read the Tensor operators matrix from the fold 
        
    Args:
        ESO_foldname (str): the fold to read the matrix of Tensor operators.
        Nstate (int): the number of the state of |JM>.
        lmax (int): the highest order of the Steven's operators.

    Returns:
        ESO (np): the matrix of Tensor operators.
    """
    Odim = int(lmax*(lmax+3)/2)
    ESO = np.empty([Nstate,Nstate,Odim], dtype = complex)
    A = np.empty([Nstate,Nstate,Odim], dtype = float)
    B = np.empty([Nstate,Nstate,Odim], dtype = float)
    filename = ESO_foldname+"/"+str(Nstate)+"_real.eso"
    with open(filename,"r") as file:
        for i in range(0,Odim):
            line = file.readline()
            for j in range(0,Nstate):
                line = file.readline()
                A[j,...,i] = line.split()
    filename = ESO_foldname+"/"+str(Nstate)+"_imag.eso"
    with open(filename,"r") as file:
        for i in range(0,Odim):
            line = file.readline()
            for j in range(0,Nstate):
                line = file.readline()
                B[j,...,i] = line.split()
    ESO = A+B*1j
    return ESO

def get_matrix_from_Orca(Mosize,ncol,foldname):
    """
    read matrix from Molforge output(.dat) from Orca

    Args:
        Mosize (int): the size of the molecule orbitals.
        ncol (int): the number of the column in the input file, in Orca is 6.
        foldname (str): the fold where the output is.

    Returns:
        Masoc (np): the matrix of Hamiltonian.
        Lx (np): the matrix of Lx operator.
        Lz (np): the matrix of Lz operator.
        Sx (np): the matrix of Sx operator.
        Sz (np): the matrix of Sz operator.
    """
    Masoc = read_Orca(Mosize,ncol,foldname+"SOCR.dat")+read_Orca(Mosize,ncol,foldname+"SOCI.dat")*1j
    Lx = read_Orca(Mosize,ncol,foldname+"Lx.dat")*1j
    Lz = read_Orca(Mosize,ncol,foldname+"Lz.dat")*1j
    Sx = read_Orca(Mosize,ncol,foldname+"Sx.dat")*0.5
    Sz = read_Orca(Mosize,ncol,foldname+"Sz.dat")*0.5
    return Masoc,Lx,Lz,Sx,Sz 

def get_E_J(Masoc,Lx,Lz,Sx,Sz,Nstate):
    """
    calculate the eigenvalue/vecor of the Hamiltonian and the Jz operator
    
    Args:
        Masoc (np): the matrix of Hamiltonian.
        Lx (np): the matrix of Lx operator.
        Lz (np): the matrix of Lz operator.
        Sx (np): the matrix of Sx operator.
        Sz (np): the matrix of Sz operator.
        Nstate (int): the number of the state of |JM>.

    Returns:
        Jvec (np): the eigenvector of Jz operator.
        enr2 (np): the energy.
    """
    eigval, eigvec = np.linalg.eigh(Masoc)
    enr2 = np.zeros([Nstate],dtype=complex)
    print("#################################################")
    print("#################################################")
    print("Nj eigevanlues of the electronic H")
    print("#################################################")
    print("#################################################")
    table = PrettyTable(["EIGstate","Energy(cm^{-1})"])
    for i in range(0,Nstate):
        enr2[i] = (eigval[i]-eigval[0])*219474.63
        table.add_row([i+1,enr2[i].real])
    print(table)
        
    #print(pd.DataFrame(Lz))
    #print(eigval[:Nstate])
    """df = pd.DataFrame(eigvec.imag[:Nstate,:Nstate])
    print(df)
    df = pd.DataFrame(eigvec.T.conjugate().imag[:Nstate,:Nstate])
    print(df)"""
    eig_state_Sx = np.einsum('pm, mn, nq->pq', eigvec.T.conjugate(), Sx, eigvec)
    eig_state_Sz = np.einsum('pm, mn, nq->pq', eigvec.T.conjugate(), Sz, eigvec)
    eig_state_Lx = np.einsum('pm, mn, nq->pq', eigvec.T.conjugate(), Lx, eigvec)
    eig_state_Lz = np.einsum('pm, mn, nq->pq', eigvec.T.conjugate(), Lz, eigvec)
    tteig_state_Jz = eig_state_Sz + eig_state_Lz
    tteig_state_Jx = eig_state_Sx + eig_state_Lx
    eig_state_Jz = tteig_state_Jz[:Nstate,:Nstate]
    Mj, Jvec = np.linalg.eigh(eig_state_Jz)
    print("#################################################")
    print("#################################################")
    print("Jz eigenvalues from electronic H NjxNj subspace")
    print("#################################################")
    print("#################################################")
    table = PrettyTable(["EIGstate","Mj"])
    for i in range(0,Nstate):
        table.add_row([i+1,Mj[i]])
    print(table)

    J_state_Jx = np.zeros([Nstate,Nstate],dtype=complex)

    #"""adjust the phase of the Jvec
    for i in range(Nstate-1):
        for k in range(Nstate):
            for l in range(Nstate):
                J_state_Jx[i,i+1] = J_state_Jx[i,i+1]+Jvec.T.conjugate()[i,k]*tteig_state_Jx[k,l]*Jvec[l,i+1]
        phi = math.atan2(J_state_Jx.imag[i,i+1],J_state_Jx.real[i,i+1])
        Jvec[...,i+1] = Jvec[...,i+1]*np.exp(-1j*phi)
    #"""
    return Jvec,enr2

def get_Bkq_lsq(ESO,J_state_SOC,lmax,Nstate):
    """
    calculate the Crystal field parameters Bkq using the least square method
    
    Args:
        ESO (np): the matrix of Tensor operators.
        J_state_SOC (np): the matrix of Hamiltonian with the basis of the eigenvector of Jz operator
        lmax (int): the highest order of the Steven's operators.
        Nstate (int): the number of the state of |JM>.

    Returns:
        Bkq (np): the output Crystal field parameters Bkq.
    """
    Odim = int(lmax*(lmax+3)/2)
    Bkq = np.zeros([Odim],dtype=complex)
    A = np.zeros([Nstate*Nstate,Odim],dtype=complex)
    B = np.zeros([Nstate*Nstate],dtype=complex)
    k = 0
    for i in range(Nstate):
        for j in range(Nstate):
            A[k,...] = ESO[i,j,...]
            B[k] = J_state_SOC[i,j]
            k = k+1

    Bkq = zgels(A,B)[1][:Odim]
    """
    Q, R = np.linalg.qr(A)
    Bkq= np.einsum('pm, mq, q->p', np.linalg.inv(R), Q.T, B)
    #QQB = np.einsum('pm, mq, q->p', Q, Q.T, B)
    #print(QQB-B)
    """

    print("#################################################")
    print("#################################################")
    print("Crystal Field Parameters:")
    print("#################################################")
    print("#################################################")
    table = PrettyTable(["k","q","Bkq"])
    i = 0
    for k in range(2,lmax+2,2):
        for q in range(-k,k+1):
            table.add_row([k,q,Bkq.real[i]])
            #print(k,q,Bkq.real[i])
            i = i+1
    print(table)

    return Bkq

def get_Bkq_proj(ESO,J_state_SOC,lmax):
    """
    calculate the Crystal field parameters Bkq using the projection method
    
    Args:
        ESO (np): the matrix of Tensor operators.
        J_state_SOC (np): the matrix of Hamiltonian with the basis of the eigenvector of Jz operator
        lmax (int): the highest order of the Steven's operators.

    Returns:
        Bkq (np): the output Crystal field parameters Bkq.
    """
    Odim = int(lmax*(lmax+3)/2)
    Bkq = np.zeros([Odim],dtype=complex)
    print("#################################################")
    print("#################################################")
    print("Crystal Field Parameters:")
    print("#################################################")
    print("#################################################")
    table = PrettyTable(["k","q","Bkq"])
    i = 0
    for k in range(2,lmax+2,2):
        for q in range(-k,k+1):
            esotemp = ESO[...,i]
            htemp = np.matmul(J_state_SOC,esotemp)
            eso2 = np.matmul(esotemp,esotemp)
            Bkq[i] = htemp.trace().real/eso2.trace().real
            table.add_row([k,q,Bkq.real[i]])
            i = i+1
    print(table)
    return Bkq

def get_expansion(Hmat,Nstate):
    """
    calculate the eigenvector of Hamiltonian using the |JM> basis
    
    Args:
        Hmat (np): the matrix of Tensor operators.
        Nstate (int): the number of the state of |JM>.

    Returns:
        eigvec (np): the eigenvector.
    """
    eigval, eigvec = np.linalg.eigh(Hmat)
    enr2 = np.zeros([Nstate],dtype=complex)

    print("#################################################")
    print("#################################################")
    print("H eigenstates in the Jz eigenstates basis")
    print("#################################################")
    print("#################################################")

    for i in range(0,Nstate):
        enr2[i] = (eigval[i]-eigval[0])
        print("State: ",i+1,"   Energy: ",enr2[i].real," cm^{-1}")
        table = PrettyTable(["M_j","Real","Imag","Percent"])
        for j in range(Nstate):
            Mj = int(2*j-Nstate+1)
            Percent = (eigvec[j,i]*eigvec[j,i].conjugate()).real*100
            table.add_row(["{x}/2".format(x=Mj),eigvec.real[j,i],eigvec.imag[j,i],'%1.2f'%Percent+"%"])
        print(table)

    return eigvec




"""Masoc = np.loadtxt("SOC.txt",dtype=complex)        
Lxt = np.loadtxt("Lx.txt",dtype=complex)
Lyt = np.loadtxt("Ly.txt",dtype=complex)
Lzt = np.loadtxt("Lz.txt",dtype=complex)
Sxt = np.loadtxt("Sx.txt",dtype=complex)
Syt = np.loadtxt("Sy.txt",dtype=complex)
Szt = np.loadtxt("Sz.txt",dtype=complex)
Out_Orca(Lxt[80:,80:].imag,6,'Lx.dat')
Out_Orca(Lzt[80:,80:].imag,6,'Lz.dat')
Out_Orca(Sxt[80:,80:].real,6,'Sx.dat')
Out_Orca(Szt[80:,80:].real,6,'Sz.dat')
Out_Orca(Masoc[80:,80:].real,6,'SOCR.dat')
Out_Orca(Masoc[80:,80:].imag,6,'SOCI.dat')"""

Mosize = 126
Nstate = 16
lmax = int(10)
Odim = int(lmax*(lmax+3)/2)
ncol = 6
ESO_foldname = "inpfile/ESO_mol"
foldname = 'easyorca4f/'
Jmax = (Nstate-1)/2

#ESO = read_ESO_matrix(ESO_foldname,Nstate,lmax)
#ESO = gen_ESO_matrix_easyspin(ESO_foldname,Nstate,lmax)
ESO = gen_ESO_matrix_molfrorge(ESO_foldname,Nstate,lmax)

Masoc,Lx,Lz,Sx,Sz = get_matrix_from_Orca(Mosize,ncol,foldname) 
Jvec,enr2 = get_E_J(Masoc,Lx,Lz,Sx,Sz,Nstate) 
J_state_SOC = np.einsum('pm, m, mq->pq', Jvec.T.conjugate(), enr2, Jvec) #calculate the Hamiltonion using the basis of the eigenstate of the Jx vector#
"""
Bkq = np.zeros([Odim],dtype=complex)
with open("Bkq.dat","r") as file:
    for i in range(Odim):
        line = file.readline()
        Bkq[i] = line.split()[2]
J_state_SOC = np.einsum('n, pqn->pq', Bkq, ESO)
"""
Bkq = get_Bkq_lsq(ESO,J_state_SOC,lmax,Nstate)
Hmat = np.einsum('n, pqn->pq', Bkq, ESO)
eigvec = get_expansion(Hmat,Nstate)





diff = Hmat-J_state_SOC
diffper = np.zeros([Nstate,Nstate],dtype=float)
for i in range(Nstate):
    for j in range(Nstate):
        diffper[i,j] = abs(Hmat[i,j])/abs(J_state_SOC[i,j])
print(Hmat[15,0],J_state_SOC[15,0])
df=pd.DataFrame(diffper)
print(df.applymap(lambda x: '%.3f'%x))
np.savetxt('test.txt',diffper,fmt='%1.2f')



# %%
from matplotlib import pyplot as plt

def is_in_notebook():
    import sys
    return 'ipykernel' in sys.modules

def clear_output():
    """
    clear output for both jupyter notebook and the console
    """
    import os
    os.system('cls' if os.name == 'nt' else 'clear')
    if is_in_notebook():
        from IPython.display import clear_output as clear
        clear()

Mosize = 126
Nstate = 16
lmax = int(8)
Odim = int(lmax*(lmax+3)/2)
ncol = 6
ESO_foldname = "inpfile/ESO"
foldname = '13/eq/'
Jmax = (Nstate-1)/2
ESO = read_ESO_matrix(ESO_foldname,Nstate,lmax)
#ESO = gen_ESO_matrix(ESO_foldname,Nstate,lmax)
theta_2 = -2/315
theta_4 = -8/135135
theta_6 = 4/3864861

Masoc,Lx,Lz,Sx,Sz = get_matrix_from_Orca(Mosize,ncol,foldname) 
Jvec,enr2 = get_E_J(Masoc,Lx,Lz,Sx,Sz,Nstate) 
J_state_SOC = np.einsum('pm, m, mq->pq', Jvec.T.conjugate(), enr2, Jvec) #calculate the Hamiltonion using the basis of the eigenstate of the Jx vector#
Bkq0 = get_Bkq_lsq(ESO,J_state_SOC,lmax,Nstate)

Bkqtt = np.zeros([Odim,11],dtype = float)
for i in range(0,10):
    foldname = '13/1/{x}/'.format(x=i+1)
    Masoc,Lx,Lz,Sx,Sz = get_matrix_from_Orca(Mosize,ncol,foldname) 
    Jvec,enr2 = get_E_J(Masoc,Lx,Lz,Sx,Sz,Nstate) 
    J_state_SOC = np.einsum('pm, m, mq->pq', Jvec.T.conjugate(), enr2, Jvec) #calculate the Hamiltonion using the basis of the eigenstate of the Jx vector#
    Bkq = get_Bkq_proj(ESO,J_state_SOC,lmax,Nstate)
    if i < 5:
        Bkqtt[...,i] = abs(Bkq)-abs(Bkq0)
    else:
        Bkqtt[...,i+1] = abs(Bkq)-abs(Bkq0)

clear_output()
X = [-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1]
for j in range(0,5):
    plt.plot(X,Bkqtt[j,...]/theta_2,label = str(j-2))
plt.xlim(-0.1,0.1)
plt.ylim(-25,25)
plt.legend(loc=0,ncol=1)
plt.title("1-projabs") 
plt.xlabel("x_j") 
plt.ylabel("B_kq(cm^-1)") 
plt.savefig("1-projabs.png".format(f=foldname,x=i+1))
plt.show()






# %%
Bkqtt= np.loadtxt("13/3.txt",dtype=float)
X = [-0.1,-0.08,-0.06,-0.04,-0.02,0,0.02,0.04,0.06,0.08,0.1]
for j in range(0,5):
    plt.plot(X,Bkqtt[j,...]/theta_2,label = str(j-2))
plt.xlim(-0.1,0.1)
plt.ylim(-100,100)
plt.legend(loc=0,ncol=1)
plt.title("3-Mol") 
plt.xlabel("x_j") 
plt.ylabel("B_kq(cm^-1)") 
plt.savefig("3-Mol.png".format(f=foldname,x=i+1))
plt.show()

# %%
from scipy.special import perm, comb

def cleb(j1,m1,j2,m2,j,m):
    factor = 0.0
    factor = comb(j1,(j1+j2-j)/2)/comb((j1+j2+j+2)/2,(j1+j2-j)/2)
    #print(comb(j1,(j1+j2-j)/2))
    factor = factor * comb(j2,(j1+j2-j)/2) / comb(j1,(j1-m1)/2)
    factor = factor / comb(j2,(j2-m2)/2) / comb(j,(j-m)/2)
    factor = np.sqrt(factor)
    
    zmin = int(max(0,j2+(j1-m1)/2-(j1+j2+j)/2,j1+(j2+m2)/2-(j1+j2+j)/2))
    zmax = int(min((j1+j2-j)/2,(j1-m1)/2,(j2+m2)/2))
    
    sum=0.0
    for z in range(zmin,zmax+1):
        par=1
        if (2*(z/2)-int(2*(z/2.0))) != 0:
            par=-1
        sum=sum+par*comb((j1+j2-j)/2,z)*comb((j1-j2+j)/2,(j1-m1)/2-z)*comb((-j1+j2+j)/2,(j2+m2)/2-z)
    #print(factor,sum)
    return factor*sum

def spintens(K,S):
    B=1
    for i in range(-K,K+1):
        B=B*(2*S+1+i)
    B=B/((2**K)*math.factorial(2*K))

    B=np.sqrt(B)
    B=B*math.factorial(K)
    return B


k = 2
q = 0
j1 = j2 = 15/2
m1 = 15/2
m2 = 15/2

cg = CG( S(2*j2)/2, S(2*m2)/2, S(2*k)/2, S(2*q)/2, S(2*j1)/2, S(2*m1)/2, ).doit()
#redelem = sto_redelem(k,j1)
redelem = spintens(k,j1)
#print(redelem)

print(complex(cg)*redelem/(np.sqrt(2*j1+1)))
redelem = sto_redelem(k,j1)
print(complex(cg)*redelem/(np.sqrt(2*j1+1)))


#print(redelem)




#print(comb(5, 2))



