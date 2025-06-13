# %%
import numpy as np
import math
#from scipy.sparse.linalg import eigs, eigsh
#from scipy.linalg.lapack import zgels, zgelss
from prettytable import PrettyTable

from sympy.physics.quantum.cg import CG
from sympy import S

def Rot_ESO(Nstate,lmax, ESO, alpha, beta, gamma):
    from sympy.physics.wigner import wigner_d
    #half = (Nstate-1)/Integer(2)
    
    Odim = np.shape(ESO)[2]
    newESO = np.empty([Nstate,Nstate,Odim], dtype = complex)
    tempESO = np.empty([Nstate,Nstate], dtype = complex)
    ESO_list = []
    
    i = 0
    for tempk in range(int(lmax/2)):
        k = int(2*(tempk+1))
        for q in range(-k, k + 1):
            val = ESO[...,i]
            ESO_list.append(StevOper(k, q, val))
            #print(i,k,q)
            i = i+1
    """
    for i in range(Odim):
        newESO[...,i] = np.einsum('pm, mn, nq->pq', rot.T.conjugate(), ESO[...,i], rot)
        print(i)
    """
    STO_list = [StevOper(op.k, op.q, tempESO) for op in ESO_list]
    newSTO_list = [StevOper(op.k, op.q, tempESO) for op in ESO_list]
    newESO_list = [StevOper(op.k, op.q, tempESO) for op in ESO_list]

    for idx, op in enumerate(ESO_list):
        new_op = STO_list[idx]
        new_op.k = op.k
        new_op.q = op.q
        
        if op.q == 0:
            new_op.val = op.val
            #print(op.val)
        elif op.q > 0:
            # 处理正q分量
            #print(op.val)
            neg_op = next(x for x in ESO_list if x.q == -op.q and x.k == op.k)
            new_op.val = op.val - neg_op.val*(-1)*1j
            new_op.val *= (-1)**op.q / np.sqrt(2)
        else:
            # 处理负q分量
            #print(op.val)
            pos_op = next(x for x in ESO_list if x.q == -op.q and x.k == op.k)
            new_op.val = pos_op.val + op.val*(-1)*1j
            new_op.val /= np.sqrt(2)
        #newESO_list[idx] = new_op
    
    # 应用旋转矩阵
    i = 0
    for tempk in range(int(lmax/2)):
        k = int(2*(tempk+1))
        rot = wigner_d(k, alpha, beta, gamma)
        #print(rot)
        #rot = np.array([[0,0,0,0,1],[0,0,0,-1,0],[0,0,1,0,0],[0,-1,0,0,0],[1,0,0,0,0]])
        for q in range(-k, k + 1):
            tempval = 0
            for mq in range(-k, k + 1):
                temp_op = next(x for x in STO_list if x.q == mq and x.k == k)
                #print(mq,temp_op.k,temp_op.q)
                tempval = tempval+rot[q+k,mq+k].conjugate()*temp_op.val
            newSTO_list[i].val = tempval
            #print(q,STO_list[i].q,newSTO_list[i].q)
            #print('new',newSTO_list[i].q,newSTO_list[i].val)
            i = i+1

    
    # 后处理
    for idx, op in enumerate(newESO_list):
        if op.q == 0:
            op.val = newSTO_list[idx].val

        elif op.q > 0:
            neg_op = next(x for x in newSTO_list if x.q == -op.q and x.k == op.k)
            term1 = ((-1)**op.q) * newSTO_list[idx].val
            term2 = neg_op.val
            op.val = (term1 + term2) * np.sqrt(2)/2
            #print(op.val)

        else:
            pos_op = next(x for x in newSTO_list if x.q == -op.q and x.k == op.k)
            term1 = ((-1)**(1-op.q)) * pos_op.val
            term2 = newSTO_list[idx].val
            op.val = (term1 + term2) * np.sqrt(2)/2 * 1j
            #print(op.val)
        #newESO_list[idx] = new_op


    for i in range(Odim):
        #print(ESO_list[i].val)
        #print(newESO_list[i].val-ESO_list[i].val)
        newESO[...,i] = newESO_list[i].val
    return newESO

def calculate_perpendicular_direction(d1, d2):
    """
    if we want to rotate the vector from d2 to d1, we need an axis and an angle. Here calculates the axis 
    """
    normal = np.cross(d1, d2)
    normal = normal / np.linalg.norm(normal)
    return normal

def calculate_angle_between_lines(d1, d2):
    """
    if we want to rotate the vector from d2 to d1, we need an axis and an angle. Here calculates the angle 
    """
    cos_angle = np.dot(d1, d2) / (np.linalg.norm(d1) * np.linalg.norm(d2))
    angle = np.arccos(cos_angle)
    return np.degrees(angle)

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

def detmumu(mui,muj):
    return mui[0,0]*muj[1,1]-mui[1,0]*muj[0,1]


def gen_ESO_matrix_TTO(Nstate,lmax,ESO_foldname=None):
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
    if ESO_foldname is None:
        i = 0
        for k in range(2,lmax+2,2):
            for q in range(-k,k+1):
                for m1 in range(Nstate):
                    Mj1 = m1-J
                    for m2 in range(Nstate):
                        Mj2 = m2-J
                        ESO[m1,m2,i] = setevens_mat(k,q,J,Mj1,J,Mj2)
                i = i+1
        return ESO
    else:
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

class StevOper:
    def __init__(self, k, q, val=0j):
        self.k = k
        self.q = q
        self.val = val

class PSH():
    def __init__(self, SISO,Nstate=16,lmax=12,TFO=1,TFB=0,ROT=False,Angle=None):
        """
        
        Args:
            Nstate (int): the number of the state of |JM>.
            lmax (int): the max rank of tensor operator.
            TFO: the type of the tensor operator. use TFO=0 for extend steven operator from easyspin; use TFO=1 for Tesseral tensor operator from Molforge
            TFB: the method to calculate B_kq, use TFB=0 for Least square method; use TFB=1 for Projection method
            ROT: whether rotate the CFP
            Angle: the zyz eular angle you want to rotate, if Angle is None, it will rotate the easy mag axis to z axiz.


        Returns:
            Jvec (np): the eigenvector of Jz operator.
            enr2 (np): the energy.
        """
        self.ang_mom = SISO.orbital_ang_mom()
        self.sang_mom = SISO.spin_ang_mom()
        self.SOC = SISO.SOC_Hamiltonian
        self.Nstate = Nstate
        self.Jmax = (Nstate-1)/2
        self.lmax = int(lmax)
        self.Odim = int(lmax*(lmax+3)/2)
        self.TFO = TFO
        self.TFB = TFB
        self.ROT = ROT
        self.Angle = Angle
        np.savetxt('SOC.txt',self.SOC,fmt='%.08f')
        np.savetxt('Lx.txt',self.ang_mom[...,0],fmt='%.08f')
        np.savetxt('Ly.txt',self.ang_mom[...,1],fmt='%.08f')
        np.savetxt('Lz.txt',self.ang_mom[...,2],fmt='%.08f')
        np.savetxt('Sx.txt',self.sang_mom[...,0],fmt='%.08f')
        np.savetxt('Sx.txt',self.sang_mom[...,1],fmt='%.08f')
        np.savetxt('Sz.txt',self.sang_mom[...,2],fmt='%.08f')

    def gen_ESO(self):
        if self.TFO==0:
            ESO_foldname = "inpfile/ESO_es"
            print("ESO")
            self.ESO = read_ESO_matrix(self.Nstate,self.lmax,ESO_foldname)
        if self.TFO==1:
            print("TTO")
            self.ESO = gen_ESO_matrix_TTO(self.Nstate,self.lmax)

    def get_E_J(self):
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
        Masoc = self.SOC
        Nstate = self.Nstate
        eigval, eigvec = np.linalg.eigh(Masoc)
        np.savetxt('eigval.txt',eigval,fmt='%.08f')
        enr2 = np.zeros([Nstate],dtype=complex)
        print(f"#################################################")
        print(f"#################################################")
        print(f"Nj eigevanlues of the electronic H")
        print(f"#################################################")
        print(f"#################################################")
        table = PrettyTable(["EIGstate","Energy(cm^{-1})"])
        for i in range(0,Nstate):
            enr2[i] = (eigval[i]-eigval[0])*219474.63
            table.add_row([i+1,enr2[i].real])
        print(f'{table}')
            
        eig_state_ang_mom = np.einsum('pm, mni, nq->pqi', eigvec.T.conjugate(), self.ang_mom, eigvec)
        eig_state_sang_mom = np.einsum('pm, mni, nq->pqi', eigvec.T.conjugate(), self.sang_mom, eigvec)
        eig_state_Sx = eig_state_sang_mom[...,0]
        eig_state_Sy = eig_state_sang_mom[...,1]
        eig_state_Sz = eig_state_sang_mom[...,2]
        eig_state_Lx = eig_state_ang_mom[...,0]
        eig_state_Ly = eig_state_ang_mom[...,1]
        eig_state_Lz = eig_state_ang_mom[...,2]
        np.savetxt('eig_Lx.txt',eig_state_ang_mom[...,0],fmt='%.08f')
        np.savetxt('eig_Ly.txt',eig_state_ang_mom[...,1],fmt='%.08f')
        np.savetxt('eig_Lz.txt',eig_state_ang_mom[...,2],fmt='%.08f')
        np.savetxt('eig_Sx.txt',eig_state_sang_mom[...,0],fmt='%.08f')
        np.savetxt('eig_Sy.txt',eig_state_sang_mom[...,1],fmt='%.08f')
        np.savetxt('eig_Sz.txt',eig_state_sang_mom[...,2],fmt='%.08f')
        tteig_state_Jz = eig_state_Sz*0.5 + eig_state_Lz
        tteig_state_Jx = eig_state_Sx*0.5 + eig_state_Lx
        tteig_mux = eig_state_Sx + eig_state_Lx
        tteig_muy = eig_state_Sy + eig_state_Ly
        tteig_muz = eig_state_Sz + eig_state_Lz
        eig_state_Jz = tteig_state_Jz[:Nstate,:Nstate]
        Mj, Jvec = np.linalg.eigh(eig_state_Jz)
        print(f"#################################################")
        print(f"#################################################")
        print(f"Jz eigenvalues from electronic H NjxNj subspace")
        print(f"#################################################")
        print(f"#################################################")
        table = PrettyTable(["EIGstate","Mj"])
        for i in range(0,Nstate):
            table.add_row([i+1,Mj[i]])
        print(f'{table}')

        J_state_Jx = np.zeros([Nstate,Nstate],dtype=complex)

        #"""adjust the phase of the Jvec
        for i in range(Nstate-1):
            for k in range(Nstate):
                for l in range(Nstate):
                    J_state_Jx[i,i+1] = J_state_Jx[i,i+1]+Jvec.T.conjugate()[i,k]*tteig_state_Jx[k,l]*Jvec[l,i+1]
            phi = math.atan2(J_state_Jx.imag[i,i+1],J_state_Jx.real[i,i+1])
            Jvec[...,i+1] = Jvec[...,i+1]*np.exp(-1j*phi)
        #"""

        J_state_SOC = np.einsum('pm, m, mq->pq', Jvec.T.conjugate(), enr2, Jvec) #calculate the Hamiltonion using the basis of the eigenstate of the Jx vector#
        self.eig_state_Jz=eig_state_Jz
        self.Jvec=Jvec
        self.enr2=enr2
        self.Mj=Mj
        self.J_state_SOC=J_state_SOC
        self.tteig_mux=tteig_mux
        self.tteig_muy=tteig_muy
        self.tteig_muz=tteig_muz
        return 
    
    def get_E_J_rot(self):
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
        Masoc = self.SOC
        Nstate = self.Nstate
        eigval, eigvec = np.linalg.eigh(Masoc)
        np.savetxt('eigval.txt',eigval,fmt='%.08f')
        enr2 = np.zeros([Nstate],dtype=complex)
        print(f"#################################################")
        print(f"#################################################")
        print(f"Nj eigevanlues of the electronic H")
        print(f"#################################################")
        print(f"#################################################")
        table = PrettyTable(["EIGstate","Energy(cm^{-1})"])
        for i in range(0,Nstate):
            enr2[i] = (eigval[i]-eigval[0])*219474.63
            table.add_row([i+1,enr2[i].real])
        print(f'{table}')
            
        eig_state_ang_mom = np.einsum('pm, mni, nq->pqi', eigvec.T.conjugate(), self.ang_mom, eigvec)
        eig_state_sang_mom = np.einsum('pm, mni, nq->pqi', eigvec.T.conjugate(), self.sang_mom, eigvec)
        eig_state_Sx = eig_state_sang_mom[...,0]
        eig_state_Sy = eig_state_sang_mom[...,1]
        eig_state_Sz = eig_state_sang_mom[...,2]
        eig_state_Lx = eig_state_ang_mom[...,0]
        eig_state_Ly = eig_state_ang_mom[...,1]
        eig_state_Lz = eig_state_ang_mom[...,2]
        np.savetxt('eig_Lx.txt',eig_state_ang_mom[...,0],fmt='%.08f')
        np.savetxt('eig_Ly.txt',eig_state_ang_mom[...,1],fmt='%.08f')
        np.savetxt('eig_Lz.txt',eig_state_ang_mom[...,2],fmt='%.08f')
        np.savetxt('eig_Sx.txt',eig_state_sang_mom[...,0],fmt='%.08f')
        np.savetxt('eig_Sy.txt',eig_state_sang_mom[...,1],fmt='%.08f')
        np.savetxt('eig_Sz.txt',eig_state_sang_mom[...,2],fmt='%.08f')
        tteig_state_Jz = eig_state_Sz*0.5 + eig_state_Lz
        tteig_state_Jy = eig_state_Sy*0.5 + eig_state_Ly
        tteig_state_Jx = eig_state_Sx*0.5 + eig_state_Lx

        roteig_state_Jx = self.gdgvec[0,0]*tteig_state_Jx+self.gdgvec[1,0]*tteig_state_Jy+self.gdgvec[2,0]*tteig_state_Jz
        roteig_state_Jy = self.gdgvec[0,1]*tteig_state_Jx+self.gdgvec[1,1]*tteig_state_Jy+self.gdgvec[2,1]*tteig_state_Jz
        roteig_state_Jz = self.gdgvec[0,2]*tteig_state_Jx+self.gdgvec[1,2]*tteig_state_Jy+self.gdgvec[2,2]*tteig_state_Jz
        eig_state_Jz = roteig_state_Jz[:Nstate,:Nstate]
        eig_state_Jx = roteig_state_Jx[:Nstate,:Nstate]
        print("Finish read")
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
                    J_state_Jx[i,i+1] = J_state_Jx[i,i+1]+Jvec.T.conjugate()[i,k]*eig_state_Jx[k,l]*Jvec[l,i+1]
            phi = math.atan2(J_state_Jx.imag[i,i+1],J_state_Jx.real[i,i+1])
            print(phi)
            Jvec[...,i+1] = Jvec[...,i+1]*np.exp(-1j*phi)
        #"""

        J_state_SOC = np.einsum('pm, m, mq->pq', Jvec.T.conjugate(), enr2, Jvec)
        self.eig_state_Jz=eig_state_Jz
        self.Jvec=Jvec
        self.enr2=enr2
        self.Mj=Mj
        self.J_state_SOC=J_state_SOC

        return Jvec,enr2,eig_state_Jz,Mj

    def get_g_tensor(self):
        Nstate = self.Nstate
        eig_muz = self.tteig_muz[:Nstate,:Nstate]
        eig_mux = self.tteig_mux[:Nstate,:Nstate]
        eig_muy = self.tteig_muy[:Nstate,:Nstate]
        Ma_G = np.empty([3,3], dtype = float)
        gclist=['gX','gY','gZ']
        cerclist=['Xm','Ym','Zm']
        for KDi in range(0,int(Nstate/2)):
            print("#################################################")
            print("g tensor for MULTIPLET {x}:".format(x=KDi+1))
            print("#################################################")
            rows = np.array([[2*KDi,2*KDi],[2*KDi+1,2*KDi+1]]) 
            cols = np.array([[2*KDi,2*KDi+1],[2*KDi,2*KDi+1]])
            #rows = np.array([[KDi,KDi],[Nstate-1-KDi,Nstate-1-KDi]]) 
            #cols = np.array([[KDi,Nstate-1-KDi],[KDi,Nstate-1-KDi]]) 
            peig_mux = eig_mux[rows,cols]  
            peig_muy = eig_muy[rows,cols]  
            peig_muz = eig_muz[rows,cols]
            #print(peig_Sz)
            peig_mu = np.stack((peig_mux,peig_muy,peig_muz),0)
            #print(peig_mu)
            for i in range(0,3):
                for j in range(0,3):
                    Ma_G[i,j] = -0.5*(detmumu(peig_mu[i,...],peig_mu[j,...])+detmumu(peig_mu[j,...],peig_mu[i,...]))
            Aval, gvec = np.linalg.eigh(Ma_G)
            gval = 2*np.sqrt(np.abs(Aval))
            table = PrettyTable(["","Value","Axes","x","y","z"])
            #print(Ma_G)
            if KDi == 0:
                self.gdgvec = gvec
            for i in range(0,3):
                table.add_row([gclist[i],gval[i],cerclist[i],gvec[0,i],gvec[1,i],gvec[2,i]])
            print(table)
            print("x , y , z  -- initial Cartesian axes")
            print("Xm, Ym, Zm -- main magnetic axes")
    
    def get_Bkq(self):
        if self.TFB==0:
            self.get_Bkq_lsq(0)
        if self.TFB==1:
            self.get_Bkq_proj(0)
        self.get_expansion(self.J_state_SOC)

    def get_rot_Bkq(self):
        if self.TFB==0:
            self.get_Bkq_lsq(0)
        if self.TFB==1:
            self.get_Bkq_proj(0)
        self.rot_Bkq()
        if self.TFB==0:
            self.get_Bkq_lsq(1)
        if self.TFB==1:
            self.get_Bkq_proj(1)
        Hmat = np.einsum('n, pqn->pq', self.Bkq, self.ESO)
        self.get_expansion(Hmat)

    def get_Bkq_lsq(self,rot=0):
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
        if rot == 0:
            ESO = self.ESO
        else:
            ESO = self.rotESO
        J_state_SOC=self.J_state_SOC
        lmax=self.lmax
        Nstate=self.Nstate
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

        Q, R = np.linalg.qr(A)
        Bkq= np.einsum('pm, mq, q->p', np.linalg.inv(R), Q.T, B)

        """
        Q, R = np.linalg.qr(A)
        Bkq= np.einsum('pm, mq, q->p', np.linalg.inv(R), Q.T, B)
        Bkq = zgels(A,B)[1][:Odim]
        """

        print(f"#################################################")
        print(f"#################################################")
        if rot == 0:
            print(f"Crystal Field Parameters:")
        else:
            print(f"Rotated Crystal Field Parameters:")
        print(f"#################################################")
        print(f"#################################################")
        table = PrettyTable(["k","q","Bkq"])
        i = 0
        for k in range(2,lmax+2,2):
            for q in range(-k,k+1):
                table.add_row([k,q,Bkq.real[i]])
                #print(k,q,Bkq.real[i])
                i = i+1
        print(f'{table}')
        self.Bkq = Bkq
        return 

    def get_Bkq_proj(self,rot=0):
        """
        calculate the Crystal field parameters Bkq using the projection method
        
        Args:
            ESO (np): the matrix of Tensor operators.
            J_state_SOC (np): the matrix of Hamiltonian with the basis of the eigenvector of Jz operator
            lmax (int): the highest order of the Steven's operators.

        Returns:
            Bkq (np): the output Crystal field parameters Bkq.
        """
        if rot == 0:
            ESO = self.ESO
        else:
            ESO = self.rotESO
        J_state_SOC=self.J_state_SOC
        lmax=self.lmax
        Odim = int(lmax*(lmax+3)/2)
        Bkq = np.zeros([Odim],dtype=complex)
        print(f"#################################################")
        print(f"#################################################")
        if rot == 0:
            print(f"Crystal Field Parameters:")
        else:
            print(f"Rotated Crystal Field Parameters:")
        print(f"#################################################")
        print(f"#################################################")
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
        print(f'{table}')
        self.Bkq = Bkq
        return 

    def get_expansion(self,Hmat):
        """
        calculate the eigenvector of Hamiltonian using the |JM> basis
        
        Args:
            Hmat (np): the matrix of Tensor operators.
            Nstate (int): the number of the state of |JM>.

        Returns:
            eigvec (np): the eigenvector.
        """
        Nstate = self.Nstate
        eigval, eigvec = np.linalg.eigh(Hmat)
        enr2 = np.zeros([Nstate],dtype=complex)

        print(f"#################################################")
        print(f"#################################################")
        print(f"H eigenstates in the Jz eigenstates basis")
        print(f"#################################################")
        print(f"#################################################")

        for i in range(0,Nstate):
            enr2[i] = (eigval[i]-eigval[0])
            print(f"State: {i+1}   Energy: {enr2[i].real} cm-1")
            table = PrettyTable(["M_j","Real","Imag","Percent"])
            for j in range(Nstate):
                Mj = int(2*j-Nstate+1)
                Percent = (eigvec[j,i]*eigvec[j,i].conjugate()).real*100
                table.add_row(["{x}/2".format(x=Mj),eigvec.real[j,i],eigvec.imag[j,i],'%1.2f'%Percent+"%"])
            print(f'{table}')
        self.expan = eigvec
        return 
    
    def rot_Bkq(self):
        from scipy.spatial.transform import Rotation
        axis=calculate_perpendicular_direction(self.gdgvec[...,2], [0,0,1])
        theta=calculate_angle_between_lines(self.gdgvec[...,2], [0,0,1])
        rotation_matrix = Rotation.from_rotvec(theta*axis,degrees=True).as_matrix()
        euler_angles = Rotation.from_matrix(rotation_matrix).as_euler('zyz',degrees=False)
        #euler_angles = euler_angles*np.pi/180
        print('the euler_angles(zyz) is',euler_angles)

        if self.Angle is None:
            self.rotESO = Rot_ESO(self.Nstate,self.lmax, self.ESO, euler_angles[0], euler_angles[1], euler_angles[2])
        else:
            self.rotESO = Rot_ESO(self.Nstate,self.lmax, self.ESO, self.Angle[0], self.Angle[1], self.Angle[2])

    def kernel(self):
        self.gen_ESO()
        self.get_E_J()
        self.get_g_tensor()
        if self.ROT:
            """
            self.get_E_J_rot()
            self.get_Bkq()
            """
            self.get_rot_Bkq()
        else:
            self.get_Bkq()
        return





if __name__ == "__main__":
    from pyscf import scf, gto
    from embed_sim import myavas, sacasscf_mixer, siso

    title = 'CoSH4'
    mol = gto.M(atom = '''
               Co     0.005942    0.005942   -0.004764
               S      0.189139    0.186135    2.381110
               S      1.656432    1.659547   -0.548211
               S      0.353936   -2.181551   -0.929305
               S     -2.183408    0.352012   -0.925639
               H     -0.746735   -0.744550    2.637044
               H      1.562463   -2.352031   -0.365709
               H     -2.350968    1.563598   -0.367756
               H      1.513199    1.510898   -1.876769
     ''',
     basis={'default':'def2tzvp'}, symmetry=0 ,spin = 3,charge = -2,verbose= 4)


    mf = scf.rohf.ROHF(mol).x2c()

    chk_fname = title + '_rohf.chk'

    mf.chkfile = chk_fname
    mf.init_guess = 'chk'
    mf.level_shift = .1
    mf.max_cycle = 1000
    mf.max_memory = 100000
    mf.kernel()

    ncas, nelec, mo = myavas.avas(mf, 'Co 3d', threshold=0.5, minao='def2tzvp' ,canonicalize=False)

    mycas = sacasscf_mixer.sacasscf_mixer(mf, ncas, nelec, statelis=[0, 40, 0, 10])
    mycas.kernel(mo)

    mysiso = siso.SISO(title, mycas)
    mysiso.kernel()

    mypsh = PSH(mysiso,Nstate=int(16),lmax=int(12),TFO=1,TFB=0,ROT=False)
    mypsh.kernel()
    

    

