import numpy as np
from pyscf import scf
from pyscf import df
from embed_sim import myavas, siso, sacasscf_mixer, rdiis
from mag_anal import psh
import time
import basis_set_exchange as bse
from pyscf.gto.basis import parse_gaussian
from pyscf import gto


Molecular_Name
Molecular_Charge
Calculate_Title

mol = gto.M(atom = m_name,
     basis={
            'Dy':parse_gaussian.parse(bse.get_basis('x2c-TZVPall-2c',elements='Dy',fmt='gaussian94')),
            'C' :parse_gaussian.parse(bse.get_basis('x2c-SVPall-2c' ,elements='C' ,fmt='gaussian94')),
            'H' :parse_gaussian.parse(bse.get_basis('x2c-SVPall-2c' ,elements='H' ,fmt='gaussian94'))
            },
     symmetry=0 ,spin = 5 ,charge = m_ch,verbose= 4)
     # basis：基函数，spin：未成对电子数，charge：总电荷，verbose：数字越大输出信息越详细

rohf_start = time.time()
auxbas = df.autoaux(mol)
mf = scf.rohf.ROHF(mol).x2c().density_fit(auxbasis = auxbas)
# 有未成对电子，所以用ROHF方法，没有自选污染，无需UHF方法；x2c：标量相对论效应修正
chk_fname = title + '_rohf.chk'
# 生成断点文件
#Hartree-Fock方法参数设置
mf.diis = rdiis.RDIIS(rdiis_prop='dS',imp_idx=mol.search_ao_label(['Dy.*']),power=0.2)
mf.chkfile = chk_fname
mf.init_guess = 'chk'
mf.level_shift = 1   # 用于收敛的参数
mf.max_cycle = 1000   # 最大迭代次数
#mf.max_memory = 100000 # 最大内存
mf.conv_tol = 1e-6   # 收敛精度
mf.kernel() # 核函数，是该函数真正执行计算的函数


rohf_end = time.time()
print('ROHF time', rohf_end - rohf_start)
cas_start = time.time()
# 用avas方法产生初猜轨道
ncas, nelec, mo = myavas.avas(mf, 'Dy 4f', minao='x2c-TZVPall-2c', openshell_option=3, threshold=0.5)
# ncas：活性空间轨道数，nelec：活性空间电子数，mo：轨道系数

mycas = sacasscf_mixer.sacasscf_mixer(mf, ncas, nelec, statelis=[0, 0, 0, 0, 0, 21]) # statelis：激发态列表，从左往右四个数字分别代表自旋单重态，双重态，三重态，四重态，其数目可用轨道和电子排列组合来计算
cas_result = mycas.kernel(mo)

# e_corr = sacasscf_mixer.sacasscf_nevpt2(mycas) # NEVPT2对能量进行校正
# mycas.fcisolver.e_states = mycas.fcisolver.e_states + e_corr
# print('CASSCF excitation_ene')
# print(np.sort(es_cas.e_states) - np.min(np.array(es_cas.e_states)))
cas_end = time.time()
print('CAS time', cas_end - cas_start)
"""
dmet_nevpt2_start = time.time()
e_corr = sacasscf_mixer.sacasscf_nevpt2(es_cas)
es_cas.fcisolver.e_states = es_cas.fcisolver.e_states + e_corr
print('DMET-NEVPT2 excitation_ene')
print(np.sort(es_cas.e_states) - np.min(np.array(es_cas.e_states)))
dmet_nevpt2_end = time.time()
print('DMET-NEVPT2 time', dmet_nevpt2_end - dmet_nevpt2_start)
"""
siso_start = time.time()
mysiso = siso.SISO(title, mycas)
mysiso.kernel()
siso_end = time.time()
print('SISO time', siso_end - siso_start)




# e_corr = sacasscf_mixer.sacasscf_nevpt2(mycas) # NEVPT2对能量进行校正
# mycas.fcisolver.e_states = mycas.fcisolver.e_states + e_corr

# 旋轨耦合会造成自旋简并度的破缺，所以需要用SISO方法进行计算 

psh_start = time.time()
mypsh = psh.PSH(mysiso,Nstate=int(16),lmax=int(12),TFO=1,TFB=0,ROT=True)
mypsh.kernel()
psh_end = time.time()
print('PSH time', psh_end - psh_start)

