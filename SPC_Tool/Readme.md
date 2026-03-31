！！！！！！！！！！！！！！！！！！！！
该文件并未及时更新，具体内容以manual为准
！！！！！！！！！！！！！！！！！！！！

这个脚本包适用于在天河超算上进行单分子磁体的弛豫时间的计算
首先你需要在inpfile文件夹里准备以下内容：

1.POSCAR
目标SMMs的原始结构，也可以使用sh文件夹里的POSCAR.sh来生成

2.opt_INCAR
用于计算结构优化的VASP的INCAR文件

3.freq_INCAR
用于计算频率的VASP的INCAR文件

4.job_vasp.sh
用于提交VASP计算的文件，可以修改#SBATCH -J以改变任务名称，以更方便地在任务列表中查询

5.cassi_tmp.py; dmet_tmp.py
用于提交磁性计算的文件

其次你需要在Pyscf文件夹中修改部分数据

1.run1_stru.sh
修改npts以决定你希望每一个扭曲模式有多少个取点；
mdp一般不需要修改，其详细意义见下文

2.run2_inp.sh
修改npts与run1_stru.sh中的一致
修改Title命名你的体系（必须）
修改Charge，其物理意义是你所选基团的电荷

3.job/job.sh
修改npts与run1_stru.sh中的一致
修改Title与run2_inp.sh中的一致

4.job/jobeq.sh
修改Title与run2_inp.sh中的一致






文件目录：
~/
 |---  vaspsh/
 |---  inpfile/
 |---  pyscfsh/
		|--- job/
	
	
使用方式：
一、VASP过程
1.运行opt.sh
该脚本会在根目录下新建tmpfile/和VASP/opt/文件夹，用于进行结构优化计算。
该过程会调用POSCAR，opt_INCAR文件，使用vaspkit生成POTCAR和KPOINTS。并根据job_vasp.sh中的设置提交计算

2.运行freq.sh
该脚本会在根目录下新建VASP/freq/文件夹，用于进行振动模计算。
该过程会首先判断VASP/opt/的计算是否收敛，若不收敛则会重复提交。若收敛，则根据优化计算的CONTCAR，使用Select_group.py选择出合适的目标基团，在inpfile中增加opt_CONTCAR，freq_POSCAR，dpxs，Atom_list和eq_opt.xyz。该过程会调用freq_INCAR文件，复制VASP/opt/中的POTCAR和KPOINTS。并根据job_vasp.sh中的设置提交计算。
Select_group.py的运算逻辑为，以指定分子(--ct Co 指以第一个Co原子为中心)为中心，以2.5A为半径进行第一轮原子集选择，之后将以每个原子集中的原子以1.8A为半径进行选择，直到原子集数量不变为止。如果有不同的选取需求请自行查看Select_group.py里的注释


二、Pyscf过程

1.run1_stru.sh
生成用于计算微扰哈密顿量的结构。主要调用gen_structures_cell.py的功能
gen_structures_cell.py的-x：沿着xyz轴扭曲，生成结构。该方法调用inpfile/eq_opt.xyz。-mdp 1指最大位移为0.1A
gen_structures_cell.py的-f：沿着简振模扭曲，生成初步结构，该方法初步生成的结构为晶格数据，需要进一步使用Replace_group.py生成用于最终计算的结构。该方法调用inpfile下的dpxs,opt_CONTCAR,Atom_list。-mdp 1指最大位移为零点振动模

2.run2_inp.sh
生成上述结构对应的计算文件，这一过程确定了计算文件中的结构读取，电荷以及体系命名，脚本会将job文件夹中的.sh脚本移动到指定位置。



