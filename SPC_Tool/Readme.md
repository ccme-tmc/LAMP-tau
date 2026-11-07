# SPC_Tool

`SPC_Tool` 是论文中用于**生成扰动结构、组织 DMET+CASSCF-SO / CASSCF 计算、提取晶体场参数（CFPs）** 的脚本集合。它与仓库根目录的 `mag_anal` 配合使用：`SPC_Tool` 负责产生 `CFPs.out`、`freq.log`、`eqpyout` 等输入，`mag_anal` 负责计算 SPC 参数和弛豫时间 τ(T)。

> English documentation: see [`../docs/workflow.md`](../docs/workflow.md), [`../docs/input_formats.md`](../docs/input_formats.md) and [`../docs/hpc.md`](../docs/hpc.md). The Chinese manual in `manual/manual.tex` is retained for historical reference; where it disagrees with `docs/`, `docs/` is authoritative.

## 目录结构

```
SPC_Tool/
├── inpfile/                 # 输入模板与示例
│   ├── opt_INCAR, freq_INCAR
│   ├── dmet_tmp.py.tmpl     # PySCF/DMET 输入模板（不能直接运行）
│   └── cassi_tmp.py.tmpl    # CASSCF 输入模板（不能直接运行）
├── opt_freq/vaspsh/         # VASP 结构优化与频率计算
├── Per_B_kp/Pyscfsh/        # 生成扰动结构、生成输入、提取 CFPs
├── sh/                      # Slurm/PBS/天河作业模板
├── Relaxation/              # 旧版弛豫示例（inpfile/CFPs_v1.yml）
└── manual/                  # 中文 manual（LaTeX + PDF）
```

## 快速开始

1. 准备 `inpfile/POSCAR`、`opt_INCAR`、`freq_INCAR`。
2. 运行 `opt_freq/vaspsh/opt.sh` 做结构优化。
3. 运行 `opt_freq/vaspsh/freq.sh` 做频率计算，并调用 `Select_group.py` 选择分子基团。
4. 在 `Per_B_kp/Pyscfsh/` 中依次运行：
   - `run1_stru.sh`：生成扰动结构；
   - `run2_Regr.sh`：从周期性结构中提取分子；
   - `run3_inp.sh`：用 `insert_xyz.py` 填充 `dmet_tmp.py.tmpl`；
   - 在集群上运行生成的 `dmet_cassiso.py`；
   - `run4_CFPs.sh`：提取 `CFPs.out`。
5. 回到仓库根目录，运行 `lamp-tau --yml <your.yml>` 计算 τ(T)。

## 重要说明

- `SPC_Tool/sh/*.sh` 是**模板**：请按你的集群修改分区、队列、模块和邮箱；VASP 路径通过 `VASP_BIN` 指定。
- `inpfile/*.tmpl` 是模板文件，不能直接执行；它们由 `insert_xyz.py` 填充。
- DMET+CASSCF-SO 依赖 [LAMP-emb](https://github.com/ZGuan97/LAMP_emb) 的 `embed_sim` 包，请单独安装并加入 `PYTHONPATH`。
- VASP 与 OpenMolcas 是第三方程序，不在本仓库中分发。
- 旧的 `Relaxation/inpfile/CFPs_v1.yml` 仅作格式示例；可运行的示例在 `../examples/relaxation/`。
