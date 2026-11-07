# LAMP-tau User Manual

> **中文摘要**：本手册面向 LAMP-tau 的使用者，按“安装 → 快速开始 → 完整 SPC 工作流 → YAML/文件格式 → 输出文件 → Python API → HPC → 故障排查”的顺序组织。可运行的 DyAC 示例在 `examples/relaxation/`；英文文档 `docs/` 与本书内容一致。若只关心弛豫时间分析，直接阅读第 4、6、8 节即可。

Version 1.0.0 — accompanying the manuscript
*Accelerating ab initio spin-phonon relaxation simulation of single-ion magnets by quantum embedding and spatial truncation* (J. Chem. Phys., 2026).

---

## 1. Overview

LAMP-tau is a Python toolkit for *ab initio* spin-phonon relaxation simulation of single-ion magnets (SIMs). It has two parts:

- **`mag_anal`** — analysis of electronic-structure results: crystal-field parameters (CFPs), effective spin Hamiltonian, *g* tensors, spin-phonon coupling (SPC) matrix elements and relaxation times τ(T).
- **`SPC_Tool`** — scripts that prepare and organise the underlying calculations: VASP phonons, displaced structures, DMET+CASSCF-SO (PySCF) or CASSCF (OpenMolcas) inputs, and CFP extraction.

The electronic-structure workflow is

```
VASP phonons ──► displaced structures ──► DMET+CASSCF-SO / CASSCF ──► CFPs
                                                                      │
                                                     SPC parameters ◄─┘
                                                            │
                                                            ▼
                                                         τ(T), Ueff, τ0
```

`mag_anal` is the part used in step 4–5; `SPC_Tool` covers steps 1–3.

### 1.1 Requirements at a glance

| Item | Requirement |
|---|---|
| Python | ≥ 3.10 (the code uses `match`/`case`) |
| Core Python packages | NumPy, SciPy, PyYAML, prettytable, Matplotlib, ASE, SymPy |
| DMET+CASSCF-SO | PySCF, `basis_set_exchange`, [LAMP-emb](https://github.com/ZGuan97/LAMP_emb) (`embed_sim`) |
| Phonons | VASP (periodic) or a molecular Hessian (e.g. ORCA) |
| Alternative electronic structure | OpenMolcas (for the Molcas path) |
| Scheduler | Slurm, PBS or Tianhe; the shipped `.sh` files are templates |

VASP, OpenMolcas and ORCA are external programs and are **not** distributed with LAMP-tau.

---

## 2. Installation

### 2.1 From a clone

```bash
git clone https://github.com/Dengyf98/LAMP_tau.git
cd LAMP_tau
python -m pip install -e .
```

This installs the `mag_anal` package and the `lamp-tau` command. If you only need the analysis part:

```bash
python -m pip install -r requirements.txt
```

### 2.2 LAMP-emb and PySCF

The DMET part uses the `embed_sim` package from LAMP-emb:

```bash
git clone https://github.com/ZGuan97/LAMP_emb.git
export PYTHONPATH=/path/to/LAMP_emb:$PYTHONPATH
python -m pip install pyscf basis_set_exchange
```

`dmet_tmp.py.tmpl` imports `embed_sim.myavas`, `ssdmet`, `siso`, `sacasscf_mixer` and `rdiis`.

### 2.3 Check the installation

```bash
lamp-tau --help
python -c "import mag_anal; print(mag_anal.__version__)"
```

On a headless machine set:

```bash
export MPLBACKEND=Agg
```

### 2.4 Optional development tools

```bash
python -m pip install -e ".[dev]"   # pytest, ruff
python -m pytest tests -q
```

---

## 3. Repository layout

```
LAMP_tau/
├── mag_anal/
│   ├── Anal_Tau.py          # τ(T) analysis; CLI entry point
│   ├── psh.py               # pseudospin Hamiltonian, g tensors, CFP projection
│   ├── gen_structures.py    # displaced structures from phonon modes
│   ├── Select_group.py      # molecular-group selection / spatial truncation
│   └── dictionary.py        # element/mass dictionaries
├── SPC_Tool/
│   ├── inpfile/             # input templates (VASP, PySCF, OpenMolcas)
│   ├── opt_freq/vaspsh/     # VASP geometry optimisation and frequencies
│   ├── Per_B_kp/Pyscfsh/    # displaced structures, inputs, CFP extraction
│   ├── sh/                  # scheduler templates
│   ├── Relaxation/          # legacy example configuration
│   └── manual/              # Chinese manual (LaTeX + PDF)
├── examples/relaxation/     # runnable DyAC example
├── docs/                    # workflow, formats, HPC and data notes
├── tests/                   # smoke test
├── pyproject.toml
├── requirements.txt
├── CITATION.cff
└── LICENSE
```

---

## 4. Quick start: relaxation analysis (DyAC)

The shipped example uses DyAC and runs from the repository root.

```bash
lamp-tau --yml examples/relaxation/CFPs_v1.yml
```

or, without installing:

```bash
python -m mag_anal.Anal_Tau --yml examples/relaxation/CFPs_v1.yml
```

What happens:

1. `Read_input` parses the YAML file and creates `Out_foldname`.
2. `Read_freq` reads `freq.log`/`POSCAR` and builds the phonon modes.
3. `Read_CFPs` reads `CFPs_py.out` and fits ∂B<sub>kq</sub>/∂Q.
4. `Read_Static_effH` reads the equilibrium output (`eqpyout`) and the tensor-operator matrices.
5. `Analyse_Relax` computes the spin-phonon matrix, transition probabilities and τ(T).

Expected console tail:

```
TestU_{eff} of Orbach is 190.1 K
log(tau) is -11.30 s
rate-determing step
```

Outputs (in `examples/relaxation/TPs_file/` by default):

| File | Content |
|---|---|
| `OB_v1.dat` | `1/T`, `-log10(τ/s)`, standard deviation |
| `TPs1.png` | `log10(τ/s)` vs `1/T` (Orbach) |
| `TPs2.png` | comparison figure |
| `VP.log` | spin-phonon matrix elements |
| `TPs/*.tps`, `TPs/*.tp` | cached transition matrices per temperature |

The full example takes a few minutes on a laptop. For a faster check, reduce `Trange1` to two temperatures.

### 4.1 Programmatic use

```python
from mag_anal.Anal_Tau import CAL_RELAX

tau = CAL_RELAX("examples/relaxation/CFPs_v1.yml")
tau.kernel()

# Orbach rate-determining analysis for the first temperature
tau.TP[0].Mod_Analyse_OB()
```

`tau.TP` is a list of `TPs_Analyze` objects, one per temperature.

---

## 5. Full SPC workflow

### 5.1 Step 1 — geometry optimisation and phonons

Prepare in `SPC_Tool/inpfile/`:

- `POSCAR` (or a `.cif` for `gen_POSCAR.py`);
- `opt_INCAR`, `freq_INCAR`;
- `job_vasp.sh` adapted to your cluster.

Run, in `SPC_Tool/opt_freq/vaspsh/`:

```bash
bash POSCAR.sh     # optional: build POSCAR from a structure
bash opt.sh        # VASP geometry optimisation
bash freq.sh       # VASP frequencies; calls Select_group.py
```

`freq.sh` checks that the optimisation reached `reached required accuracy`, then creates `VASP/freq/`, runs `Select_group.py` and copies `freq_POSCAR`, `eq_opt.xyz`, `dpxs` and `Atom_list` into `inpfile/`.

`Select_group.py` options:

| Option | Meaning |
|---|---|
| `--cry` | crystal structure file |
| `--ct` | centre element (e.g. `Dy`) |
| `--num` | centre atom index (optional) |
| `--out` | output folder |
| `-f`, `--freq` | also write a `Selective dynamics` `freq_POSCAR` |

### 5.2 Step 2 — displaced structures

In `SPC_Tool/Per_B_kp/Pyscfsh/`:

```bash
bash run1_stru.sh
bash run2_Regr.sh   # only for the periodic-cell path
```

`run1_stru.sh` calls either:

```bash
python3 gen_structures_xyz.py --pos ../../inpfile/eq_opt.xyz -np 10 -mdp 0.02
# or
python3 gen_structures_cell.py --out ../../inpfile/freq.log --pos ../../inpfile/opt_CONTCAR -np 10 -mdp 0.2
```

| Script | Options |
|---|---|
| `gen_structures_xyz.py` | `--pos`, `-np/--npts`, `-mdp/--maxdp` |
| `gen_structures_cell.py` | `--out` (VASP frequency file), `--pos`, `-np`, `-mdp`, `-v` |
| `Replace_group.py` | `-i/--inp`, `-o/--out`, `--inpfile` |

`gen_structures_*` writes `perb/<mode>/<point>/mol_<mode>_<point>.xyz`. `Replace_group.py` extracts the molecular group from the periodic cell using `inpfile/Atom_list` and `inpfile/dpxs`.

### 5.3 Step 3 — electronic-structure inputs

```bash
bash run3_inp.sh
```

`run3_inp.sh` fills the template for the equilibrium structure and every displacement:

```bash
python insert_xyz.py --tmp ../../inpfile/dmet_tmp.py.tmpl \
    --xyz mol_<mode>_<point>.xyz -o dmet_cassiso.py -c <charge> -t <title>
```

`insert_xyz.py` options:

| Option | Meaning |
|---|---|
| `--tmp` | template file (`.tmpl`) |
| `--xyz` | structure file |
| `-o/--out` | output Python input |
| `-c/--charge` | total charge |
| `-t/--title` | job title |

The resulting `dmet_cassiso.py` still needs a **complete basis block**. The template contains Dy/C/H entries as an example; add every element present in your system.

For OpenMolcas, use `SPC_Tool/inpfile/molcas_cas.in` and the Molcas job templates instead.

### 5.4 Step 4 — run the electronic-structure jobs

The generated layout is:

```
perb/
├── eq/
│   ├── dmet_cassiso.py
│   └── jobeq.sh
└── <mode>/<point>/
    ├── dmet_cassiso.py
    └── pyjob.sh
```

Submit `jobeq.sh` first, then the displaced jobs. `SPC_Tool/sh/` contains templates:

| Script | Purpose |
|---|---|
| `job_pyscf.sh` | one displaced PySCF job |
| `jobeq_pyscf.sh` | equilibrium PySCF job |
| `job_molcas.sh`, `jobeq_molcas.sh` | OpenMolcas equivalents |
| `check_pyscf.sh`, `check_molcas.sh` | check completion and trim output |
| `clean_pyscf.sh` | remove checkpoint files |
| `yhjob_pyscf.sh`, `yhjob_molcas.sh` | batch submission |

These files are **templates**; edit partition/queue, e-mail, modules and `VASP_BIN`. See §10.

### 5.5 Step 5 — extract CFPs

```bash
bash run4_CFPs.sh
```

which is equivalent to:

```bash
python3 Extract_CFPs.py -el perb/eq/pyout -nd 201 -np 10
```

| Option | Meaning |
|---|---|
| `-el/--eqlog` | equilibrium magnetic output |
| `-nd/--nmds` | number of vibrational modes |
| `-np/--npts` | number of displacement points |
| `--prefix` | directory containing `perb` (default `perb`) |
| `--stdout` | per-point output name (default `stdout1`) |
| `--out` | output CFP file (default `CFPs.out`) |
| `-v` | verbose |

The result `CFPs.out` is the `CFPs_filename` consumed by `mag_anal`.

### 5.6 Step 6 — relaxation time

Create a YAML file pointing to `CFPs.out`, the equilibrium output, `freq.log`, `POSCAR`, `Atom_list` and the tensor-operator folder, then run:

```bash
lamp-tau --yml my_system.yml
```

See §6 for all keys and §8 for outputs.

---

## 6. YAML reference

A complete example is `examples/relaxation/CFPs_v1.yml`.

### 6.1 Required keys

| Key | Type | Meaning |
|---|---|---|
| `CFPs_filename` | path | `CFPs.out` from `Extract_CFPs.py` |
| `mag_filename` | path | equilibrium electronic-structure output |
| `Ab_Type` | int | `0` Molcas, `1` OpenMolcas, `3` PySCF |
| `freq_filename` | path | VASP frequency output |
| `POSCAR` | path | equilibrium structure |
| `TO_foldname` | path | tensor-operator folder |
| `VP_filename` | path | spin-phonon matrix output/prefix |
| `TPs_foldname` | path | cache folder for transition matrices |
| `Out_foldname` | path | folder for τ(T) outputs |
| `CFPs_X` | list | displacement grid; must contain `0` |
| `Correction_slope`, `Corrected_intercept` | float | frequency correction `ν → a·ν + b` |
| `Fitting_order` | int | polynomial order of the CFP fit |
| `Linewidth_parameter` | float | broadening (cm⁻¹) |
| `If_Data` | bool | overlay experimental τ(T) |
| `If_Perb` | bool | apply the perturbation treatment |
| `Perb_T` | float | perturbation field (T) |
| `g_z` | float | ground-state `g_z` |
| `If_uniformization` | bool | use the uniformization method |
| `CFPs_Type` | int | see §6.2 |
| `If_get_picture_Tt` | bool | write τ(T) figures |
| `If_Orbach`, `If_Raman` | bool | enable relaxation channels |

### 6.2 `CFPs_Type`

| Value | Meaning | Extra keys |
|---|---|---|
| `0` | normal-mode CFPs | none |
| `1` | Cartesian (xyz) CFPs | `Atom_list` |
| `2` | xyz + spatial reduction | `Molfile`, plus `Atom_radius` **or** `Atom_label` |
| `3` | Lorentz random CFPs | `llzA`, `llzw` |

For `CFPs_Type=1` and `2`, `Atom_list` points to the atom indices retained by the spatial-truncation selection. The default is `inpfile/Atom_list`; specify it explicitly for portable examples.

### 6.3 Optional keys

| Key | Meaning |
|---|---|
| `Atom_list` | atom indices for the mode projection (default `inpfile/Atom_list`) |
| `First_list` | first-shell indices for `CFPs_Type=2` |
| `Molfile`, `ct_idx`, `Atom_radius`, `Atom_label` | spatial-reduction options |
| `llzA`, `llzw` | Lorentz random-CFP parameters |
| `If_select_sigma`, `select_sigma` | restrict modes to a frequency window (cm⁻¹) |
| `If_T_step`, `nTrange`, `Trange1..n` | temperature grid |
| `ndata`, `Tr1..`, `logtime1..`, `label1..` | experimental data overlay |

### 6.4 Temperature grid

- `If_T_step: True`: each `TrangeN` is `min,max,step`.
- `If_T_step: False`: each `TrangeN` is `min,max,count`.

Example:

```yaml
If_T_step: True
nTrange: 1
Trange1: 10,50,5      # 10, 15, ..., 50 K
```

### 6.5 Experimental overlay

With `If_Data: True`, provide:

```yaml
ndata: 1
Tr1: [ ... ]          # 1/T
logtime1: [ ... ]     # log10(τ/s)
label1: Exp
```

---

## 7. File formats

### 7.1 `CFPs.out`

```
The Perturbation CFPs
<blank line>
Total modes 201
Total points 10
Mode 1 Point 1
<5 numbers>
<9 numbers>
<13 numbers>
...
```

The blank line is optional in recent versions (the reader searches for `Total modes`/`Total points`), but `Extract_CFPs.py` writes it.

### 7.2 `freq.log`

VASP frequency output with blocks:

```
   1 f  =  102.701176 THz ... 3425.742380 cm-1 ...
             X         Y         Z           dx          dy          dz
      ...
```

The mode order must match `CFPs.out`.

### 7.3 `eqpyout`

PySCF output containing the `Rotated Crystal Field Parameters:` table and spin-orbit eigenstates. For Molcas/OpenMolcas use `Ab_Type=0/1`.

### 7.4 `.eso`

`TO_foldname` contains `<Nstate>_real.eso` and `<Nstate>_imag.eso`:

```
k = 2,q = -2
<16 × 16 real numbers>
...
```

The shipped `examples/relaxation/ESO_TTO/` files are for 16 states (Dy³⁺).

### 7.5 Templates

`SPC_Tool/inpfile/*.tmpl` are **not runnable**. They contain three placeholder lines:

```
Molecular_Name
Molecular_Charge
Calculate_Title
```

`insert_xyz.py` replaces these lines. Complete the basis block before running.

---

## 8. Output files

### 8.1 `Out_foldname`

| File | Content |
|---|---|
| `OB_v1.dat` | columns `1/T`, `-log10(τ/s)`, standard deviation |
| `RM_v1.dat` | same for the Raman channel |
| `TPs1.png` | `log10(τ/s)` vs `1/T` |
| `TPs2.png` | comparison figure |
| `uniform.dat` | uniformization results (`1/T`, `-log10(rate)` for three estimates) |
| `OBTrans.eps`, `RMTrans.eps`, `RMPie.eps` | mode-resolved analysis (after `Mod_Analyse_OB`/`RM`) |

The Y-axis of `TPs1.png` is `log10(τ/s)`; column 2 of `OB_v1.dat` is its negative (`-log10(τ/s)`). Negate column 2 to recover the plotted quantity.

### 8.2 `TPs_foldname`

| File | Content |
|---|---|
| `OB<T>.tps` | Orbach transition matrix at `T` K (text) |
| `OB<T>.tp` | summed Orbach matrix |
| `LM<T>.tps` | Raman transition matrix |
| `LM<T>.tp` | summed Raman matrix |

If a cached file exists, it is read instead of recomputed. Delete it to force recalculation. The `If_select_sigma`/`select_sigma` header is checked by `check_select_sigma`.

### 8.3 `VP_filename`

The spin-phonon matrix. A `.dat` file is written first and renamed on completion; an interrupted run leaves `*.dat` and can be removed safely.

### 8.4 `psh.py`

`PSH` writes `SOC.txt`, `Lx/Ly/Lz.txt`, `Sx/Sy/Sz.txt`, `eigval.txt` and `eig_*.txt` to `outdir` (default: current directory).

---

## 9. Python API

### 9.1 `CAL_RELAX`

```python
from mag_anal.Anal_Tau import CAL_RELAX

tau = CAL_RELAX("my_system.yml")
tau.kernel()
```

`kernel()` calls `Read_input`, `Read_freq`, `Read_CFPs`, `Read_Static_effH`, `Analyse_Relax`.

Useful attributes after `kernel()`:

| Attribute | Meaning |
|---|---|
| `tau.T_list` | temperature grid |
| `tau.TP` | list of `TPs_Analyze` (one per T) |
| `tau.CFPs` | `CFPs_Analyze` object |
| `tau.Eigen` | `Eigenstates` object |
| `tau.Out_foldname` | output folder |

### 9.2 `CFPs_Analyze`

```python
from mag_anal.Anal_Tau import CFPs_Analyze

cfps = CFPs_Analyze("CFPs.out", CFPs_X, atom_list="Atom_list")
cfps.propare_CFPs(mol, CFPs_Type=1, order=1)
```

Methods: `read_CFPs`, `read_CFPs_xyz`, `read_gtensors`, `propare_CFPs`, `gen_slope`, `gen_slope_xyz`, `gen_slope_reduced`, `transB_sl`, `transB_sl_select`, `write_CFPs`.

### 9.3 `Eigenstates`

Reads energies/wavefunctions from PySCF (`read_energies_Pyscf`, `read_wavefunction_Pyscf`), Molcas (`_M`) or OpenMolcas (`_OM`), builds the tensor-operator matrices and the spin-phonon matrix.

### 9.4 `TPs_Analyze`

| Method | Purpose |
|---|---|
| `read_OBTPs`, `read_RMTPs` | read cached matrices |
| `sum_OBRMTPs` | sum channels |
| `Mod_Analyse_OB`, `Mod_Analyse_RM` | rate-determining-step analysis |
| `save_TP`, `save_TPP` | write matrices |

### 9.5 `PSH` (pseudospin Hamiltonian)

```python
from mag_anal.psh import PSH

psh = PSH(mysiso, Nstate=16, lmax=12, TFO=1, TFB=0, ROT=True, outdir="psh_out")
psh.kernel()
```

| Argument | Meaning |
|---|---|
| `SISO` | `embed_sim.siso.SISO` object |
| `Nstate` | number of `|J M>` states |
| `lmax` | maximum tensor rank |
| `TFO` | `0` EasySpin ESO, `1` tesseral tensor operators (Molforge style) |
| `TFB` | `0` least squares, `1` projection |
| `ROT` | rotate CFPs to the easy axis |
| `Angle` | explicit ZYZ Euler angles; `None` uses the easy axis |
| `outdir` | directory for diagnostic text files |

`psh.kernel()` runs `gen_ESO`, `get_E_J`, `get_g_tensor` and `get_Bkq`/`get_rot_Bkq`.

### 9.6 `gen_structures` / `Select_group`

```python
from mag_anal.gen_structures import mol_freq
from mag_anal.Select_group import choose_group
```

`mol_freq(outfile, posfile)` reads VASP phonons; `choose_group(cryfile, center, outfold, r1, r2)` selects a molecular group.

---

## 10. HPC notes

The files under `SPC_Tool/sh/` are templates. Before submitting:

- replace `YOUR_PARTITION` with your partition/queue;
- replace `YOUR_EMAIL` or remove the mail line;
- set `VASP_BIN` instead of a hard-coded path;
- replace `module load ...` / `# source activate ...` with the commands that provide PySCF, LAMP-emb and ASE;
- make `num_freq`, `num_pt` and `name` consistent across `run1_stru.sh`, `run2_Regr.sh`, `run3_inp.sh`, `run4_CFPs.sh` and the job scripts.

Slurm example:

```bash
#!/bin/bash
#SBATCH -o stdout
#SBATCH --partition=YOUR_PARTITION
#SBATCH -J lamp_tau
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --time=100:00:00

module purge
module load pyscf
export OMP_NUM_THREADS=7
export PYTHONPATH=/path/to/LAMP_emb:$PYTHONPATH
python dmet_cassiso.py
```

PBS example:

```bash
#!/bin/bash
#PBS -N lamp_tau
#PBS -q YOUR_QUEUE
#PBS -l nodes=1:ppn=28
#PBS -l walltime=100:00:00
#PBS -o stdout
#PBS -j oe

cd "$PBS_O_WORKDIR"
python dmet_cassiso.py
```

The original scripts also used Tianhe `yhbatch`; replace it with the current scheduler commands.

---

## 11. Troubleshooting

| Symptom | Cause / fix |
|---|---|
| `ModuleNotFoundError: embed_sim` | Install LAMP-emb and add it to `PYTHONPATH`. |
| `ModuleNotFoundError: mag_anal` | Install with `pip install -e .` or run from the repository root. |
| `FileNotFoundError: .../Atom_list` | Set `Atom_list` in the YAML or create `inpfile/Atom_list`. |
| `NotImplementedError: Ab_Type=2` | ORCA output is not implemented; use `Ab_Type=0/1/3`. |
| `FileNotFoundError: CFPs.out` | Run `run4_CFPs.sh` first. |
| `KeyError: Missing required key` | Add the missing key to the YAML (the message lists it). |
| `ValueError: Could not find 'Total modes'` | The CFP file is not in the expected format; see §7.1. |
| Plots do not open on a server | Set `MPLBACKEND=Agg`. |
| `TPs` values look stale | Delete the cached `OB<T>.tps`/`LM<T>.tps` and rerun. |
| Shell script does nothing | The scheduler header/modules still belong to the author's cluster; edit them (§10). |
| `vasp_std: command not found` | Set `export VASP_BIN=/path/to/vasp_std`. |

---

## 12. Validation

The shipped DyAC example was used as a smoke test:

- `python -m compileall mag_anal SPC_Tool` passes;
- `python -m mag_anal.Anal_Tau --yml examples/relaxation/CFPs_v1.yml` completes;
- the 201-mode, 10-50 K run reports `Ueff ≈ 190 K` and `log(τ/s) ≈ -11.3`;
- `tests/test_smoke.py` checks that `OB_v1.dat` is produced (uses two temperatures).

Run:

```bash
MPLBACKEND=Agg python -m pytest tests -q
```

---

## 13. Citation, license and support

If you use LAMP-tau, please cite the manuscript and the archived release:

- manuscript: *Accelerating ab initio spin-phonon relaxation simulation of single-ion magnets by quantum embedding and spatial truncation*, J. Chem. Phys. (2026);
- software: see `CITATION.cff`; the Zenodo DOI is added after the first release.

License: MIT (see `LICENSE`). Third-party programs keep their own licenses:

- PySCF — Apache-2.0;
- LAMP-emb — not redistributed here; install it separately;
- VASP, OpenMolcas, ORCA — their own licenses;
- ASE, NumPy, SciPy, Matplotlib, SymPy, prettytable, PyYAML — their own licenses.

For questions, open an issue on GitHub. For scientific questions, contact the corresponding authors of the paper.
