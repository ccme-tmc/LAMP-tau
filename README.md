# LAMP-tau

Python tools for *ab initio* spin-phonon relaxation simulation of single-ion magnets (SIMs), built on top of [PySCF](https://pyscf.org/) and intended to be used together with [LAMP-emb](https://github.com/ZGuan97/LAMP_emb).

This repository accompanies the manuscript

> **Accelerating ab initio spin-phonon relaxation simulation of single-ion magnets by quantum embedding and spatial truncation**
> Y. Deng et al., *J. Chem. Phys.* (2026)

and contains the modules used to compute crystal-field parameters (CFPs), magnetic anisotropy *g* tensors, spin-phonon coupling (SPC) parameters and spin-relaxation times τ(T).

> 中文简介：本仓库包含论文中使用的 DMET+CASSCF-SO / CFP / 自旋-声子耦合 / 弛豫时间计算脚本。`mag_anal` 负责分析，`SPC_Tool` 负责生成扰动结构、提取 CFPs 与组织 PySCF/VASP 计算。安装与最小示例见下文。

## What is included

| Component | Purpose |
|---|---|
| `mag_anal/` | Analysis of *ab initio* results: effective spin Hamiltonian, *g* tensors, CFPs, SPC parameters and τ(T). |
| `SPC_Tool/` | Scripts for the SPC workflow: VASP phonons, displaced structures, PySCF/OpenMolcas inputs, CFP extraction. |
| `examples/relaxation/` | A complete, runnable DyAC example for the relaxation analysis. |
| `docs/` | Workflow, input-format and HPC notes; complete user [manual](docs/manual.md) ([HTML](docs/manual.html)). |

## Installation

Requirements: **Python ≥ 3.10** (the code uses `match`/`case`), NumPy, SciPy, PyYAML, prettytable, Matplotlib, ASE and SymPy.

```bash
git clone https://github.com/Dengyf98/LAMP_tau.git
cd LAMP_tau
python -m pip install -e .
# or, if you only want the analysis part:
python -m pip install -r requirements.txt
```

The DMET+CASSCF-SO part additionally needs:

- [PySCF](https://pyscf.org/) (`pip install pyscf`) and `basis_set_exchange`;
- [LAMP-emb](https://github.com/ZGuan97/LAMP_emb) (provides the `embed_sim` package). Install it separately and make sure it is on `PYTHONPATH`.

The periodic-DFT part of `SPC_Tool` uses **VASP**; the OpenMolcas path uses **OpenMolcas**. Both are external programs and are not bundled here.

## Quick start: relaxation analysis (DyAC example)

From the repository root:

```bash
lamp-tau --yml examples/relaxation/CFPs_v1.yml
```

or

```bash
python -m mag_anal.Anal_Tau --yml examples/relaxation/CFPs_v1.yml
```

The example reads:

- `examples/relaxation/DyAC/CFPs_py.out` — per-displacement CFPs (201 Cartesian modes × 10 points);
- `examples/relaxation/DyAC/eqpyout` — equilibrium PySCF/CASSCF-SO output;
- `examples/relaxation/DyAC/freq.log` and `POSCAR` — VASP phonon frequencies/eigenvectors and structure;
- `examples/relaxation/DyAC/Atom_list` — atoms retained by the spatial-truncation selection;
- `examples/relaxation/ESO_TTO/` — tesseral tensor-operator matrices for the 16 states.

Outputs are written to `examples/relaxation/TPs_file/` (τ(T) data and figures). On a headless machine use `MPLBACKEND=Agg`.

A minimal programmatic call is:

```python
from mag_anal.Anal_Tau import CAL_RELAX

tau = CAL_RELAX("examples/relaxation/CFPs_v1.yml")
tau.kernel()
```

## Quick start: SPC workflow

The full workflow is described in `docs/workflow.md` and in `SPC_Tool/manual/manual.tex` (Chinese). In short:

1. **Phonons (VASP):** optimize the structure and compute the Hessian/frequencies; `SPC_Tool/opt_freq/vaspsh/` contains the helper scripts.
2. **Displaced structures:** `SPC_Tool/Per_B_kp/Pyscfsh/run1_stru.sh` generates the displacements, `run2_Regr.sh` extracts the molecular group.
3. **Inputs:** `run3_inp.sh` fills the templates in `SPC_Tool/inpfile/` and creates one computation per displacement.
4. **Electronic structure:** run the generated DMET+CASSCF-SO (PySCF) or CASSCF (OpenMolcas) jobs on your cluster.
5. **CFPs:** `run4_CFPs.sh` extracts the perturbed CFPs into `CFPs.out`.
6. **Relaxation:** run `lamp-tau` with a YAML file pointing to the CFP file, the equilibrium output, the phonon data and the tensor-operator folder.

The `.sh` files are **templates**: edit the partition/queue, the number of modes/points and the system name for your own cluster (see `docs/hpc.md`).

## Repository layout

```
LAMP_tau/
├── mag_anal/                  # analysis package
│   ├── Anal_Tau.py            # τ(T), SPC and relaxation analysis (CLI: lamp-tau)
│   ├── psh.py                 # pseudospin Hamiltonian, g tensors, CFP projection
│   ├── gen_structures.py      # displaced structures from phonon modes
│   └── Select_group.py        # molecular-group selection / spatial truncation
├── SPC_Tool/
│   ├── inpfile/               # templates and input examples
│   ├── opt_freq/vaspsh/       # VASP optimization/frequency helpers
│   ├── Per_B_kp/Pyscfsh/      # displaced structures, inputs, CFP extraction
│   ├── sh/                    # scheduler templates
│   ├── Relaxation/            # legacy relaxation example
│   └── manual/                # Chinese manual (LaTeX + PDF)
├── examples/relaxation/       # runnable DyAC example
├── docs/                      # workflow, input formats, HPC notes
├── pyproject.toml
├── requirements.txt
└── LICENSE
```

## Input formats

A complete user manual is available as [Markdown](docs/manual.md) and [HTML](docs/manual.html).

See `docs/input_formats.md` for the YAML keys and the file formats of `CFPs_py.out`, `freq.log`, `eqpyout` and the `.eso` tensor-operator files.

## Data and reproducibility

The `examples/relaxation/DyAC` data are copied from the calculations used in the paper (DyAC). The full input/output data set (structures, phonons, CFPs, SPC parameters, raw τ(T) and figure scripts) will be archived on Zenodo; the DOI will be added here and in `CITATION.cff` once the release is created.

## Supported platforms

- The Python analysis (`mag_anal`) is platform-independent and was tested on Linux and Windows.
- `SPC_Tool` targets HPC clusters: the scheduler headers are templates for SLURM/PBS/Tianhe; edit them for your machine.
- VASP and OpenMolcas are commercial/third-party codes and are not distributed here.

## Known limitations

- The `Ab_Type=2` (ORCA) reader is **not implemented**; use Molcas/OpenMolcas (`0`/`1`) or PySCF (`3`).
- `SPC_Tool` phonon support is currently VASP-only.
- The DMET implementation itself lives in [LAMP-emb](https://github.com/ZGuan97/LAMP_emb), not in this repository.
- The `.sh` scripts assume a Unix shell and a batch scheduler.

## Citation

If you use this code, please cite the paper above and the archived release. See `CITATION.cff`; the Zenodo DOI will be inserted after the first release.

## License

MIT — see `LICENSE`. Third-party programs (PySCF, LAMP-emb, VASP, OpenMolcas, ASE, EasySpin) keep their own licenses; LAMP-emb is **not** redistributed here.

## Contact

Please open an issue on GitHub for questions about the code. For the scientific results, contact the corresponding authors of the paper.
