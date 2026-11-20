# Data

## Shipped example

`examples/relaxation/inpfile/` contains the data for the DyAC complex used in the paper:

| File | Content |
|---|---|
| `CFPs_py.out` | Per-displacement CFPs, 201 Cartesian modes × 10 points. |
| `eqpyout` | Equilibrium PySCF/CASSCF-SO output. |
| `freq.log` | VASP phonon frequencies and displacement vectors. |
| `POSCAR` | Equilibrium structure. |
| `Atom_list` | Atoms retained in the spatial-truncation example. |
| `dpxs` | Lattice-displacement indices (for the periodic-cell path). |
| `eq_opt.xyz` | Equilibrium molecular geometry. |

`examples/relaxation/inpfile/TTO/` contains the 16-state tesseral tensor-operator matrices.

The data were copied from the calculation archive of the paper. Author-specific paths in `eqpyout` were replaced by `<path>`; no numerical data were changed.

## Full data set

The complete archive (all three compounds, all displaced structures, phonons, CFP/SPC files, raw τ(T) and figure scripts) is published as a separate repository:

> **<https://github.com/ccme-tmc/SMMRT-DMET-ST-data>** — 2345 files, 297.8 MiB, tag `v1.0.0`.

| Directory there | System | Runnable configuration |
|---|---|---|
| `Dyacac3_cell_Analyse/` | DyAC, periodic cell | `inpfile/CFPs_v1.yml` |
| `DybbBr/` | DybBr | `inpfile/CFPs_Re_test2.yml` |
| `DyCp_Orcapyscf/` | DyCp, ORCA vs. PySCF | `inpfile/4f_ratecaulcu_input.yml` |
| `Spacetruncation/DyAC/` | DyAC, truncation series | `inpfile/CFPs_Re_test2.yml` |
| `Spacetruncation/DybBr/` | DybBr, truncation series | `inpfile/CFPs_Re_test2.yml` |
| `Spacetruncation/DyCp/` | DyCp, truncation series | `inpfile/4f_ratecaulcu_input.yml` |

Two configurations in `Spacetruncation/DyAC/inpfile/` (`CFPs.yml`, `CFPs_Re.yml`) read Molcas equilibrium logs that the authors no longer keep; they are archived unchanged as a record of the protocol and need `mag_filename: inpfile/eqpyout` (or a fresh Molcas run) to be usable.

The data repository also carries:

- `manifest.csv` — `path,system,kind,bytes,sha256` for every archived file;
- `tools/build_archive.py` (`--check` verifies a clone, `--source <dir>` rebuilds the archive from an origin folder) and `tools/rules.yaml` (the exclusion rules);
- `EXCLUDED.md` — every excluded group with its size and reason (scheduler logs, `*.chk` checkpoints, a stray PyMOL wheel, duplicate `VASP.log`, legacy `VP*.log` Molcas dumps, and four unused Molcas equilibrium logs).

All archived files are byte-identical copies of the originals; author cluster paths were **not** rewritten there (unlike the sanitized `examples/relaxation/inpfile` copy above), which keeps the checksums meaningful. Edit those paths before re-running on another machine.

A Zenodo DOI is optional; if one is created, add it here and in `CITATION.cff`.

## Large files

If you add new data, prefer a Zenodo/figshare DOI over committing large binaries to Git. The `.gitignore` excludes common workflow outputs (`perb/`, `TPs/`, `TPs_file/`, `*.chk`, `*.h5`). The data repository follows the same policy: nothing above 50 MiB, enforced by `tools/rules.yaml`.
