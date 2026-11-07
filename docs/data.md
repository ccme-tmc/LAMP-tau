# Data

## Shipped example

`examples/relaxation/DyAC/` contains the data for the DyAC complex used in the paper:

| File | Content |
|---|---|
| `CFPs_py.out` | Per-displacement CFPs, 201 Cartesian modes × 10 points. |
| `eqpyout` | Equilibrium PySCF/CASSCF-SO output. |
| `freq.log` | VASP phonon frequencies and displacement vectors. |
| `POSCAR` | Equilibrium structure. |
| `Atom_list` | Atoms retained in the spatial-truncation example. |
| `dpxs` | Lattice-displacement indices (for the periodic-cell path). |
| `eq_opt.xyz` | Equilibrium molecular geometry. |

`examples/relaxation/ESO_TTO/` contains the 16-state tesseral tensor-operator matrices.

The data were copied from the calculation archive of the paper. Author-specific paths in `eqpyout` were replaced by `<path>`; no numerical data were changed.

## Full data set

The complete archive (all three compounds, all displaced structures, phonons, CFP/SPC files, raw τ(T) and figure scripts) is intended for Zenodo. Add the DOI here and in `CITATION.cff` after the first release.

## Large files

If you add new data, prefer a Zenodo/figshare DOI over committing large binaries to Git. The `.gitignore` excludes common workflow outputs (`perb/`, `TPs/`, `TPs_file/`, `*.chk`, `*.h5`).
