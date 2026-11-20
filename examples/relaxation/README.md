# Relaxation example (DyAC)

This directory mirrors the layout of the paper's `SPC_Tool/Relaxation/`
folder:

```
examples/relaxation/
├── inpfile/                 # inputs for one compound
│   ├── CFPs_v1.yml          # YAML control file
│   ├── CFPs_py.out          # per-displacement CFPs
│   ├── eqpyout              # equilibrium PySCF/CASSCF-SO output
│   ├── freq.log             # VASP phonon frequencies/eigenvectors
│   ├── POSCAR               # equilibrium structure
│   ├── Atom_list            # atoms retained by the truncation selection
│   ├── dpxs                 # lattice-displacement indices
│   ├── eq_opt.xyz           # molecular geometry
│   └── TTO/                 # 16-state tensor-operator matrices
├── run.ipynb                # step-by-step notebook
├── expected/                # reference τ(T) data
├── TPs/                     # cache (created on first run)
└── TPs_file/                # outputs (created on first run)
```

The data are for DyAC (201 Cartesian modes, 10 displacement points,
±0.1 Å). The folder is self-contained: every path in `CFPs_v1.yml` is
relative to this directory.

## Run from the command line

```bash
cd examples/relaxation
lamp-tau --yml inpfile/CFPs_v1.yml
# or, without installing:
PYTHONPATH=../.. python -m mag_anal.Anal_Tau --yml inpfile/CFPs_v1.yml
```

On a headless machine:

```bash
MPLBACKEND=Agg lamp-tau --yml inpfile/CFPs_v1.yml
```

## Run from the notebook

Open `run.ipynb` with Jupyter and execute the cells in order. The notebook

1. adds the repository root to `sys.path`;
2. runs `CAL_RELAX("inpfile/CFPs_v1.yml")`;
3. calls `tau.TP[0].Mod_Analyse_OB()` for the rate-determining-step table;
4. compares the generated `TPs_file/OB_v1.dat` with `expected/OB_v1.dat`.

## Expected result

For the shipped DyAC data (10–50 K, Orbach only):

```
TestU_{eff} of Orbach is 190.1 K
log(tau) is -11.30 s
```

`TPs_file/OB_v1.dat` has three columns: `1/T`, `-log10(τ/s)`, standard
deviation. Negate column 2 to recover the `log10(τ/s)` plotted in
`TPs1.png`.

## Using your own compound

Replace the files under `inpfile/` with your own:

- `CFPs_py.out` from `SPC_Tool/Per_B_kp/Pyscfsh/Extract_CFPs.py`;
- `eqpyout` from the equilibrium DMET+CASSCF-SO/PySCF run;
- `freq.log`/`POSCAR` from VASP (or a molecular Hessian);
- `Atom_list` from the truncation selection;
- `TTO/<Nstate>_{real,imag}.eso` for the number of states used.

Then edit `inpfile/CFPs_v1.yml` (paths, `g_z`, `Linewidth_parameter`,
`CFPs_X`, temperature range).

## Notes

- `TPs/` is a cache. Delete it to force recomputation; an interrupted run
  can leave `TPs_file/VP.log.dat`, which is safe to remove.
- `If_Raman: True` enables the Raman channel. It is much more expensive
  than the Orbach-only example (the reference DyCp run writes ~90 MB
  `LM*.tps` files), so it is disabled here.
- The example uses DyAC; the paper's DyCp folder has the same layout, so
  you can swap `inpfile/` contents without changing the commands.
