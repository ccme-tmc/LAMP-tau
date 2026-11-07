# Workflow

The workflow follows the four steps described in the paper.

```
phonons (VASP)  →  displaced geometries  →  DMET+CASSCF-SO / CASSCF  →  CFPs
                                                                        ↓
                                                    SPC parameters  →  τ(T)
```

## Step 1 — geometry optimisation and phonons

Use the VASP helpers in `SPC_Tool/opt_freq/vaspsh/`:

- `POSCAR.sh` / `gen_POSCAR.py` create the starting `POSCAR`;
- `opt.sh` optimises the structure;
- `freq.sh` computes the Hessian/frequencies and calls `Select_group.py` to choose the molecular group.

Outputs: `freq.log`, `POSCAR`, `eq_opt.xyz`, `dpxs`, `Atom_list`.

## Step 2 — displaced structures

```bash
cd SPC_Tool/Per_B_kp/Pyscfsh
# Cartesian displacements (used for the CFP derivatives)
bash run1_stru.sh
# If the phonon eigenvectors come from a periodic cell, extract the molecular group:
bash run2_Regr.sh
```

Edit the number of modes/points in the scripts first; see `docs/hpc.md`.

## Step 3 — electronic-structure inputs

```bash
bash run3_inp.sh
```

This fills `SPC_Tool/inpfile/dmet_tmp.py.tmpl` for every displaced structure and creates `perb/<mode>/<point>/dmet_cassiso.py`. It also copies the scheduler template.

The DMET+CASSCF-SO calculation itself uses [LAMP-emb](https://github.com/ZGuan97/LAMP_emb) (`embed_sim`) and PySCF:

```bash
python dmet_cassiso.py
```

For the OpenMolcas path, use the corresponding Molcas templates and `sh/job_*_molcas.sh`.

## Step 4 — CFP extraction

```bash
bash run4_CFPs.sh
```

This writes `CFPs.out`, which is the `CFPs_filename` expected by `mag_anal`.

## Step 5 — relaxation time

```bash
lamp-tau --yml examples/relaxation/CFPs_v1.yml
```

The YAML file must point to `CFPs.out`, the equilibrium output, `freq.log`, `POSCAR`, `Atom_list` and the tensor-operator folder. Outputs (including `OB_v1.dat`, `RM_v1.dat` and `TPs1.png`) go to `Out_foldname`.

## Script-to-step map

| Step | Scripts |
|---|---|
| phonons | `SPC_Tool/opt_freq/vaspsh/{POSCAR,opt,freq}.sh` |
| selection | `SPC_Tool/opt_freq/vaspsh/Select_group.py`, `mag_anal/Select_group.py` |
| displacements | `SPC_Tool/Per_B_kp/Pyscfsh/run1_stru.sh`, `gen_structures_xyz.py`, `gen_structures_cell.py`, `Replace_group.py` |
| inputs | `run3_inp.sh`, `insert_xyz.py`, `inpfile/*.tmpl` |
| run | `SPC_Tool/sh/job_*_pyscf.sh`, `jobeq_pyscf.sh`, `check_pyscf.sh` |
| CFPs | `run4_CFPs.sh`, `Extract_CFPs.py` |
| relaxation | `mag_anal/Anal_Tau.py` |
