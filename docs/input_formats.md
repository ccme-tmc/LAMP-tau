# Input formats

This page documents the files consumed by `mag_anal.Anal_Tau` and by the SPC workflow. All paths in the YAML file are resolved relative to the directory **from which you run the command**; the shipped example is run from the repository root.

## YAML parameter file

A complete example is `examples/relaxation/CFPs_v1.yml`.

| Key | Type | Meaning |
|---|---|---|
| `CFPs_filename` | path | Per-displacement CFP file produced by `Extract_CFPs.py` (`CFPs.out`). |
| `mag_filename` | path | Equilibrium electronic-structure output (`eqpyout`, Molcas/OpenMolcas `*.log`, …). |
| `Ab_Type` | int | `0` Molcas, `1` OpenMolcas, `3` PySCF. `2` (ORCA) is not implemented. |
| `freq_filename` | path | VASP frequency output containing `cm-1` frequencies and displacement vectors. |
| `POSCAR` | path | Equilibrium structure used for the phonon calculation. |
| `Atom_list` | path, optional | Atom indices retained by the spatial-truncation selection (default `inpfile/Atom_list`). |
| `First_list` | path, optional | First-shell atom indices for `CFPs_Type=2` (default `inpfile/First_list`). |
| `TO_foldname` | path | Folder with the tensor-operator matrices `<J1M1|O_k^q|J2M2>`. |
| `VP_filename` | path | Output/prefix for the spin-phonon matrix elements. |
| `TPs_foldname` | path | Folder used to cache the transition-probability matrices. |
| `Out_foldname` | path | Folder for τ(T) data and figures. |
| `CFPs_X` | list | Displacement grid (Å for Cartesian modes, zero-point amplitudes for normal modes). Must contain `0`. |
| `CFPs_Type` | int | `0` normal-mode CFPs, `1` Cartesian (xyz) CFPs, `2` xyz + spatial reduction, `3` Lorentz random CFPs. |
| `Correction_slope`, `Corrected_intercept` | float | Linear frequency correction `ν → a·ν + b`. |
| `Fitting_order` | int | Polynomial order for the CFP-vs-displacement fit. |
| `Linewidth_parameter` | float | Gaussian broadening of the delta function (cm⁻¹). |
| `If_Perb`, `Perb_T`, `g_z` | bool/float | Optional perturbation treatment and the ground-state `g_z`. |
| `If_Orbach`, `If_Raman` | bool | Enable the Orbach / Raman channels. |
| `If_uniformization` | bool | Use the uniformization method. |
| `If_get_picture_Tt` | bool | Produce τ(T) figures. |
| `If_Data` | bool | Overlay experimental τ(T) data (`Tr1`, `logtime1`, `label1`, …). |
| `If_T_step`, `nTrange`, `Trange1..n` | bool/int/list | Temperature grid. `If_T_step=True`: `min,max,step`; otherwise `min,max,count`. |
| `If_select_sigma`, `select_sigma` | bool/float, optional | Restrict the mode list to modes within `select_sigma` cm⁻¹ of a transition. |

## `CFPs_py.out` (`Extract_CFPs.py` output)

The first three lines are free-form; the next two are

```
Total modes 201
Total points 10
```

Then, for every mode and point:

```
Mode <i> Point <j>
<5 numbers: k=2>
<9 numbers: k=4>
<13 numbers: k=6>
```

Each number is `B_kq` at that displacement, already shifted by the equilibrium value. `mag_anal` inserts the zero-displacement point automatically.

## `freq.log`

A VASP frequency output with blocks of the form

```
   1 f  =  102.701176 THz ...  3425.742380 cm-1 ...
             X         Y         Z           dx          dy          dz
      ...
```

Only the first frequency block is used when several are present. The order of the modes must match the order used in `CFPs_py.out`.

## `eqpyout`

A PySCF output containing the `Rotated Crystal Field Parameters:` table and the spin-orbit eigenstates. Molcas/OpenMolcas outputs are selected with `Ab_Type=0/1`.

## `.eso` tensor-operator files

`TO_foldname` must contain `<Nstate>_real.eso` and `<Nstate>_imag.eso`, each starting with `k = 2,q = -2` followed by `Nstate × Nstate` real numbers. The shipped `ESO_TTO/` files are for `Nstate=16` (Dy³⁺, eight Kramers doublets).

## Templates

`SPC_Tool/inpfile/dmet_tmp.py.tmpl` and `cassi_tmp.py.tmpl` are **not runnable templates**. `insert_xyz.py` replaces the lines

```
Molecular_Name
Molecular_Charge
Calculate_Title
```

with the structure file, total charge and job title, then writes `dmet_cassiso.py`. Complete the basis block for the elements in your system.
