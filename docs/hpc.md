# HPC notes

The files under `SPC_Tool/sh/` are **templates**. They were written for a Slurm cluster and for Tianhe (`yhbatch`); they will not run unchanged on another machine. Edit the scheduler header and the environment setup before submitting.

## Checklist

- **Partition/queue:** replace `YOUR_PARTITION` with the partition of your cluster.
- **E-mail:** replace `YOUR_EMAIL` or remove the `--mail-user` line.
- **VASP binary:** set `VASP_BIN` (for example `export VASP_BIN=/path/to/vasp_std`) instead of hard-coding it.
- **Modules/conda:** replace `module load ...` / `# source activate ...` with the commands that provide PySCF, LAMP-emb and ASE on your machine.
- **Number of modes/points:** `num_freq`, `num_pt` and `name` must agree across `run1_stru.sh`, `run2_Regr.sh`, `run3_inp.sh`, `run4_CFPs.sh`, `check_pyscf.sh` and the job scripts.

## Slurm template

```bash
#!/bin/bash
#SBATCH -o stdout
#SBATCH --partition=YOUR_PARTITION
#SBATCH -J lamp_tau
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --time=100:00:00

module purge
module load pyscf           # or: source activate pyscf
export OMP_NUM_THREADS=7
export PYTHONPATH=/path/to/LAMP_emb:$PYTHONPATH

python dmet_cassiso.py
```

## PBS template

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

## Tianhe (`yhbatch`)

The original scripts used `yhbatch` and a `TH_LONG` partition. Replace these with the current Tianhe scheduler commands and the modules available on your account.

## VASP

`SPC_Tool/opt_freq/vaspsh/` expects `vaspkit` for `POTCAR`/`KPOINTS`. Set `VASP_BIN` in `SPC_Tool/sh/job_vasp.sh`; the script uses `${VASP_BIN:-vasp_std}`.

## Headless operation

Set `MPLBACKEND=Agg` before running `lamp-tau` on a login node or in a batch job; the code does not open interactive windows by default.
