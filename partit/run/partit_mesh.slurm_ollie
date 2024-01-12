#!/bin/sh
#SBATCH --nodes=1 
#SBATCH --job-name=partit
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH -p smp
#SBATCH --mem=30G
###SBATCH -w prod-0199
#SBATCH -o slurm-out.out
#SBATCH -e slurm-err.out


module purge
module load intel.compiler
module load intel.mpi
module load netcdf
#module load centoslibs
module list

export OMP_NUM_THREADS=1

EXE=../bin/fvom_ini.x

date
srun ${EXE} 
date

