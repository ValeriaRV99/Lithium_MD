#!/bin/bash
#SBATCH --job-name=CPU1
#SBATCH --nodes=1
#SBATCH --ntasks=10 
#SBATCH --cpus-per-task=1
#SBATCH --partition=price-pi
#SBATCH --mem=124000
#SBATCH --time=10:00:00
#SBATCH --output=slurm.%N.%j.out
#SBATCH --error=slurm.%N.%j.err

export XDG_RUNTIME_DIR=/scratch/vr371/tmp/

module use /projects/community/modulefiles
module load intel/19.1.1
module load python/3.8.5-gc563
module load mvapich2/2.2

srun -n 10 --mpi=pmi2 python3 1.py > "1".out

