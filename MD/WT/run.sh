#!/bin/bash

#SBATCH --job-name=NVTFvW
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --partition=price-pi
#SBATCH --mem=200GB
#SBATCH --exclusive
##SBATCH --constraint=hal
#SBATCH --time=3-00:00:00
#SBATCH --output=slrum.%N.%j.out
#SBATCH --error=slrum.%N.%j.err

export XDG_RUNTIME_DIR=/scratch/vr371/tmp/  ##needed for jupyter writting temporary files

module use /projects/community/modulefiles
module load intel/19.1.1
# module load python/3.8.5-gc563
module load mvapich2/2.2
export OMP_NUM_THREADS=1
source /projectsn/mp1009_1/Valeria/pyenv/lpps/bin/activate

srun --mpi=pmi2 -N 1 -n 64 python3 MD.py > "MD".out
