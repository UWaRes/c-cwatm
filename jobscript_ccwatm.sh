#!/bin/bash
#
#SBATCH --job-name=ccwatm                     # Specify job name
#SBATCH --partition=shared          # Specify partition name for job execution
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=00:09:00                       # Set a limit on the total run time
#SBATCH --account=ch0636                     # Charge resources on this project account
#SBATCH --output=ccwatm.o
#SBATCH --error=ccwatm.e

ulimit -s unlimited
module purge

# Load OpenMPI first (same one used when building mpi4py)
module load openmpi/4.1.2-intel-2021.5.0

# Activate your conda environment
source ~/.bashrc  # or source /path/to/conda/etc/profile.d/conda.sh
conda activate mpi_env
# pyoasis variables
source /home/g/g300116/test_oasis/oasis3-mct/INSTALL_OASIS.levante/python/init.sh

#time mpirun --mca opal_common_ucx_opal_mem_hooks 1 -np $nproc_exe1 $exe1 : -np $nproc_exe2 $exe2
#time mpirun -np 1 python3 run_cwatm.py settings_CCWatM_5min_example.ini

time mpirun --mca opal_common_ucx_opal_mem_hooks 1 -np 1 python3 run_cwatm.py settings_CCWatM_5min_example.ini : -np 1 python3 run_oasis_dummy.py

#time mpirun -np 1 python3 run_cwatm.py settings_CCWatM_5min_example.ini
#time mpirun -np 1 python3 run_oasis_dummy.py
