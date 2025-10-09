#!/bin/bash
#
#SBATCH --job-name=ccwatm                     # Specify job name
#SBATCH --partition=shared         # Specify partition name for job execution
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:05:00                       # Set a limit on the total run time
#SBATCH --account=ch0636                     # Charge resources on this project account
#SBATCH --output=ccwatm.o
#SBATCH --error=ccwatm.e

ulimit -s unlimited
module purge

source ~/.bashrc 
conda activate mpi_env

settingsfile="settings_CCWatM_5min_example.ini"

time srun python3 run_cwatm.py "$settingsfile"
