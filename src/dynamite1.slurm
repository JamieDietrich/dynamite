#!/bin/bash
#SBATCH --job-name=DYNAMITE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=94
#SBATCH --mem=4gb                    
#SBATCH --time=01:00:00   
#SBATCH --partition=standard
#SBATCH --account=apai   
 
# SLURM Inherits your environment. cd $SLURM_SUBMIT_DIR not needed
python dynamite.py
