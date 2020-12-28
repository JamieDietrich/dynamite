#!/bin/bash
#SBATCH --job-name=DYNAMITE
#SBATCH --nodes=5
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=94
#SBATCH --mem=4gb                    
#SBATCH --time=01:00:00   
#SBATCH --partition=standard
#SBATCH --account=apai   
 
# SLURM Inherits your environment. cd $SLURM_SUBMIT_DIR not needed
python dynamite.py dynamite_config.txt 1 &
python dynamite.py dynamite_config.txt 2 &
python dynamite.py dynamite_config.txt 3 &
python dynamite.py dynamite_config.txt 4 &
python dynamite.py dynamite_config.txt 5 &
wait

#SBATCH --job-name=DYNAMITE
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb                    
#SBATCH --time=01:00:00   
#SBATCH --partition=standard
#SBATCH --account=apai 

python dynamite_merge.py dynamite_config.txt
