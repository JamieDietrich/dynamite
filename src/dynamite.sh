#! /bin/bash

jid1=$(sbatch dynamite5.slurm)
A=$(echo $jid1 | cut -d' ' -f4)
echo $A
jid2=$(sbatch --dependency=afterany:$A dynamite_merge.slurm)
echo $jid2 | cut -d' ' -f4
