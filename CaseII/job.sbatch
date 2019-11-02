#!/bin/bash
## This job file is specific to Cartesius (the Dutch supercomputer cluster)
#SBATCH -n 256
#SBATCH --job-name=X.5W1.5
#SBATCH -t 96:00:00
#SBATCH -p normal
##SBATCH -t 1:00:00
##SBATCH -p short

ulimit -s 10240
module load openmpi
mpirun -np 256 ./dropOnDropImpact
