#!/bin/bash

#SBATCH --qos=short
#SBATCH --job-name=expIII_new
#SBATCH --account=coen
#SBATCH --output=expIII_new-%j-%N.out
#SBATCH --error=expIII_new-%j-%N.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mail-type=END
#SBATCH --mail-user=strenge@control.tu-berlin.de

echo "------------------------------------------------------------"
echo "SLURM JOB ID: $SLURM_JOBID"
echo "$SLURM_NTASKS tasks"
echo "------------------------------------------------------------"

module load julia/1.1.0
module load hpc/2015
julia expIII_new.jl $SLURM_NTASKS
