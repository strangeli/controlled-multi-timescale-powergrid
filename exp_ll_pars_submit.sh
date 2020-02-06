#!/bin/bash

#SBATCH --qos=medium
#SBATCH --job-name=exp_ll_pars
#SBATCH --account=coen
#SBATCH --output=exp_ll_pars-%j-%N.out
#SBATCH --error=exp_ll_pars-%j-%N.err
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
julia exp_ll_pars.jl $SLURM_NTASKS
