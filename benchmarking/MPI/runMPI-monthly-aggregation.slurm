#!/bin/bash
#SBATCH --job-name=mpi_Job-day-level
#SBATCH --partition=high_mem
##SBATCH --nodes=5
#SBATCH --ntasks-per-node=16
#SBATCH --exclusive
#SBATCH --qos=normal+
#SBATCH --mem=320GB
#SBATCH --output=slurm-%x-%j-%u-%N.out
#SBATCH --error=slurm-%x-%j-%u-%N.out

##sbatch command example: sbatch -N2 runMPI-monthly-aggregation.slurm a.py

echo "node number: $SLURM_NNODES"
echo "nodes: $(echo $SLURM_NODELIST | cut -d'[' -f 2 | cut -d']' -f 1)"

python_file=$1
echo "mpi command: mpirun -N $((16*$SLURM_NNODES)) python3 $python_file"
mpirun -N $((16*$SLURM_NNODES)) python3 $python_file

