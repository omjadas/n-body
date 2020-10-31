#!/bin/bash
#SBATCH --error=n-body.err
#SBATCH --output=n-body.out
#SBATCH --partition=snowy
#SBATCH --time=0:10:00
#SBATCH --nodes=8-8
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=32G
#SBATCH --job-name=n-body

module load gcc/8.3.0
module load openmpi/3.1.4

mpicxx -O3 -o n-body n-body.cpp
cat in.dat | mpiexec ./n-body
