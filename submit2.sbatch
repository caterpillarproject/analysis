#!/bin/bash
#SBATCH -p HyperNodes
#SBATCH -J SHMF
#SBATCH -o output2.o
#SBATCH -e output2.o

echo "SHMF"
python run_analysis.py 2
