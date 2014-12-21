#!/bin/bash
#SBATCH -p HyperNodes
#SBATCH -J Nvmax
#SBATCH -o output1.o
#SBATCH -e output1.o

echo "Nvmax"
python run_analysis.py 1
