#!/bin/bash
#SBATCH -p HyperNodes
#SBATCH -J SubPhase
#SBATCH -o output5.o
#SBATCH -e output5.o

echo "Subhalo phase space"
python run_analysis.py 5
