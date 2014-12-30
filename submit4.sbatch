#!/bin/bash
#SBATCH -p HyperNodes
#SBATCH -J MassAccr
#SBATCH -o output4.o
#SBATCH -e output4.o

echo "Mass Accretion History"
python run_analysis.py 4
