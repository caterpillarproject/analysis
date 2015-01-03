#!/bin/bash
#SBATCH -p HyperNodes
#SBATCH -J SubProf
#SBATCH -o output6.o
#SBATCH -e output6.o

echo "Sub Profiles"
python ../run_analysis.py 6
