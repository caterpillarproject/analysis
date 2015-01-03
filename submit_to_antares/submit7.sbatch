#!/bin/bash
#SBATCH -p HyperNodes
#SBATCH -J SubRad
#SBATCH -o output7.o
#SBATCH -e output7.o

echo "SubRadial"
python run_analysis.py 7
