#!/bin/bash
#SBATCH -p HyperNodes
#SBATCH -J SubRMass
#SBATCH -o output8.o
#SBATCH -e output8.o

echo "SubRadialMass (takes a long time)"
python run_analysis.py 8
