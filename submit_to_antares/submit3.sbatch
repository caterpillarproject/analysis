#!/bin/bash
#SBATCH -p HyperNodes
#SBATCH -J HostProf
#SBATCH -o output3.o
#SBATCH -e output3.o

echo "Host Profiles"
python run_analysis.py 3
