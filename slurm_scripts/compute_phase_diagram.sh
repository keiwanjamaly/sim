#!/bin/bash
#SBATCH --job-name=phase_diagram
#SBATCH --partition=general1
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=512   
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=FAIL

cd ..
python3.11 compute_phase_diagram.py -N $1 -L $2
