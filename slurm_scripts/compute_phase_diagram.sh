#!/bin/bash
#SBATCH --job-name=phase_diagram
#SBATCH --partition=general1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=80
#SBATCH --mem-per-cpu=512   
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=FAIL

cd ..
python3.11 -m python_files.phase_diagram -N $1 -L $2
