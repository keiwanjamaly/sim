#!/bin/bash
#SBATCH --job-name=foo
#SBATCH --partition=general1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1024   
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=FAIL

cd ..
echo $1
echo $2
python3.11 -m python_files.gross_neveu.couplings -N $1 -L $2
