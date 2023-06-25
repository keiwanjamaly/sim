#!/bin/bash
#SBATCH --job-name=foo
#SBATCH --partition=general1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=512   
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=FAIL

echo $1
echo $2
python3.11 compute_couplings.py -N $1 -L $2
