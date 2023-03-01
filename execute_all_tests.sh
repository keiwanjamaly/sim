#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=600
#SBATCH --time=72:00:00
#SBATCH --partition=itp
#SBATCH -o log/output_%j.txt
#SBATCH -e log/error_%j.txt

eval "$(conda shell.bash hook)"
conda activate flow-master
python3 calculate_numeric_tests.py