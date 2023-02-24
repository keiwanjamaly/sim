#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=600
#SBATCH --time=24:00:00
#SBATCH --partition=itp
#SBATCH -o log/output_%j.txt
#SBATCH -e log/error_%j.txt

eval "$(conda shell.bash hook)"
conda activate flow-master
time python3 Flow.py