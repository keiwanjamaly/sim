#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=500
#SBATCH --time=4:00:00
#SBATCH --partition=itp
#SBATCH -o log/output_%j.txt
#SBATCH -e log/error_%j.txt

eval "$(conda shell.bash hook)"
conda activate flow-master
python3 calculate_phase_diagram.py