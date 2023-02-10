#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=500
#SBATCH --time=00:02:00
#SBATCH --job-name=test
#SBATCH --partition=itp
#SBATCH -o log/output_%j.txt
#SBATCH -e log/error_%j.txt

eval "$(conda shell.bash hook)"
conda activate flow-master
python3 Flow.py -N $1 -L $2 -kir $3 -sigmaMax $4 -grid $5 -mu $6 -T $7 -o $8