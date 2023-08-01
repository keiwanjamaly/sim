#!/bin/bash
#SBATCH --job-name=foo
#SBATCH --partition=general1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1024   
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL

cd ..
echo $1
echo $2
echo $3
python3.11 -m python_files.phase_diagram -N $1 -L $2 --boundary --save ${3}/${1}_${2}_phase_boundary.txt --save_LP ${3}/lp/lp_${1}_${2}.txt
