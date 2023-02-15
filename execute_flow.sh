#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=500
#SBATCH --time=00:02:00
#SBATCH --partition=itp
#SBATCH -o log/output_%j.txt
#SBATCH -e log/error_%j.txt

N=`sed -n 1p ./tasks/$SLURM_ARRAY_TASK_ID`
Lambda=`sed -n 2p ./tasks/$SLURM_ARRAY_TASK_ID`
kir=`sed -n 3p ./tasks/$SLURM_ARRAY_TASK_ID`
sigmaMax=`sed -n 4p ./tasks/$SLURM_ARRAY_TASK_ID`
grid=`sed -n 5p ./tasks/$SLURM_ARRAY_TASK_ID`
mu=`sed -n 6p ./tasks/$SLURM_ARRAY_TASK_ID`
T=`sed -n 7p ./tasks/$SLURM_ARRAY_TASK_ID`
outFile=`sed -n 8p ./tasks/$SLURM_ARRAY_TASK_ID`


eval "$(conda shell.bash hook)"
conda activate flow-master
python3 Flow.py -N $N -L $Lambda -kir $kir -sigmaMax $sigmaMax -grid $grid -mu $mu -T $T -o $outFile

rm -f ./tasks/$SLURM_ARRAY_TASK_ID