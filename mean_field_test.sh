#!/bin/bash

Kir=0.0001
N=-1
SigmaMax=6.0
grid=1000
Lambdas=(
# 100000
# 100
# 10
# 5
# 4
# 3
# 2
1
)

for Lambda in "${Lambdas[@]}"; do
    for ((x=0;x<=100;x++)); do
        mu=0`echo "var=$x.;var*=0.008;var" | bc -l`
        for ((y=0;y<=100;y++)); do
            T=0`echo "var=$y.;var*=0.006;var+=0.01;var" | bc -l`
            echo $Lambda $mu $T
            sbatch execute_flow.sh $N $Lambda $Kir $SigmaMax $grid $mu $T mean_field_test/Lambda=$Lambda
        done
    done
done