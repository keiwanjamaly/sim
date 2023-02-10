#!/bin/bash

Lambda=100000
Kir=0.0001
T=0.0125
mu=0.6
N=2
deltaSigma=0.006

for ((i=2;i<=26;i++)); do
    z=`echo "var=$i.;var/=2;var" | bc -l`
    grid=`echo "var=$z;var/=$deltaSigma;var" | bc`
    sbatch execute_flow.sh $N $Lambda $Kir $z $grid $mu $T spatial_domain_test
done