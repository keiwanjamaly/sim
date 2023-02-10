#!/bin/bash

Lambda=100000
Kir=0.0001
T=0.35
mu=0.0
N=2
deltaSigma=0.006

for ((i=2;i<=23;i++)); do
    z=`echo "var=$i.;var/=2;var" | bc -l`
    grid=`echo "var=$z;var/=$deltaSigma;var" | bc`
    sbatch execute_flow.sh $N $Lambda $Kir $z $grid $mu $T spatial_domain_test/T=${T}_mu=$mu
done