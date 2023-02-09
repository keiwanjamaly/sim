#!/bin/bash

Lambda=100000
Kir=0.0001
T=0.0125
mu=0.6
N=2
SigmaMax=6.0
deltaSigmaMin=0.0002
deltaSigmaMax=0.1
nums=(
    0.002
    0.00245725
    0.00301904
    0.00370927
    0.0045573
    0.00559922
    0.00687934
    0.00845213
    0.0103845
    0.01275867
    0.01567562
    0.01925946
    0.02366267
    0.02907255
    0.03571927
    0.04388561
    0.05391897
    0.06624621
    0.08139178
    0.1
)
for i in "${nums[@]}"; do
    # echo "$i"
    grid=$(echo "$SigmaMax/$i" | bc)
    # echo $grid
    python3 Flow.py -N $N -L $Lambda -kir $Kir -sigmaMax $SigmaMax -grid $grid -mu $mu -T $T -o spatial_resolution_test/T=${T}_mu=$mu
done