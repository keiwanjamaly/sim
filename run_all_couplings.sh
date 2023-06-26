#!/bin/bash

list=(10 200 500 1000)

for ((i = 14; i < 17; i++)); do
  for j in "${list[@]}"; do
      sbatch run_coupling.sh $i $j
  done
done

