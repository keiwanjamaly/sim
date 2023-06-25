#!/bin/bash

list=(10) # 200 500 1000)

for ((i = 2; i < 3; i++)); do
  for j in "${list[@]}"; do
      sbatch run_coupling.sh $i $j
  done
done

