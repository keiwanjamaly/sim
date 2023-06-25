#!/bin/bash

list=(10 200 500 1000)

for ((i = 0; i < 15; i++)); do
  for j in "${list[@]}"; do
      sbatch run_coupling.sh $i $j
  done
done

