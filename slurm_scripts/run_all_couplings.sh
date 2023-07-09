#!/bin/bash

# list=(10 200 500 1000)
list=(20 30 40 50)

for ((i = 11; i < 17; i++)); do
  for j in "${list[@]}"; do
	echo $i $j
      sbatch run_coupling.sh $i $j
  done
done

