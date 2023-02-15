from create_job import create_job, get_starting_index, print_last_task_index
import numpy as np
# #!/bin/bash

# Kir=0.0001
# N=-1
# SigmaMax=6.0
# grid=1000
# Lambdas=(
# # 100000
# # 100
# # 10
# # 5
# # 4
# # 3
# # 2
# 1
# )

# for Lambda in "${Lambdas[@]}"; do
#     for ((x=0;x<=100;x++)); do
#         mu=0`echo "var=$x.;var*=0.008;var" | bc -l`
#         for ((y=0;y<=100;y++)); do
#             T=0`echo "var=$y.;var*=0.006;var+=0.01;var" | bc -l`
#             echo $Lambda $mu $T
#             sbatch execute_flow.sh $N $Lambda $Kir $SigmaMax $grid $mu $T mean_field_test/Lambda=$Lambda
#         done
#     done
# done

Lambda = 1
Kir = 1e-4
N = -1
grid = 1000
sigmaMax = 6.0
index_offset = get_starting_index()

mu_array = np.linspace(0.0, 0.8, 30)
T_array = np.linspace(0.01, 0.7, 30)

for mu in mu_array:
    for T in T_array:
        create_job(N, Lambda, Kir, sigmaMax, grid, mu, T,
                   f'mean_field_test/Lambda={Lambda}', index_offset)
        index_offset += 1
