from create_job import create_job, get_starting_index, print_last_task_index
import numpy as np

Lambda = 1e5
Kir = 1e-4
N = 2
grid = 1000
sigmaMax = 6.0
index_offset = get_starting_index()

mu_array = np.linspace(0.0, 0.8, 30)
T_array = np.linspace(0.01, 0.7, 30)

for mu in mu_array:
    for T in T_array:
        create_job(N, Lambda, Kir, sigmaMax, grid, mu, T,
                   f'phase_diagram', index_offset)
        index_offset += 1
