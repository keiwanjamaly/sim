from create_job import create_job, get_starting_index, print_last_task_index
import numpy as np


def run_with(T, mu):
    Lambda = 100000
    Kir = 0.0001
    N = 2
    deltaSigma = 0.006
    sigmaMaxArray = np.arange(1, 22, 0.5)
    N_array = list(np.rint(sigmaMaxArray/deltaSigma))
    index_offset = get_starting_index()

    for i, (sigmaMax, grid) in enumerate(np.column_stack([sigmaMaxArray, N_array])):
        create_job(N, Lambda, Kir, sigmaMax, int(grid), mu, T,
                   f'spatial_domain_test/T={T}_mu={mu}', i+index_offset)


run_with(0.35, 0.0)
run_with(0.35, 0.4)
run_with(0.0125, 0.6)
run_with(0.0125, 0.0)

print_last_task_index()
