from create_job import create_job, get_starting_index, print_last_task_index
import numpy as np


def run_with(T, mu):
    LambdaArray = np.geomspace(1, 1e6, 20)
    Kir = 0.0001
    N = 2
    grid = 1000
    sigmaMax = 6.0
    index_offset = get_starting_index()

    for i, Lambda in enumerate(LambdaArray):
        create_job(N, Lambda, Kir, sigmaMax, grid, mu, T,
                   f'cutoff_test/T={T}_mu={mu}', i+index_offset)


run_with(0.35, 0.0)
run_with(0.35, 0.4)
run_with(0.0125, 0.6)
run_with(0.0125, 0.0)
