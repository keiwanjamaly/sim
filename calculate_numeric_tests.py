from joblib import Parallel, delayed
from Flow import Flow
import numpy as np
import os
from itertools import chain


def execute_flow(Lambda, kir, sigma_max, n_grid, mu, T, n_flavor, path, tolerance):
    flow = Flow(Lambda, kir, sigma_max, n_grid, mu,
                T, n_flavor, tolerance=tolerance, number_of_observables=1000)
    flow.compute()
    flow.get_observables_for_all_positions()
    flow.save(path)


def generate_calculate_sigma_max_test(T, mu):
    Lambda = 1e5
    kir = 1e-4
    n_flavor = 2
    delta_sigma = 0.006

    path = f'tests/sigma_max_test/mu={mu}_T={T}'
    if not os.path.exists(path):
        os.makedirs(path)

    sigma_max_array = np.arange(1, 12, 0.5)
    N_array = [int(elem)
               for elem in list(np.rint(sigma_max_array/delta_sigma))]
    for i, (sigma_max, grid) in enumerate(np.column_stack([sigma_max_array, N_array])):
        yield delayed(execute_flow)(Lambda, kir, sigma_max, int(grid), mu, T, n_flavor, path)


def generate_calculate_delta_sigma_test(T, mu):
    Lambda = 1e5
    kir = 1e-4
    n_flavor = 2
    sigma_max = 6.0

    path = f'tests/delta_sigma_test/mu={mu}_T={T}'
    if not os.path.exists(path):
        os.makedirs(path)

    delta_sigma_array = np.geomspace(0.002, 0.1, 20)
    N_array = [int(elem)
               for elem in list(np.rint(sigma_max/delta_sigma_array))]
    for grid in N_array:
        yield delayed(execute_flow)(Lambda, kir, sigma_max, int(grid), mu, T, n_flavor, path)


def generate_calculate_cutoff_test(T, mu):
    LambdaArray = np.logspace(2, 6, 20)
    kir = 1e-4
    n_flavor = 2
    sigma_max = 6.0
    grid = 1000

    path = f'tests/cutoff_test/mu={mu}_T={T}'
    if not os.path.exists(path):
        os.makedirs(path)

    for Lambda in LambdaArray:
        yield delayed(execute_flow)(Lambda, kir, sigma_max, int(grid), mu, T, n_flavor, path)


def generate_calculate_tolerance_test(T, mu):
    Lambda = 2000
    kir = 5e-2
    n_flavor = 2
    sigma_max = 6.0
    grid = 1000
    toleranceArray = np.geomspace(1e-8, 1e-12)

    path = f'tests/tolerance_test/mu={mu}_T={T}'
    if not os.path.exists(path):
        os.makedirs(path)

    for tolerance in toleranceArray:
        yield delayed(execute_flow)(Lambda, kir, sigma_max, int(grid), mu, T, n_flavor, path, tolerance)


if __name__ == "__main__":
    job_list = iter(())
    # # sigma_max tests
    # job_list = chain(job_list, generate_calculate_sigma_max_test(0.35, 0.0))
    # job_list = chain(job_list, generate_calculate_sigma_max_test(0.35, 0.4))
    # job_list = chain(job_list, generate_calculate_sigma_max_test(0.0125, 0.6))
    # job_list = chain(job_list, generate_calculate_sigma_max_test(0.0125, 0.0))

    # # delta_sigma tests
    # job_list = chain(job_list, generate_calculate_delta_sigma_test(0.35, 0.0))
    # job_list = chain(job_list, generate_calculate_delta_sigma_test(0.35, 0.4))
    # job_list = chain(
    #     job_list, generate_calculate_delta_sigma_test(0.0125, 0.6))
    # job_list = chain(
    #     job_list, generate_calculate_delta_sigma_test(0.0125, 0.0))

    # lambda tests
    # job_list = chain(job_list, generate_calculate_cutoff_test(0.35, 0.0))
    # job_list = chain(job_list, generate_calculate_cutoff_test(0.35, 0.4))
    # job_list = chain(
    #     job_list, generate_calculate_cutoff_test(0.0125, 0.6))
    # job_list = chain(
    #     job_list, generate_calculate_cutoff_test(0.0125, 0.0))

    # tolerance tests
    job_list = chain(
        job_list, generate_calculate_tolerance_test(0.0125, 0.0))
    Parallel(n_jobs=-1)(job_list)
