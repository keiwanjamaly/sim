import numpy as np


def create_inhomogenious_grid_from_number_of_cells(sigma_max: float, N: int):
    pivot_point = 1.0
    balance_number = 2
    N_first_array = N-balance_number
    N_second_array = balance_number
    # spacing_first_array = pivot_point/N
    # first_array = np.linspace(0.0, pivot_point, N, endpoint=False)
    second_array = np.geomspace(pivot_point, sigma_max, N_second_array)
    while second_array[1] - second_array[0] > pivot_point/N_first_array:
        balance_number += 1
        N_first_array = N-balance_number
        N_second_array = balance_number
        second_array = np.geomspace(pivot_point, sigma_max, N_second_array)
    balance_number -= 1
    first_array = np.linspace(0.0, pivot_point, N_first_array, endpoint=False)
    second_array = np.geomspace(pivot_point, sigma_max, N_second_array)
    return np.concatenate((first_array, second_array))


def create_homogenious_grid(sigma_max: float, N: int):
    return np.linspace(0.0, sigma_max, N)


def create_inhomogenious_grid_from_cell_spacing(sigma_max: float, delta_sigma: float):
    pivot_point = 1.0
    first_array = np.arange(0.0, pivot_point, delta_sigma)
    N_second_array = 2
    second_array = np.geomspace(
        first_array[-1] + delta_sigma, sigma_max, N_second_array)
    while second_array[1] - second_array[0] > delta_sigma:
        N_second_array += 1
        second_array = np.geomspace(
            first_array[-1] + delta_sigma, sigma_max, N_second_array)
    N_second_array -= 1
    second_array = np.geomspace(
        first_array[-1] + delta_sigma, sigma_max, N_second_array)
    return np.concatenate((first_array, second_array))
