import numpy as np
import Gross_Neveu


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


def main():
    import matplotlib.pyplot as plt
    grid_points = create_inhomogenious_grid_from_number_of_cells(200, 1000)
    plt.hist(grid_points, bins=400)
    plt.title(
        f'smallest cell size is {(grid_points[1] - grid_points[0]):.3f} with {len(grid_points)} grid points')
    plt.show()
    # pass
    grid_points = create_inhomogenious_grid_from_cell_spacing(2000, 0.006)
    plt.hist(grid_points, bins=400)
    plt.title(
        f'smallest cell size is {(grid_points[1] - grid_points[0]):.3f} with {len(grid_points)} grid points')
    plt.show()


def compare_grid_points():
    import matplotlib.pyplot as plt
    from scipy import interpolate
    sigma_max = 200.0
    delta_sigma = 0.006
    mu = 0.0
    T = 0.01
    samples = 3
    N_Flavor = 2
    Lambda = 100
    kir = 1e-1
    grid_points = np.arange(0.0, sigma_max, delta_sigma)
    model_linear = Gross_Neveu.GN_2p1(0.0, Lambda, kir, 1, samples,
                                      mu, T, N_Flavor, h=1.0, sigma_0=1.0, grid_points=grid_points)
    y_linear = model_linear.return_data.solution[-1]
    x_linear = model_linear.return_data.grid
    f_linear = interpolate.interp1d(x_linear, y_linear)

    grid_points = create_inhomogenious_grid_from_cell_spacing(
        sigma_max, delta_sigma)
    model_non_linear = Gross_Neveu.GN_2p1(0.0, Lambda, kir, 1, samples,
                                          mu, T, N_Flavor, h=1.0, sigma_0=1.0, grid_points=grid_points)

    y_non_linear = model_non_linear.return_data.solution[-1]
    x_non_linear = model_non_linear.return_data.grid
    f_non_linear = interpolate.interp1d(x_non_linear, y_non_linear)

    error = np.abs(f_linear(x_linear) - f_non_linear(x_linear))
    rel_error = error/np.abs(f_linear(x_linear))
    plt.loglog(x_linear, error, label="absolute")
    plt.loglog(x_linear, rel_error, label="absolute")
    plt.legend()
    plt.grid()
    # plt.xscale("log")
    plt.show()


if __name__ == "__main__":
    # main()
    compare_grid_points()
