from scipy import interpolate
import numpy as np
from python_files.gross_neveu.Gross_Neveu import GN_2p1
from python_files.grid_creator import create_inhomogenious_grid_from_cell_spacing
import matplotlib.pyplot as plt
from timeit import default_timer as timer


def test_inhomogenious_grid(mu, T, save):
    sigma_max = 2000.0
    delta_sigma = 0.006
    samples = 3
    N_Flavor = 2
    Lambda = 1000
    kir = 1e-2
    grid_points = np.arange(0.0, sigma_max, delta_sigma)
    grid_points = np.append(grid_points, sigma_max + delta_sigma)
    print(f'grid points linear are {len(grid_points)}')
    start = timer()
    model_linear = GN_2p1(grid_points, Lambda, kir, samples,
                          mu, T, N_Flavor, h=1.0, sigma_0=1.0)
    print(timer() - start)
    y_linear = model_linear.return_data.solution[-1]
    x_linear = model_linear.return_data.grid

    grid_points = create_inhomogenious_grid_from_cell_spacing(
        sigma_max, delta_sigma)
    print(f'grid points non-linear are {len(grid_points)}')
    start = timer()
    model_non_linear = GN_2p1(grid_points, Lambda, kir, samples,
                              mu, T, N_Flavor, h=1.0, sigma_0=1.0)
    print(timer() - start)

    y_non_linear = model_non_linear.return_data.solution[-1]
    x_non_linear = model_non_linear.return_data.grid

    f_linear = interpolate.interp1d(x_linear, y_linear)(x_non_linear)
    f_non_linear = interpolate.interp1d(
        x_non_linear, y_non_linear)(x_non_linear)

    error_abs = np.abs(f_non_linear - f_linear)
    error_rel = np.abs(f_non_linear - f_linear)/np.abs(f_linear)

    error_abs = error_abs[1:-3]
    error_rel = error_rel[1:-3]
    x_non_linear = x_non_linear[1:-3]

    if save is None:
        plt.plot(x_non_linear, error_abs, label='abs error')
        plt.plot(x_non_linear, error_rel, label='rel error')
        plt.xscale('log')
        plt.yscale('log')
        plt.legend()
        plt.show()
    else:
        for sigma, rel, abs in zip(x_non_linear, error_rel, error_abs):
            print(sigma, rel, abs, sep='\t', file=save)

    # np.testing.assert_allclose(f_non_non_linear, f_linear, atol=1.4e-2, rtol=0)
    # np.testing.assert_allclose(f_non_non_linear, f_linear, atol=0, rtol=1e-2)
