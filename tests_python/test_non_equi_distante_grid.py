from scipy import interpolate
import numpy as np
from python_files.gross_neveu.Gross_Neveu import GN_2p1
from python_files.grid_creator import create_inhomogenious_grid_from_cell_spacing


def test_Gross_Neveu():
    sigma_max = 200.0
    delta_sigma = 0.006
    mu = 0.0
    T = 0.01
    samples = 3
    N_Flavor = 2
    Lambda = 100
    kir = 1e-1
    grid_points = np.arange(0.0, sigma_max, delta_sigma)
    model_linear = GN_2p1(grid_points, Lambda, kir, samples,
                          mu, T, N_Flavor, h=1.0, sigma_0=1.0)
    y_linear = model_linear.return_data.solution[-1]
    x_linear = model_linear.return_data.grid
    f_linear = interpolate.interp1d(x_linear, y_linear)(x_linear)

    grid_points = create_inhomogenious_grid_from_cell_spacing(
        sigma_max, delta_sigma)
    model_non_linear = GN_2p1(grid_points, Lambda, kir, samples,
                              mu, T, N_Flavor, h=1.0, sigma_0=1.0)

    y_non_linear = model_non_linear.return_data.solution[-1]
    x_non_linear = model_non_linear.return_data.grid
    f_non_linear = interpolate.interp1d(x_non_linear, y_non_linear)(x_linear)

    np.testing.assert_allclose(f_non_linear, f_linear, atol=1.4e-2, rtol=0)
    np.testing.assert_allclose(f_non_linear, f_linear, atol=0, rtol=1e-2)
