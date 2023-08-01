from python_files.gross_neveu.compute_observable import sigma
from python_files.grid_creator import create_inhomogenious_grid_from_cell_spacing
from python_files.gross_neveu.Gross_Neveu import get_model
from python_files.gross_neveu.compute_observable import compute_u as compute_u_gn
from python_files.gross_neveu.mean_field_computation import compute_potential_2p1
import numpy as np


def compute_sigma(one_over_g2, dimension, mu, T, sigma_max,
                  Lambda, kir, delta_sigma, N_Flavor, h, sigma_0):
    result = sigma(one_over_g2, dimension, mu, T, sigma_max,
                   Lambda, kir, delta_sigma, N_Flavor, h, sigma_0)

    return mu, T, result


def compute_u(one_over_g2, dimension, mu, T, sigma_max,
              Lambda, kir, delta_sigma, N_Flavor, h, sigma_0):
    if np.isinf(Lambda) and np.isinf(N_Flavor) and dimension == 2:
        x = np.arange(0, sigma_0 + 1, 0.00006)
        result = compute_potential_2p1(x, mu, T, h, sigma_0)
    else:
        x, result = compute_u_gn(one_over_g2, dimension, mu, T, sigma_max,
                                 Lambda, kir, delta_sigma, N_Flavor, h, sigma_0)

    return mu, T, x, result
