from python_files.gross_neveu.compute_observable import sigma
from python_files.grid_creator import create_inhomogenious_grid_from_cell_spacing
from python_files.gross_neveu.Gross_Neveu import get_model
from python_files.gross_neveu.compute_observable import compute_u as compute_u_gn


def compute_sigma(one_over_g2, dimension, mu, T, sigma_max,
                  Lambda, kir, delta_sigma, N_Flavor, h, sigma_0):
    result = sigma(one_over_g2, dimension, mu, T, sigma_max,
                   Lambda, kir, delta_sigma, N_Flavor, h, sigma_0)

    return mu, T, result


def compute_u(one_over_g2, dimension, mu, T, sigma_max,
              Lambda, kir, delta_sigma, N_Flavor, h, sigma_0):
    x, result = compute_u_gn(one_over_g2, dimension, mu, T, sigma_max,
                             Lambda, kir, delta_sigma, N_Flavor, h, sigma_0)

    return mu, T, x, result
