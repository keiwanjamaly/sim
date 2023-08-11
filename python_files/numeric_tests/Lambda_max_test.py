import numpy as np
from python_files.gross_neveu.Gross_Neveu import get_model
from python_files.gross_neveu.couplings.couplings_io import get_exact_coupling_from_file
from python_files.gross_neveu.couplings.couplings_io import generate_filename
from python_files.phase_diagram.computation_function import compute_u
from python_files.data_class import Potential
from python_files.numeric_tests.observables_evaluation import plot_or_save_observables


def compute_point(Lambda, N_Flavor, mu, T):
    dimension = 2
    delta_sigma = 0.006
    if np.isinf(N_Flavor):
        one_over_g2 = get_model(dimension).calculate_one_g2(
            h=1.0, sigma_0=1.0, Lambda=Lambda)
    else:
        filename = generate_filename(Lambda, N_Flavor, "./data")
        one_over_g2 = get_exact_coupling_from_file(filename)
    sigma_max = 2000
    kir = 1e-2
    h = 1.0
    sigma_0 = 1.0
    mu, T, x, result = compute_u(one_over_g2, dimension, mu, T, sigma_max,
                                 Lambda, kir, delta_sigma, N_Flavor, h, sigma_0)
    potential = Potential(x, result)
    return potential.sigma, potential.mass_square


def Lambda_max_test(mu, T, save):
    Lambda_array = np.flip(np.array(
        [10.0, 20.0, 30.0, 40.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 1500.0, 2000.0]))
    N_Flavor = 2
    solution_array_mass_square = []
    solution_array_sigma_0 = []
    for Lambda in Lambda_array:
        print(f'computing with Lambda {Lambda}')
        sigma_0, mass_square = compute_point(
            Lambda, N_Flavor, mu, T)
        solution_array_sigma_0.append(sigma_0)
        solution_array_mass_square.append(mass_square)
    plot_or_save_observables(
        Lambda_array, solution_array_mass_square, solution_array_sigma_0, save)
