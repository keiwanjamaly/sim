import numpy as np
from python_files.grid_creator import create_inhomogenious_grid_from_cell_spacing
from python_files.gross_neveu.Gross_Neveu import get_model
from python_files.gross_neveu.couplings.couplings_io import get_exact_coupling_from_file
from python_files.gross_neveu.couplings.couplings_io import generate_filename
from python_files.data_class import Potential
from python_files.numeric_tests.observables_evaluation import plot_or_save_observables
import matplotlib.pyplot as plt


def kir_max_test(mu, T, save=None):
    Lambda = 1000.0
    sigma_max = 2000
    delta_sigma = 0.006
    N_Flavor = 2
    samples = 1000
    kir = 1e-2
    dimension = 2
    sigma_0 = 1.0
    h = 1.0
    grid_points = create_inhomogenious_grid_from_cell_spacing(
        sigma_max, delta_sigma)
    model = get_model(dimension)
    filename = generate_filename(Lambda, N_Flavor, "./data")
    one_over_g2 = get_exact_coupling_from_file(filename)
    model = model(grid_points, Lambda, kir, samples,
                  mu, T, N_Flavor, h, one_over_g2, sigma_0)

    if model.computation_status_code == 0:
        pass
    elif model.computation_status_code == -4:
        print(
            f'Calculation exited for mu = {mu}, T = {T}. Using -1 as a result')
        return grid_points, np.zeros(len(grid_points))
    else:
        raise RuntimeError(
            f'error code {model.computation_status_code} is not handled.')

    y = model.return_data.solution
    x = model.return_data.grid
    times = model.return_data.time
    potentials = [Potential(x, u) for u in y]

    sigma_array = [potential.sigma for potential in potentials]
    mass_square_array = [potential.mass_square for potential in potentials]
    k_array = [Lambda * np.exp(-t) for t in times]
    if save is None:
        fig, (axsigma, axmasssquare) = plt.subplots(1, 2)
        axsigma.plot(k_array, sigma_array)
        axsigma.set_xscale('log')
        axmasssquare.plot(k_array, mass_square_array)
        axmasssquare.set_xscale('log')
        axmasssquare.set_yscale('log')
        plt.show()
    else:
        print('k', 'sigma', 'masssquare', sep='\t', file=save)
        for k, sigma, mass_square in zip(k_array, sigma_array, mass_square_array):
            print(k, sigma, mass_square, sep='\t', file=save)
