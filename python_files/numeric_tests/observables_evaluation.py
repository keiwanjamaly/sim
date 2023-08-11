import numpy as np
import matplotlib.pyplot as plt


def plot_or_save_observables(x_range, solution_mass_square, solution_sigma_0, save=None):
    solution_array_mass_square = np.array(solution_mass_square)
    solution_array_sigma_0 = np.array(solution_sigma_0)
    best_mass_square = solution_array_mass_square[0]
    best_sigma_0 = solution_array_sigma_0[0]
    error_sigma_0 = np.abs(
        solution_array_sigma_0[1:] - best_sigma_0)/best_sigma_0
    error_mass_square = np.abs(
        solution_array_mass_square[1:] - best_mass_square)/best_mass_square

    if save is None:
        fig, ((ax_sigma, ax_sigma_error), (ax_mass_square,
              ax_mass_square_error)) = plt.subplots(2, 2)
        ax_sigma.plot(x_range, solution_array_sigma_0, 'ro')
        ax_sigma.set_xscale('log')
        ax_mass_square.plot(x_range, solution_array_mass_square, 'ro')
        ax_mass_square.set_xscale('log')
        if best_sigma_0 != 0:
            ax_sigma_error.plot(x_range[1:], error_sigma_0, 'bo')
            ax_sigma_error.set_xscale('log')
            ax_sigma_error.set_yscale('log')
        ax_mass_square_error.plot(x_range[1:], error_mass_square, 'bo')
        ax_mass_square_error.set_xscale('log')
        ax_mass_square_error.set_yscale('log')
        plt.show()
    else:
        print('xval', 'sigma0', 'errorsigma0', 'masssquare',
              'errormasssquare', sep='\t', end='\n', file=save)
        for i, x_val in enumerate(x_range[1:]):
            print(x_val, solution_array_sigma_0[i],
                  error_sigma_0[i], solution_array_mass_square[i], error_mass_square[i], sep='\t', end='\n', file=save)
