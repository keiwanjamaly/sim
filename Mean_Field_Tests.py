from Gross_Neveu import GN_2p1
from compute_observables import DataClass
import grid_creator

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from joblib import Parallel, delayed
import matplotlib.pyplot as plt


def analytic_2p1_MF(mu: float, T: float, h: float = 1, sigma_0: float = 1) -> float:
    intermediate = np.exp(h*sigma_0/T)/2 - 1
    if intermediate >= 1:
        mu_crit = T*np.arccosh(intermediate)
    else:
        mu_crit = -np.infty

    if mu <= mu_crit:
        return T/(h*sigma_0) * np.arccosh(np.exp(h*sigma_0/T)/2 - np.cosh(mu/T))
    else:
        return 0


def calculate_sigma(mu, T, Lambda):
    sigma_max = 6.0
    kir = 1e-4
    N_Grid = 1000
    samples = 50

    grid_points = grid_creator.create_homogenious_grid(sigma_max, N_Grid)
    sol = GN_2p1(grid_points, Lambda, kir, samples, mu, T)
    y = sol.return_data.solution
    x = sol.return_data.grid
    time = sol.return_data.time

    data_class = DataClass(x, time, y)
    sigma_0_ir = data_class.sigma_0[-1]

    return mu, T, sigma_0_ir


def get_error(mu, T, Lambda):
    mu, T, sigma_0_ir = calculate_sigma(mu, T, Lambda)
    return mu, T, sigma_0_ir - analytic_2p1_MF(mu, T)


def main():
    mu_max = 1.2
    T_max = 1.0
    n_grid = 50
    Lambda = 10
    mu_array = np.linspace(0, mu_max, n_grid)
    T_array = np.linspace(0.01, T_max, n_grid)

    job_list = []

    for mu in mu_array:
        for T in T_array:
            job_list.append(delayed(calculate_sigma)(mu, T, Lambda))
            # job_list.append(delayed(get_error)(mu, T, Lambda))

    result = Parallel(n_jobs=-1)(job_list)
    result = np.array(result)
    print("done")

    # plot the results
    points = list(zip(result[:, 0], result[:, 1]))
    values = list(result[:, 2])
    interpolation = LinearNDInterpolator(points, values)
    X, Y = np.meshgrid(np.linspace(0.0, mu_max, 1000),
                       np.linspace(0.1, T_max, 1000))
    Z = interpolation(X, Y)
    # print(job_list)
    plt.pcolormesh(X, Y, Z, shading='auto')
    plt.colorbar()
    plt.title(f'{Lambda=} with the largest error of {max(values):.2e}')
    plt.xlabel(r'$\mu$')
    plt.ylabel(r'$T$')

    plt.show()


if __name__ == "__main__":
    main()
