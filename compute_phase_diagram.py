from python_files.gross_neveu.Gross_Neveu import GN_2p1
from python_files.data_class import DataClass
from python_files.grid_creator import create_inhomogenious_grid_from_cell_spacing
from python_files.gross_neveu.couplings.couplings_io import get_exact_coupling_from_file
from python_files.gross_neveu.couplings.couplings_io import generate_filename

import numpy as np
from joblib import Parallel, delayed
import argparse
import h5py


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


def calculate_sigma(mu, T, Lambda, N_Flavor, one_over_g2):
    sigma_max = 1000.0
    kir = 1e-2
    delta_sigma = 0.006
    samples = 5

    grid_points = create_inhomogenious_grid_from_cell_spacing(
        sigma_max, delta_sigma)
    sol = GN_2p1(grid_points, Lambda, kir, samples, mu, T,
                 N_Flavor=N_Flavor, one_over_g2=one_over_g2)
    y = sol.return_data.solution
    x = sol.return_data.grid
    time = sol.return_data.time

    data_class = DataClass(x, time, y)
    sigma_0_ir = data_class.sigma_0[-1]

    print(f'mu = {mu}, T = {T} - done')

    return mu, T, sigma_0_ir


def get_error(mu, T, Lambda):
    mu, T, sigma_0_ir = calculate_sigma(mu, T, Lambda)
    return mu, T, sigma_0_ir - analytic_2p1_MF(mu, T)


def main():
    mu_max = 1.2
    T_max = 1.0
    mu_array = np.arange(0, mu_max, 0.02)
    T_array = np.arange(0.01, T_max, 0.02)

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-N', type=int, help='Set the number of flavors', required=True, default=None)
    parser.add_argument('-L', type=float,
                        help='Set the UV Cutoff', required=True, default=None)

    args = parser.parse_args()
    if args.N == -1:
        N_Flavor = np.inf
    else:
        N_Flavor = args.N
    Lambda = args.L

    print(
        f'computing phase diagram for N_Flavor = {N_Flavor} and Lambda = {Lambda}.')
    print(
        f'with {len(T_array)} points in T direction and {len(mu_array)} points in mu direction')
    print(f'total number of points is {len(T_array) * len(mu_array)}')

    filename = generate_filename(Lambda, N_Flavor, "./data")
    one_over_g2 = get_exact_coupling_from_file(filename)

    job_list = []

    for mu in mu_array:
        for T in T_array:
            job_list.append(delayed(calculate_sigma)(
                mu, T, Lambda, N_Flavor, one_over_g2))
            # job_list.append(delayed(get_error)(mu, T, Lambda))

    result = Parallel(n_jobs=-1)(job_list)
    result = np.array(result)
    print("done")

    with h5py.File(f'./data/phase_diagram/phase_diagram_Lambda_{Lambda}_N_Flavor_{N_Flavor}.hdf5', "w") as f:
        f.attrs["coupling"] = one_over_g2
        f.create_dataset("phase_diagram", data=result)


if __name__ == "__main__":
    main()
