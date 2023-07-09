from multiprocessing import Pool, Manager
import numpy as np
import h5py
import argparse
from python_files.gross_neveu.couplings.couplings_io import get_exact_coupling_from_file
from python_files.gross_neveu.couplings.couplings_io import generate_filename
from python_files.phase_diagram.computation_function import compute_sigma_spread
from python_files.gross_neveu.Gross_Neveu import get_model


def main():
    mu_max = 1.2
    T_max = 1.0

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-N', type=int, help='Set the number of flavors', required=True, default=None)
    parser.add_argument('-L', type=float,
                        help='Set the UV Cutoff', required=True, default=None)
    parser.add_argument('--delta', type=float,
                        help='sets the resolution for the phase diagram', default=0.01)

    args = parser.parse_args()
    phase_diagram_resolution = args.delta
    mu_array = np.arange(0, mu_max, phase_diagram_resolution)
    T_array = np.arange(0.01, T_max, phase_diagram_resolution)

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

    if args.N == -1:
        one_over_g2 = get_model(2).calculate_one_g2(
            h=1.0, sigma_0=1.0, Lambda=Lambda)
    else:
        filename = generate_filename(Lambda, N_Flavor, "./data")
        one_over_g2 = get_exact_coupling_from_file(filename)
    dimension = 2
    sigma_max = 2000
    kir = 1e-2
    delta_sigma = 0.006
    h = 1.0
    sigma_0 = 1.0

    job_list = []
    manager = Manager()
    lock = manager.Lock()

    for mu in mu_array:
        for T in T_array:
            job_list.append([one_over_g2, dimension, mu, T, sigma_max,
                             Lambda, kir, delta_sigma, N_Flavor, h, sigma_0, lock])

    with Pool() as p:
        result = p.map(compute_sigma_spread, job_list)

    with h5py.File(f'./data/phase_diagram/phase_diagram_Lambda_{Lambda}_N_Flavor_{N_Flavor}.hdf5', "w") as f:
        f.attrs["coupling"] = one_over_g2
        f.create_dataset("phase_diagram", data=result)


if __name__ == "__main__":
    main()
