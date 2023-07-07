from multiprocessing import Pool
import numpy as np
import h5py
import argparse
from python_files.gross_neveu.couplings.couplings_io import get_exact_coupling_from_file
from python_files.gross_neveu.couplings.couplings_io import generate_filename
from python_files.phase_diagram.computation_function import compute_simga


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
    dimension = 2
    sigma_max = 2000
    kir = 1e-2
    delta_sigma = 0.006
    h = 1.0
    sigma_0 = 1.0

    job_list = []

    for mu in mu_array:
        for T in T_array:
            job_list.append([mu, T])

    with Pool() as p:
        result_async = [p.apply_async(compute_simga, (one_over_g2, dimension, mu, T, sigma_max,
                                                      Lambda, kir, delta_sigma, N_Flavor, h, sigma_0)) for (mu, T) in job_list]
        result = [async_obj.get() for async_obj in result_async]

    combined_result = []
    for i in range(len(job_list)):
        mu, T = job_list[i]
        combined_result.append([mu, T, result[i]])

    with h5py.File(f'./data/phase_diagram/phase_diagram_Lambda_{Lambda}_N_Flavor_{N_Flavor}.hdf5', "w") as f:
        f.attrs["coupling"] = one_over_g2
        f.create_dataset("phase_diagram", data=combined_result)


if __name__ == "__main__":
    main()
