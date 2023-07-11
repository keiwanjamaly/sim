import numpy as np
from python_files.gross_neveu.Gross_Neveu import get_model
from python_files.gross_neveu.couplings.couplings_io import get_exact_coupling_from_file
from python_files.gross_neveu.couplings.couplings_io import generate_filename
from python_files.phase_diagram.computation_function import compute_sigma_spread, compute_sigma
import h5py
from joblib import Parallel, delayed


def phase_diagram_computationa(Lambda, N_Flavor, phase_diagram_resolution):
    mu_max = 1.2
    T_max = 1.0

    mu_array = np.arange(0, mu_max, phase_diagram_resolution)
    T_array = np.arange(0.01, T_max, phase_diagram_resolution)

    print(
        f'computing phase diagram for N_Flavor = {N_Flavor} and Lambda = {Lambda}.')
    print(
        f'with {len(T_array)} points in T direction and {len(mu_array)} points in mu direction')
    print(f'total number of points is {len(T_array) * len(mu_array)}')

    if np.isinf(N_Flavor):
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
    # manager = Manager()
    # lock = manager.Lock()

    for mu in mu_array:
        for T in T_array:
            job_list.append([one_over_g2, dimension, mu, T, sigma_max,
                             Lambda, kir, delta_sigma, N_Flavor, h, sigma_0])

    result = Parallel(n_jobs=-1, verbose=10)(
        delayed(compute_sigma_spread)(x) for x in job_list)
    # with Pool() as p:
    #     # result = p.map(compute_sigma_spread, job_list)
    #     multiple_results = [p.apply_async(
    #         compute_sigma, job) for job in job_list]
    #     result = [res.get() for res in multiple_results]

    with h5py.File(f'./data/phase_diagram/phase_diagram_Lambda_{Lambda}_N_Flavor_{N_Flavor}.hdf5', "w") as f:
        f.attrs["coupling"] = one_over_g2
        f.create_dataset("phase_diagram", data=result)
