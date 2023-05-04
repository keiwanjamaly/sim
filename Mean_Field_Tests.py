from Gross_Neveu import GN_1p1
from compute_observables import DataClass

import numpy as np
from joblib import Parallel, delayed


def calculate_sigma(mu, T):
    sigma_max = 6.0
    Lambda = 300
    kir = 1e-4
    N_Grid = 1000
    samples = 100

    sol = GN_1p1(sigma_max, Lambda, kir, N_Grid, samples, mu, T)
    y = sol.return_data.solution
    x = sol.return_data.grid
    time = sol.return_data.time

    data_class = DataClass(x, time, y)
    sigma_0_ir = data_class.sigma_0[-1]

    return mu, T, sigma_0_ir


def main():
    mu_array = np.linspace(0, 1.2, 100)
    T_array = np.linspace(0.01, 1.2, 100)

    job_list = []

    for mu in mu_array:
        for T in T_array:
            job_list.append(delayed(calculate_sigma)(mu, T))

    result = Parallel(n_jobs=-1)(job_list)

    print(result)



if __name__ == "__main__":
    main()


