import numpy as np
from joblib import Parallel, delayed
from scipy.optimize import newton
import Gross_Neveu
from compute_observables import DataClass
import matplotlib.pyplot as plt


def calculate_sigma(one_over_g2: float, Model, sigma_max, Lambda, kir, N_Grid,
                    N_Flavor, h, sigma_0):
    mu = 0.0
    T = 0.01
    samples = 3
    model = Model(sigma_max, Lambda, kir, N_Grid, samples,
                  mu, T, N_Flavor, h, one_over_g2, sigma_0)

    y = model.return_data.solution
    x = model.return_data.grid
    time = model.return_data.time

    data_class = DataClass(x, time, y)
    sigma_0_ir = data_class.sigma_0[-1]

    result = sigma_0_ir - sigma_0
    # plt.plot(x, y[-1])
    # plt.show()
    print(result)
    return result


def calculate_parameter(N_Flavor):
    Lambda = 100
    kir = 1e-2
    h = 1
    sigma_0 = 1

    # spatial dimension = 1
    model = Gross_Neveu.GN_2p1
    sigma_max = 12.0
    N_Grid = 2000
    one_over_g2 = model.calculate_one_g2(h, sigma_0, Lambda)
    result = newton(calculate_sigma, one_over_g2, args=(
        model, sigma_max, Lambda, kir, N_Grid, N_Flavor, h, sigma_0))
    print(f'N = {N_Flavor}, 1/g^2 = {result}')
    return (N_Flavor, result)


def main():
    N_Flavor_list = range(2, 16)

    job_list = []
    for N_Flavor in N_Flavor_list:
        job_list.append(delayed(calculate_parameter)(N_Flavor))

    result = Parallel(n_jobs=-1)(job_list)
    result = np.array(result)
    print(result)


if __name__ == "__main__":
    main()
