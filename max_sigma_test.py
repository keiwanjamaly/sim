from Flow import Flow
import numpy as np


def main():
    Lambda = 1e3
    kir = 1e-4
    n_flavor = 2
    # n_flavor = np.Inf
    sigma_max_array = np.geomspace(3, 12, 10)
    delta_sigma = 0.006

    mu = 0.0
    T = 0.3
    path = './sigma_max_test/d_2'

    for sigma_max in sigma_max_array:
        n_grid = int(sigma_max / delta_sigma)
        flow = Flow(Lambda, kir, sigma_max, n_grid, mu, T,
                    n_flavor, save_flow_flag=True, console_logging=True, number_of_observables=500, tolerance=1e-12)

        flow.compute()
        flow.get_observables_for_all_positions()
        if not (path is None):
            flow.save(path)


if __name__ == "__main__":
    main()
