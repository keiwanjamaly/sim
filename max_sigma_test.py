from Flow import Flow
import Grid
import numpy as np


def main():
    spatial_dimension = 2
    Lambda = 1e3
    kir = 1e-4
    n_flavor = 2
    # n_flavor = np.Inf
    sigma_max_array = np.geomspace(6, 1000, 10)
    # delta_sigma = 0.003
    extrapolation_oder = 2
    n_grid = 2000

    mu = 0.0
    T = 0.3
    path = f'./sigma_max_test_{n_grid}/d_{spatial_dimension}_{extrapolation_oder}'

    for sigma_max in sigma_max_array:
        grid = Grid.LinearLogGridWithContiniousTransition(
            sigma_max, n_grid, extrapolation_oder)
        flow = Flow(spatial_dimension, Lambda, kir, grid, mu, T,
                    n_flavor, save_flow_flag=True, console_logging=True, number_of_observables=500, tolerance=1e-12)
        print(flow)

        flow.compute()
        flow.get_observables_for_all_positions()
        if not (path is None):
            flow.save(path)


if __name__ == "__main__":
    main()
