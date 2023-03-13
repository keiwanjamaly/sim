from Flow import Flow
import Grid
import numpy as np


def main():
    spatial_dimension = 2
    kir = 1e-4
    n_flavor = 2
    # n_flavor = np.Inf
    Lambda_array = np.geomspace(1e3, 1e4, 10)
    sigma_max = 1000
    extrapolation_oder = 1
    n_grid = 1000

    mu = 0.0
    T = 0.3
    path = f'./cutoff_test/d_{spatial_dimension}_{extrapolation_oder}'

    for Lambda in Lambda_array:
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
