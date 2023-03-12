from Grid import Adaptive
from Flow import Flow


def main():
    pass


if __name__ == "__main__":
    main()
    Lambda = 1e3
    kir = 1e-4
    n_flavor = 2
    # n_flavor = np.Inf
    n_grid = 1000
    sigma_max = 1000
    mu = 0.0
    T = 0.3
    path = './refinement'

    grid = Adaptive(sigma_max, 1e-3)
    flow = Flow(Lambda, kir, grid, mu, T,
                n_flavor, save_flow_flag=True, console_logging=True, number_of_observables=1000, tolerance=1e-12)
    print(flow)

    flow.compute()

    while grid.remesh(flow.solution.y, flow.solution.t):
        flow = Flow(Lambda, kir, grid, mu, T,
                    n_flavor, save_flow_flag=True, console_logging=True, number_of_observables=1000, tolerance=1e-12)
        print(flow)

        flow.compute()

    flow.get_observables_for_all_positions()
    flow.save(path)
