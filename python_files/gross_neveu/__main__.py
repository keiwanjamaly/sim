from python_files.gross_neveu.Gross_Neveu import get_model
from python_files.grid_creator import create_homogenious_grid, create_inhomogenious_grid_from_cell_spacing
import argparse
import numpy as np
import matplotlib.pyplot as plt
from python_files.data_class import Potential


def main():
    parser = argparse.ArgumentParser(prog='couplings calculator')
    parser.add_argument(
        '-N', type=int, help='Set the number of flavors', default=None)
    parser.add_argument('-L', type=float,
                        help='Set the UV Cutoff', required=True, default=None)
    parser.add_argument('-k', type=float,
                        help='Set the IR Cutoff', default=1e-2)
    parser.add_argument(
        '-T', type=float, help='Temperature of the flow', required=True)
    parser.add_argument(
        '-mu', type=float, help='chemical potential of the flow', required=True)
    parser.add_argument(
        '--max', type=float, help='sigma max', default=6.0
    )
    parser.add_argument(
        '-d', type=int, help='sigma max', default=1
    )

    args = parser.parse_args()
    sigma_max = args.max
    delta_sigma = 0.006
    Lambda = args.L
    kir = args.k
    samples = 1000
    mu = args.mu
    T = args.T
    N_Flavor = args.N if args.N is not None else np.inf
    grid_points = np.arange(0.0, sigma_max, delta_sigma)
    N_Grid = len(grid_points)
    if N_Grid > 1000:
        print(
            f'Number of gridpoints is {N_Grid}, this is larger than 1000 so using unevenly spaced grid points instead')
        grid_points = create_inhomogenious_grid_from_cell_spacing(
            sigma_max, delta_sigma)

    sol = get_model(args.d)(grid_points, Lambda, kir,
                            samples, mu, T, N_Flavor=N_Flavor)

    xmax_show = 2
    x = np.array(sol.return_data.grid)
    y_show = np.array(sol.return_data.solution[-1])[x <= xmax_show]
    x_show = x[x <= xmax_show]

    potentials = [Potential(x, y) for time, y in zip(
        sol.return_data.time, sol.return_data.solution)]
    potential = potentials[-1]

    fig, (u_plot, U_plot, sigma_0_plot) = plt.subplots(1, 3, tight_layout=True)
    u_plot.plot(x_show, y_show)
    u_plot.set_xlabel(r'$\sigma$')
    u_plot.set_ylabel(r'u')
    u_plot.grid()

    U_plot.plot(x_show, potential.U(x_show))
    U_plot.set_xlabel(r'$\sigma$')
    U_plot.set_ylabel(r'U')

    sigma_0_plot.plot(Lambda * np.exp(-np.array(sol.return_data.time)), [
                      elem.sigma for elem in potentials])
    sigma_0_plot.set_xscale('log')
    sigma_0_plot.set_xlabel('t')
    sigma_0_plot.set_ylabel(r'$\sigma_0$')

    plt.show()


if __name__ == "__main__":
    main()
