from python_files.gross_neveu.Gross_Neveu import GN_2p1, GN_1p1
from python_files.grid_creator import create_homogenious_grid, create_inhomogenious_grid_from_number_of_cells
import argparse
import numpy as np
import matplotlib.pyplot as plt


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
        '-G', type=int, help='Numer of grid points', default=1000
    )
    parser.add_argument(
        '-max', type=float, help='sigma max', default=6.0
    )
    parser.add_argument(
        '-d', type=int, help='sigma max', default=1
    )

    args = parser.parse_args()
    sigma_max = args.max
    N_Grid = args.G
    Lambda = args.L
    kir = args.k
    samples = 10
    mu = args.mu
    T = args.T
    N_Flavor = args.N if args.N is not None else np.inf
    if N_Grid > 1000:
        print(
            f'Number of gridpoints is {N_Grid}, this is larger than 1000 so using unevenly spaced grid points instead')
        grid_points = create_inhomogenious_grid_from_number_of_cells(
            sigma_max, N_Grid)
    else:
        grid_points = create_homogenious_grid(sigma_max, N_Grid)

    if args.d == 1:
        sol = GN_1p1(grid_points, Lambda, kir,
                     samples, mu, T, N_Flavor=N_Flavor)
    elif args.d == 2:
        sol = GN_2p1(grid_points, Lambda, kir,
                     samples, mu, T, N_Flavor=N_Flavor)
    else:
        raise RuntimeError(
            "there is no three dimensional Gross-Neveu model defined")

    xmax_show = 1.2
    x = np.array(sol.return_data.grid)
    y = np.array(sol.return_data.solution[-1])[x <= xmax_show]
    x = x[x <= xmax_show]

    plt.plot(x, y)
    plt.xlabel(r'$\sigma$')
    plt.ylabel(r'u')
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()
