from python_files.gross_neveu.Gross_Neveu import get_model
from python_files.grid_creator import create_homogenious_grid, create_inhomogenious_grid_from_cell_spacing
from python_files.gross_neveu.couplings.couplings_io import get_closest_coupling_approximation_or_MF_coupling
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
        '--max', type=float, help='sigma max', default=2000.0
    )
    parser.add_argument(
        '-d', type=int, help='sigma max', default=2
    )
    parser.add_argument(
        '-g2_mf', help='forces_inital_mean_field_coupling', action='store_true')
    parser.add_argument(
        '-u_out', help='write the results to stdout', action='store_true')
    parser.add_argument(
        '-U_out', help='write the results to stdout', action='store_true')
    parser.add_argument(
        '-diff_out', help='write the results to stdout', action='store_true')
    parser.add_argument(
        '-source_out', help='write the results to stdout', action='store_true')

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
        # print(
        #     f'Number of gridpoints is {N_Grid}, this is larger than 1000 so using unevenly spaced grid points instead')
        grid_points = create_inhomogenious_grid_from_cell_spacing(
            sigma_max, delta_sigma)

    if args.g2_mf:
        one_over_g2 = get_model(args.d).calculate_one_g2(
            h=1, sigma_0=1.0, Lambda=Lambda)
    else:
        one_over_g2 = get_closest_coupling_approximation_or_MF_coupling(
            Lambda, N_Flavor, './data', args.d)

    sol = get_model(args.d)(grid_points, Lambda, kir,
                            samples, mu, T, N_Flavor=N_Flavor, one_over_g2=one_over_g2)

    x = np.array(sol.return_data.grid)
    potentials = [Potential(x, y) for time, y in zip(
        sol.return_data.time, sol.return_data.solution)]
    potential = potentials[-1]
    U = potential.U(x)
    max_diff = np.max(sol.return_data.diffusion, axis=0)
    max_source = np.max(sol.return_data.source, axis=0)

    xmax_show = 1.2
    y_show = np.array(sol.return_data.solution[-1])[x <= xmax_show]
    x_show = x[x <= xmax_show]
    U_show = np.array(U)[x <= xmax_show]

    if not (args.u_out or args.U_out or args.diff_out or args.source_out):
        fig, ((u_plot, U_plot, sigma_0_plot), (plot_source, plot_diffusion, _)) = plt.subplots(
            2, 3, tight_layout=True)
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

        element = 1
        plot_source.plot(x, sol.return_data.source[element])
        plot_diffusion.plot(x, sol.return_data.diffusion[300])
        plot_diffusion.set_xscale('log')
        plot_diffusion.set_yscale('symlog')
        plot_diffusion.set_title(sol.return_data.time[300])
        plot_diffusion.grid()

        plt.show()
    else:
        for sigma, u, U, diff, source in zip(x, sol.return_data.solution[-1], U, max_diff, max_source):
            # if sigma <= 2.0:
            # print(sigma, u, U, sep='\t')

            print(sigma, end='')
            if args.u_out:
                print(end='\t')
                print(u, end='')
            if args.U_out:
                print(end='\t')
                print(U, end='')
            if args.diff_out:
                print(end='\t')
                print(diff, end='')
            if args.source_out:
                print(end='\t')
                print(np.abs(source), end='')
            print(end='\n')


if __name__ == "__main__":
    main()
