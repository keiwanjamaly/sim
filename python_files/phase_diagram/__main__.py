# from multiprocessing import Pool, Manager
import numpy as np
import argparse
from python_files.phase_diagram.phase_diagram_computation import phase_diagram_computationa
from python_files.phase_diagram.plot_phase_diagram import density_plot, get_result
from python_files.phase_diagram.phase_boundary import compute_phase_boundary, plot_phase_boundary, plot_tolerance, liftschitz_point
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-N', type=int, help='Set the number of flavors', required=True, default=None)
    parser.add_argument('-L', type=float,
                        help='Set the UV Cutoff', required=True, default=None)
    parser.add_argument('--delta', type=float,
                        help='sets the resolution for the phase diagram', default=0.01)
    parser.add_argument(
        '--plot', help='plot the phase diagram instead of calculating it', action='store_true')
    parser.add_argument(
        '--tol', help='set the tolerance for the phase diagram computation or set the tolerance for the liftischt point computation', type=float, default=None)
    parser.add_argument(
        '--tol_min', help='set the minimum tolerance for the liftschitz point array', type=float, default=None)
    parser.add_argument(
        '--tol_max', help='set the maximum tolerance for the liftschitz point array', type=float, default=None)
    parser.add_argument(
        '--liftschitz', help='compute the liftschitz point', action='store_true')
    args = parser.parse_args()

    if args.N == -1:
        N_Flavor = np.inf
    else:
        N_Flavor = args.N
    Lambda = args.L

    if args.liftschitz:
        result = get_result(Lambda, N_Flavor)
        if args.tol is not None:
            mu, T = liftschitz_point(result, args.tol)
            print(mu, T)
        elif args.tol_min is not None and args.tol_max is not None:
            tol_array = np.arange(args.tol_min, args.tol_max, 0.1)
            plot_tolerance(result, tol_array)
        else:
            raise RuntimeError(
                "Tolerance has to be defined by --tol when computing the liftschitz point")
    elif args.plot:
        result = get_result(Lambda, N_Flavor)
        density_plot(result)
        if args.tol is not None:
            boundary_first_order, boundary_second_order = compute_phase_boundary(
                result, args.tol)
            plot_phase_boundary(boundary_first_order, boundary_second_order)

        # other plot parameters
        plt.title(f'{Lambda=}, {N_Flavor=}')
        plt.xlabel(r'$\mu$')
        plt.ylabel(r'$T$')
        plt.show()
    else:
        phase_diagram_resolution = args.delta
        phase_diagram_computationa(Lambda, N_Flavor, phase_diagram_resolution)


if __name__ == "__main__":
    main()
