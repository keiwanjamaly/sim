# from multiprocessing import Pool, Manager
import numpy as np
import argparse
from python_files.phase_diagram.phase_diagram_computation import phase_diagram_computationa
# from python_files.phase_diagram.phase_boundary import compute_phase_boundary, plot_phase_boundary, plot_tolerance, liftschitz_point
from python_files.phase_diagram.plot_phase_diagram import plot_observables
from python_files.phase_diagram.compute_observables import compute_observables
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
        '--compute', help='compute observables of phase diagram', action='store_true')
    args = parser.parse_args()

    if args.N == -1:
        N_Flavor = np.inf
    else:
        N_Flavor = args.N
    Lambda = args.L

    if args.compute:
        compute_observables(Lambda, N_Flavor)

    if args.plot:
        plot_observables(Lambda, N_Flavor)
    else:
        phase_diagram_resolution = args.delta
        phase_diagram_computationa(Lambda, N_Flavor, phase_diagram_resolution)


if __name__ == "__main__":
    main()
