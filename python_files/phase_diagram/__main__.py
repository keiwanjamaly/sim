# from multiprocessing import Pool, Manager
import numpy as np
import argparse
from python_files.phase_diagram.plot_phase_diagram import plot_observables
from python_files.phase_diagram.compute_observables import compute_observables
from python_files.phase_diagram.plots_for_poster import plots_poster
from python_files.phase_diagram.plot_single_potential import plot_potential
from python_files.phase_diagram.phase_diagram_computation import phase_diagram_computation
from python_files.phase_diagram.save_phase_diagram import save_final_phase_diagram_to_file
from python_files.phase_diagram.compute_phase_boundary import compute_boundary
from python_files.phase_diagram.compute_phase_boundary import compute_boundary_mu0
from python_files.phase_diagram.compute_phase_boundary import compute_boundary_reference
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-N', type=int, help='Set the number of flavors', required=True, default=None)
    parser.add_argument('-L', type=float,
                        help='Set the UV Cutoff', required=True, default=None)
    parser.add_argument('--delta', type=float,
                        help='sets the resolution for the phase diagram', default=0.01)
    parser.add_argument('-kir', type=float,
                        help='set the IR cutoff for the phase boundary calculation', default=1e-2)
    parser.add_argument(
        '--plot', help='plot the phase diagram instead of calculating it', action='store_true')
    parser.add_argument(
        '--save', help='set the output file type', type=argparse.FileType('w', encoding='latin-1'))
    parser.add_argument(
        '--save_LP', help='file of the LP to save', type=argparse.FileType('w', encoding='latin-1'))
    parser.add_argument(
        '--compute', help='compute observables of phase diagram', action='store_true')
    parser.add_argument(
        '--poster', help='create plots for the poster', action='store_true')
    parser.add_argument(
        '--potential', help='plot a single potential', nargs=2, type=float)
    parser.add_argument(
        '--mf_contour', help='calculate the contour of the mean field phase diagram', action='store_true')
    parser.add_argument(
        '--boundary', help='calculate the phase boundary', action='store_true')
    parser.add_argument(
        '--boundary_mu0', help='compute the phase boundary for mu = 0 and multiple flavors', action='store_true')
    parser.add_argument(
        '-f', help='filename of the h5py file', type=str)
    args = parser.parse_args()

    if args.N == -1:
        N_Flavor = np.inf
    else:
        N_Flavor = args.N
    if args.L == -1:
        Lambda = np.inf
    else:
        Lambda = args.L

    skip_computation = False

    if args.f is not None and args.save is not None:
        save_final_phase_diagram_to_file(args.f, args.save)
        skip_computation = True
    # if args.mf_contour:
    #     compute_contour()
    #     skip_computation = True

    if args.potential:
        mu, T = args.potential
        plot_potential(mu, T, Lambda, N_Flavor)
        skip_computation = True

    if args.poster:
        plots_poster()
        skip_computation = True

    if args.compute:
        compute_observables(Lambda, N_Flavor)
        skip_computation = True

    if args.plot:
        plot_observables(Lambda, N_Flavor, args.save)
        skip_computation = True

    if args.boundary:
        compute_boundary(Lambda, N_Flavor, args.save, args.save_LP, args.kir)
        skip_computation = True

    if args.boundary_mu0:
        compute_boundary_mu0(Lambda, args.save)
        skip_computation = True

    if not skip_computation:
        phase_diagram_resolution = args.delta
        phase_diagram_computation(Lambda, N_Flavor, phase_diagram_resolution)


if __name__ == "__main__":
    main()
