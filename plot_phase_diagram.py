import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import LinearNDInterpolator
import argparse
import h5py


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-N', type=int, help='Set the number of flavors', required=True, default=None)
    parser.add_argument('-L', type=float,
                        help='Set the UV Cutoff', required=True, default=None)
    args = parser.parse_args()
    Lambda = args.L
    if args.N != -1:
        N_Flavor = args.N
    else:
        N_Flavor = np.inf
    filename = f'./data/phase_diagram/phase_diagram_Lambda_{Lambda}_N_Flavor_{N_Flavor}.hdf5'
    with h5py.File(filename, "r") as f:
        result = f["phase_diagram"][:]
        mu_max = np.max(result[:, 0])
        T_max = np.max(result[:, 1])
        # plot the results
        points = list(zip(result[:, 0], result[:, 1]))
        values = list(result[:, 2])
        interpolation = LinearNDInterpolator(points, values)
        X, Y = np.meshgrid(np.linspace(0.0, mu_max, 1000),
                           np.linspace(0.01, T_max, 1000))
        Z = interpolation(X, Y)
        plt.pcolormesh(X, Y, Z, shading='auto')
        plt.colorbar()
        plt.title(f'{Lambda=}, {N_Flavor=}')
        plt.xlabel(r'$\mu$')
        plt.ylabel(r'$T$')

        plt.show()


if __name__ == "__main__":
    main()
