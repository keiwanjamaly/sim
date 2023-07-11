import h5py
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator
import numpy as np


def get_result(Lambda, N_Flavor):
    filename = f'./data/phase_diagram/phase_diagram_Lambda_{Lambda}_N_Flavor_{N_Flavor}.hdf5'
    with h5py.File(filename, "r") as f:
        result = f["phase_diagram"][:]
    return result


def density_plot(result):
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
