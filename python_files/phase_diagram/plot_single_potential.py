from python_files.phase_diagram.file_io import get_filename
import h5py
import numpy as np
import matplotlib.pyplot as plt
from python_files.data_class import Potential


def plot_potential(mu, T, Lambda, N_Flavor):
    filename = get_filename(Lambda, N_Flavor)
    with h5py.File(filename, "r") as file:
        mu_array = file["phase_diagram"][:, 0]
        T_array = file["phase_diagram"][:, 1]
        index = np.where(np.logical_and(T == T_array, mu == mu_array))[0]
        # mu_index = set(list(np.where(mu == mu_array)))
        grid = file["grid"][:]
        u = file["u"][index][0]

        x_plot = 1.2
        grid_plot = grid[grid <= x_plot]
        u_plot = u[grid <= x_plot]

        potential = Potential(grid, u)
        min = potential.sigma

        plt.grid()
        plt.plot([min], [0], 'rx')
        plt.title(f'$\\mu = {mu}, T = {T}$')
        plt.plot(grid_plot, u_plot)
        plt.show()
        # print(mu_index, T_index)
