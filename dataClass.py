import h5py
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import numpy as np


class FlowData:
    def __init__(self, file_name):
        self.file_name = file_name

    def __enter__(self):
        self.file = h5py.File(self.file_name, 'r')
        self.Lambda = self.file.attrs["Lambda"]
        self.N_flavor = self.file.attrs["NFlavor"]
        self.mu = self.file.attrs["mu"]
        self.T = self.file.attrs["T"]
        self.N_grid = self.file.attrs["NGrid"]
        self.sigma_max = self.file.attrs["sigmaMax"]
        self.spatial_dimension = self.file.attrs["spatial_dimension"]
        self.extrapolation = "linear" if self.file.attrs["extrapolation_order"] == 1 else "quadratic"
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        del self.file

    def __str__(self) -> str:
        print_list = [f'Λ = {self.Lambda:.1e} -> {self.file.attrs["kir"]:.1e}',
                      f'N_f = {self.N_flavor}',
                      f'NGrid = {self.file.attrs["NGrid"]}; sigma_max = {self.file.attrs["sigmaMax"]}',
                      f'mu, T = {self.mu}, {self.T}',
                      f'tolerances = {self.N_grid}; time = {self.file.attrs["computation_time"]:.2f}']
        return '\n'.join(print_list)

    def enter_plot_title(self, fig: Figure):
        fig.suptitle(r'$N_f='+f'{self.N_flavor}'+r'\quad$' +
                     r'$\mu=' + f'{self.mu}' + r'\quad$' +
                     r'$T=' + f'{self.T}' + r'$;' + f'\nspatial dimension = {self.spatial_dimension}, extrapolation = {self.extrapolation}')

    def plot_Q(self, ax: plt.Axes, part="Q"):
        x = self.file["grid"]
        y = self.file["t"]
        X, Y = np.meshgrid(x, y)
        Z = self.file[part]

        # ax.set_title("sdf")

        img = ax.pcolormesh(X, Y, Z)  # , vmin=-0.5, vmax=1.0)
        ax.set_ylabel("t")
        ax.set_xlabel(r'$\sigma$')
        clb = plt.colorbar(img, location='top')
        if part == "Q":
            clb.set_label(r'$\partial_{\sigma}\ Q$')
        elif part == "S":
            clb.set_label(r'$S$')

    def plot_grid(self, ax: plt.Axes):
        x = self.file["grid"][:]
        ax.plot((x[:-1] + x[1:])/2, np.ediff1d(x))
        ax.set_xlabel(r'$\sigma$')
        ax.set_title(
            f'Cell spacing across grid\nfor {self.N_grid} cells in total.')

    def plot_max_Q_of_x(self, ax: plt.Axes, label="", color=None):
        max_Q = np.empty_like(self.file["grid"])
        Q = self.file["Q"][:, :]
        for i in range(len(self.file["grid"])):
            max_Q[i] = np.max(Q[:, i])

        ax.plot(self.file["grid"][:], max_Q,
                label=r'$\sigma_{max}$=' + label, c=color)
        ax.set_xlabel(r'$\sigma$')
        ax.set_title(r'max$_t(\partial_{\sigma}\ Q(t, \sigma))$')

    def plot_sigma_and_msquare_flow(self, ax: plt.Axes):
        k = self.file["k"][:]
        sigma = self.file["sigma"][:]
        massSquare = self.file["massSquare"][:]

        ax.plot(k, sigma, label=r'$\sigma$')
        ax.plot(k, massSquare, label=r'$m_{\sigma}^2$')
        ax.set_xlabel(r'$k$')
        ax.set_xscale("log")
        ax.set_ylim(0.0, 1.0)
        ax.set_title(
            r'flow of expectation values')
        ax.legend()

    def add_massSquare_at_pos(self, ax: plt.Axes, pos=-1, color=None):
        ax.scatter([self.sigma_max], [
                   self.file["massSquare"][pos]], color=color, s=10)
        ax.set_xlabel(r'$\sigma_{max}$')
        ax.set_ylabel(r'$m_{\sigma}^2$')
