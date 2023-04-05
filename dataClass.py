import h5py
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import numpy as np


class FlowData:
    def __init__(self, file_name):
        self.file_name = file_name

    def __enter__(self):
        self.file = h5py.File(self.file_name, 'r')
        self.keys = self.file.attrs.keys()
        self.obj = {}
        for key, item in self.file.attrs.items():
            self.obj[key] = item
        # print(self.obj)
        # self.Lambda = self.file.attrs["Lambda"]
        # self.N_flavor = self.file.attrs["NFlavor"]
        # self.mu = self.file.attrs["mu"]
        # self.T = self.file.attrs["T"]
        # self.NGrid = self.file.attrs["NGrid"]
        # self.sigmaMax = self.file.attrs["sigmaMax"]
        # self.spatial_dimension = self.file.attrs["spatial_dimension"]
        # self.extrapolation_order = "linear" if self.file.attrs["extrapolation_order"] == 1 else "quadratic"
        # self.computation_time = self.file.attrs["computation_time"]
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is not None:
            print(f"An error occurred: {exc_type} - {exc_val}")
        self.file.close()

    def create_parameter_list(self):
        main_list = [f'Λ = {self.obj["Lambda"]:.1e} -> {self.obj["kir"]:.1e}',
                     f'NGrid = {self.obj["NGrid"]}',
                     f'sigma_max = {self.obj["sigmaMax"]:.1f}',
                     f'tolerances = {self.obj["tolerance"]}',
                     f'time = {self.obj["computation_time"]:.2f}']

        parameter_list = []
        # append the other parameters in self.obj
        for key, item in self.obj.items():
            if key not in ["Lambda", "kir", "NGrid", "sigmaMax", "computation_time", "tolerance"]:
                # if item is a float, round it to 2 digits
                if isinstance(item, float):
                    parameter_list.append(f'{key} = {item:.2f}')
                else:
                    parameter_list.append(f'{key} = {item}')

        return main_list, parameter_list

    def __str__(self) -> str:
        # use the create_parameter_list method to create a string
        main_list, parameter_list = self.create_parameter_list()
        # join both lists, that the main_list is on top and the parameter_list below and return the string
        main_list_print = '; '.join(main_list)
        parameter_list_print = '; '.join(parameter_list)
        return f'{main_list_print}\nAdditional parameters: {parameter_list_print}'

    def enter_plot_title(self, fig: Figure):
        # set the suptitle of the figure using the __str__ method
        fig.suptitle(self.__str__())

    def plot_Q(self, ax: plt.Axes, part="Q"):
        x = self.file["grid"]
        y = self.file["t"]
        X, Y = np.meshgrid(x, y)
        Z = self.file[part]

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
            f'Cell spacing across grid\nfor {self.obj["NGrid"]} cells in total.')

    def plot_max_Q_of_x(self, ax: plt.Axes, label="", color=None):
        max_Q = np.empty_like(self.file["grid"])
        Q = self.file["Q"][:, :]
        for i in range(len(self.file["grid"])):
            max_Q[i] = np.max(Q[:, i])

        ax.plot(self.file["grid"][:], max_Q,
                label=label, c=color)
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
        ax.set_ylim(-0.1, 1.0)
        ax.grid()
        ax.set_title(
            r'flow of expectation values')
        ax.legend()

    def add_massSquare_at_pos(self, ax: plt.Axes, pos=-1, color=None):
        ax.scatter([self.sigmaMax], [
                   self.file["massSquare"][pos]], color=color, s=10)
        ax.set_xlabel(r'$\sigma_{max}$')
        ax.set_ylabel(r'$m_{\sigma}^2$')

    def add_massSquare_vs_Lambda_at_pos(self, ax: plt.Axes, pos=-1, color=None):
        ax.scatter([self.Lambda], [
                   self.file["massSquare"][pos]], color=color, s=10)
        ax.set_xlabel(r'$\Lambda$')
        ax.set_ylabel(r'$m_{\sigma}^2$')

    def add_computation_time(self, ax: plt.Axes, color=None):
        ax.scatter([self.Lambda], [self.computation_time], color=color, s=10)
        ax.set_xlabel(r'$\Lambda$')
        ax.set_ylabel('time in [s]')
