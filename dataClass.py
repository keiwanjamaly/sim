import h5py
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np


class FlowData:
    def __init__(self, file_name):
        self.file_name = file_name

    def __enter__(self):
        self.file = h5py.File(self.file_name, 'r')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        del self.file

    def __str__(self) -> str:
        print_list = [f'Λ = {self.file.attrs["Lambda"]:.1e} -> {self.file.attrs["kir"]:.1e}',
                      f'N_f = {self.file.attrs["NFlavor"]}',
                      f'NGrid = {self.file.attrs["NGrid"]}; sigma_max = {self.file.attrs["sigmaMax"]}',
                      f'mu, T = {self.file.attrs["mu"]}, {self.file.attrs["T"]}',
                      f'tolerances = {self.file.attrs["tolerance"]}; time = {self.file.attrs["computation_time"]:.2f}']
        return '\n'.join(print_list)

    def plot_Data(self, data, xlim=None):
        x = self.file["grid"][:]

        if not (xlim is None):
            filter_function = np.logical_and(x >= xlim[0], x <= xlim[-1])
        else:
            filter_function = np.full_like(x, True, dtype=bool)

        # fig, ax = plt.subplots()
        # print(x[filter_function])
        y = self.file["t"]
        X, Y = np.meshgrid(x[filter_function], y)
        Z = self.file[data][:, filter_function]
        # , vmin=-0.5, vmax=1.0)
        img = plt.pcolor(X, Y, Z)
        plt.colorbar(img)
        plt.xlabel(r'$\sigma$')
        plt.ylabel(r'$t$')
        plt.title(
            r'$\partial_{\sigma}Q$ for $\Lambda='+f'{self.file.attrs["Lambda"]:.1e}'+r'$'+f' d={self.file_name[-6]}')
        plt.show()

    def get_max_Q_of_x(self):
        max_Q = np.empty_like(self.file["grid"])
        Q = self.file["Q"][:, :]
        for i in range(len(self.file["grid"])):
            max_Q[i] = np.max(Q[:, i])

        return self.file["grid"][:], max_Q

    def get_sigma_max(self):
        return self.file.attrs["sigmaMax"]
