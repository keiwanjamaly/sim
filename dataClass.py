import h5py
import matplotlib.pyplot as plt
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

    def plot_Q(self):
        x = self.file["grid"]
        y = self.file["t"]
        X, Y = np.meshgrid(x, y)
        Z = self.file["Q"]

        # fig, ax = plt.subplots()
        img = plt.pcolor(X, Y, Z)  # , vmin=-0.5, vmax=1.0)
        plt.xscale("log")
        plt.colorbar(img)
        plt.show()

    def get_max_Q_of_x(self):
        max_Q = np.empty_like(self.file["grid"])
        Q = self.file["Q"][:, :]
        for i in range(len(self.file["grid"])):
            max_Q[i] = np.max(Q[:, i])

        return self.file["grid"][:], max_Q

    def get_sigma_max(self):
        return self.file.attrs["sigmaMax"]
