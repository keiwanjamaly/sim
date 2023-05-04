from ctypes import CDLL
import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple
from math import sqrt, atanh

from c_Grid import Grid_Interface
from c_Return_Data import Return_Data_Interface
from c_Computation_Data import Computation_Data_Interface
from c_Computation import Computation_Interface
from c_Gross_Neveu import Gross_Neveu_Interface, Gross_Neveu_Interface_Vacuum
from compute_observables import DataClass


class GN_1p1():
    d = 1
    d_gamma = 2
    tolerances = 1e-10

    def __init__(self, sigma_max, Lambda, kir, N_Grid, samples, mu, T,
                 N_Flavor=np.infty, h=1, one_over_g2=None, sigma_0=1) -> None:

        if T == 0:
            if mu != 0:
                raise RuntimeError(
                    f'There is no implementation for T=0 and mu != 0, mu is {mu}')
            self.lib = CDLL("./build/src/libgross_neveu_vacuum.dylib")
        else:
            self.lib = CDLL("./build/src/libgross_neveu.dylib")

        print(self.dimension(), self.dimension_gamma())

        grid_points = np.linspace(0, sigma_max, N_Grid)
        grid = Grid_Interface(self.lib, grid_points)
        self.return_data = Return_Data_Interface(
            self.lib, samples, grid.pointer)

        # calculate 1/g**2 if not set by user
        if (one_over_g2 is None):
            self.one_over_g2 = self.calculate_one_g2(h, sigma_0, Lambda)
        else:
            self.one_over_g2 = one_over_g2

        print(self.one_over_g2)

        if T == 0:
            physics_data = Gross_Neveu_Interface_Vacuum(
                self.lib, h, self.one_over_g2, self.dimension(), self.dimension_gamma(), N_Flavor)
        else:
            physics_data = Gross_Neveu_Interface(
                self.lib, h, self.one_over_g2, self.dimension(),
                self.dimension_gamma(), mu, T, N_Flavor)

        computation_data = Computation_Data_Interface(
            self.lib, Lambda, kir, grid.pointer, physics_data.pointer,
            self.tolerances)

        Computation_Interface(
            self.lib, computation_data.pointer, self.return_data.pointer)

    @classmethod
    def calculate_one_g2(cls, h, sigma_0, Lambda) -> float:
        tmp = 1/sqrt(1.0 + (h * sigma_0 / Lambda)**2)
        return cls.d_gamma * (1 / (2 * np.pi)) * (atanh(tmp) - tmp)

    @classmethod
    def dimension(cls):
        return cls.d

    @classmethod
    def dimension_gamma(cls):
        return cls.d_gamma


class GN_2p1(GN_1p1):
    d = 2
    d_gamma = 4
    tolerances = 1e-12

    @classmethod
    def calculate_one_g2(cls, h, sigma_0, Lambda) -> float:
        print('this is getting called')
        tmp = sqrt((h * sigma_0)**2 + Lambda**2)
        tmp = cls.d_gamma * (Lambda**2 + 2 * h * sigma_0 * (h * sigma_0 - tmp)) \
            / (8 * np.pi * tmp)
        return tmp


def main():
    # import json
    #     fig, ((ax1_1, ax1_2, ax1_3), (ax2_1, ax2_2, ax2_3)) = plt.subplots(2, 3)
    #     fig.tight_layout()
    #
    sigma_max = 6.0
    Lambda = 1e5
    kir = 1e-4
    N_Grid = 1000
    samples = 100
    #
    #     def compute_1_1_sol(mu: float, T: float, N_Flavor: int) -> Tuple[float, float]:
    #         sol = GN_1p1(sigma_max, Lambda, kir, N_Grid,
    #                      samples, mu, T, N_Flavor=N_Flavor)
    #         y = sol.return_data.solution
    #         x = sol.return_data.grid
    #         time = sol.return_data.time
    #         return x, y, time
    #
    #     # compare MF results
    #     mu = 0.1
    #     T = 0.1
    #     N_Flavor = np.infty
    #     x, y, time = compute_1_1_sol(mu, T, N_Flavor)
    #     data_class_MF = DataClass(x, time, y)
    #
    #     # test mean field computation for 1+1 dimension
    #     with open('1+1_test_files/flow_MF_T=0.1,mu=0.1.json') as user_file:
    #         parsed_json = json.load(user_file)
    #         x_ref = np.array(parsed_json["graphs"][-1]["x"])
    #         y_ref = np.array(parsed_json["graphs"][-1]["y"])
    #
    #     if (x_ref != x).all():
    #         raise RuntimeWarning("For MF the grids do not match")
    #     ax1_1.plot(x, y[-1], label="computation")
    #     ax1_1.plot(x_ref, y_ref, label="reference", linestyle="dashed")
    #     ax1_1.legend()
    #     ax1_1.grid()
    #     ax1_1.set_xlabel(r'$\sigma$')
    #     difference_y = np.abs(y[-1] - y_ref)
    #     rel_error_y = difference_y / np.abs(y[-1] + y_ref)
    #
    #     ax1_2.plot(x[:-1], difference_y[:-1], label="abs_error")
    #     ax1_2.plot(x[:-1], rel_error_y[:-1], label="rel_error")
    #     ax1_2.legend()
    #     ax1_2.grid()
    #     ax1_2.set_xlabel(r'$\sigma$')
    #
    #     sigma_array = data_class_MF.sigma_0
    #     ax1_3.plot(time, sigma_array)
    #     ax1_3.grid()
    #     ax1_3.set_xlabel(r'$t$')
    #     ax1_3.set_ylabel(r'$\sigma_0$')
    #
    #     # compare N_Flavor = 2 results
    #     mu = 0.1
    #     T = 0.1
    #     N_Flavor = 2
    #     x, y, time = compute_1_1_sol(mu, T, N_Flavor)
    #     data_class_MF = DataClass(x, time, y)
    #
    #     # test N_Flavor = 2 for 1+1 dimension
    #     with open('1+1_test_files/flow_N=2,T=0.1,mu=0.1.json') as user_file:
    #         parsed_json = json.load(user_file)
    #         x_ref = np.array(parsed_json["graphs"][-1]["x"])
    #         y_ref = np.array(parsed_json["graphs"][-1]["y"])
    #
    #     if (x_ref != x).all():
    #         raise RuntimeWarning("For N_flavor=2 the grids do not match")
    #     ax2_1.plot(x, y[-1], label="computation")
    #     ax2_1.plot(x_ref, y_ref, label="reference", linestyle="dashed")
    #     ax2_1.legend()
    #     ax2_1.grid()
    #     ax2_1.set_xlabel(r'$\sigma$')
    #     difference_y = np.abs(y[-1] - y_ref)
    #     rel_error_y = difference_y / np.abs(y[-1] + y_ref)
    #
    #     ax2_2.plot(x[:-1], difference_y[:-1], label="abs_error")
    #     ax2_2.plot(x[:-1], rel_error_y[:-1], label="rel_error")
    #     ax2_2.legend()
    #     ax2_2.grid()
    #     ax2_2.set_xlabel(r'$\sigma$')
    #
    #     sigma_array = data_class_MF.sigma_0
    #     ax2_3.plot(time, sigma_array)
    #     ax2_3.grid()
    #     ax2_3.set_xlabel(r'$t$')
    #     ax2_3.set_ylabel(r'$\sigma_0$')
    #
    #     plt.show()

    # test vacuum solution
    def compute_2_1_sol(mu: float, T: float, N_Flavor: int) -> Tuple[float, float]:
        sol = GN_2p1(sigma_max, Lambda, kir, N_Grid,
                     samples, mu, T, N_Flavor=N_Flavor)
        y = sol.return_data.solution
        x = sol.return_data.grid
        time = sol.return_data.time
        return x, y, time

    mu = 0.0
    T = 0.0
    N_Flavor = np.infty
    x, y, time = compute_2_1_sol(mu, T, N_Flavor)

    # plt.plot(x, y[0], label="initial")
    plt.plot(x, y[-1], label="final")
    plt.xlim([0, 1.2])
    plt.ylim([-.5, .5])
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == "__main__":
    main()
