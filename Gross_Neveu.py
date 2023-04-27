from ctypes import CDLL
import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple

from c_Grid import Grid_Interface
from c_Return_Data import Return_Data_Interface
from c_Computation_Data import Computation_Data_Interface
from c_Computation import Computation_Interface
from c_Gross_Neveu import Gross_Neveu_Interface


class GN_1p1():
    d = 1
    d_gamma = 2
    tolerances = 1e-10

    def __init__(self, sigma_max, Lambda, kir, N_Grid, samples, mu, T,
                 N_Flavor=np.infty, h=1, sigma_0=1) -> None:
        self.lib = CDLL("./build/src/libgross_neveu.dylib")

        self.dk = self.dimension()

        grid_points = np.linspace(0, sigma_max, N_Grid)
        grid = Grid_Interface(self.lib, grid_points)
        self.return_data = Return_Data_Interface(self.lib, samples, grid.pointer)

        physics_data = Gross_Neveu_Interface(
            self.lib, Lambda, h, sigma_0, self.dimension(),
            self.dimension_gamma(), mu, T, N_Flavor)

        computation_data = Computation_Data_Interface(
            self.lib, Lambda, kir, grid.pointer, physics_data.pointer,
            self.tolerances)

        Computation_Interface(
            self.lib, computation_data.pointer, self.return_data.pointer)

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


def main():
    fig, ((ax1, ax2), (ax3)) = plt.subplots(2, 2)

    sigma_max = 6.0
    Lambda = 1e5
    kir = 1e-4
    N_Grid = 1000
    samples = 100
    N_Flavor = 2

    def compute_sol(mu: float, T: float) -> Tuple[float, float]:
        sol = GN_1p1(sigma_max, Lambda, kir, N_Grid, samples, mu, T, N_Flavor=N_Flavor)
        y = sol.return_data.solution[-1]
        x = sol.return_data.grid
        return x, y

    mu = 0.1
    T = 0.1
    x, y = compute_sol(mu, T)

    ax1.plot(x, y)
    ax1.grid()

    plt.show()


    # print("hello")
    # print(sol)


if __name__ == "__main__":
    main()
