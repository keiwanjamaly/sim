from python_files.interfaces.c_Grid import Grid_Interface
from python_files.interfaces.c_Gross_Neveu import Gross_Neveu_Interface, Gross_Neveu_Interface_Vacuum
from python_files.interfaces.c_Computation import Computation_Interface
from python_files.interfaces.c_Computation_Data import Computation_Data_Interface
from python_files.interfaces.c_Return_Data import Return_Data_Interface
from ctypes import CDLL
import numpy as np
from math import sqrt, atanh
import platform


class GN_1p1():
    d = 1
    d_gamma = 2
    tolerances = 1e-10

    def __init__(self, grid_points, Lambda, kir, samples, mu, T,
                 N_Flavor=np.infty, h=1, one_over_g2=None, sigma_0=1) -> None:

        if T == 0:
            if mu != 0:
                raise RuntimeError(
                    f'There is no implementation for T=0 and mu != 0, mu is {mu}')
            if platform.system() == 'Linux':
                self.lib = CDLL("./build/src/libgross_neveu_vacuum.so")
            elif platform.system() == 'Darwin':
                self.lib = CDLL("./build/src/libgross_neveu_vacuum.dylib")
            else:
                raise RuntimeError('Platform not defined')
        else:
            if platform.system() == 'Linux':
                self.lib = CDLL("./build/src/libgross_neveu.so")
            elif platform.system() == 'Darwin':
                self.lib = CDLL("./build/src/libgross_neveu.dylib")
            else:
                raise RuntimeError('Platform not defined')

        grid = Grid_Interface(self.lib, grid_points)
        self.return_data = Return_Data_Interface(
            self.lib, samples, grid.pointer)

        # calculate 1/g**2 if not set by user
        if (one_over_g2 is None):
            self.one_over_g2 = self.calculate_one_g2(h, sigma_0, Lambda)
        else:
            self.one_over_g2 = one_over_g2

        if T == 0:
            physics_data = Gross_Neveu_Interface_Vacuum(
                self.lib, h, self.one_over_g2, self.dimension(),
                self.dimension_gamma(), N_Flavor)
        else:
            physics_data = Gross_Neveu_Interface(
                self.lib, h, self.one_over_g2, self.dimension(),
                self.dimension_gamma(), mu, T, N_Flavor)

        computation_data = Computation_Data_Interface(
            self.lib, Lambda, kir, grid.pointer, physics_data.pointer,
            self.tolerances)

        self.computation = Computation_Interface(
            self.lib, computation_data.pointer, self.return_data.pointer)

        self.computation_status_code = self.computation.computation_status

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
        tmp = sqrt((h * sigma_0)**2 + Lambda**2)
        tmp = cls.d_gamma * (Lambda**2 + 2 * h * sigma_0 * (h * sigma_0 - tmp)) \
            / (8 * np.pi * tmp)
        return tmp


def get_model(dimension: int):
    if dimension == 1:
        model = GN_1p1
    elif dimension == 2:
        model = GN_2p1
    else:
        raise RuntimeError(
            "Gross Neveu Model is not implemented for more than 2 spatial dimensions")
    return model
