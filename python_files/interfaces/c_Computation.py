from ctypes import CDLL, Structure, POINTER, c_int
import numpy as np
import platform

from python_files.interfaces.c_Grid import Grid_Interface
from python_files.interfaces.c_Return_Data import Return_Data, Return_Data_Interface
from python_files.interfaces.c_Computation_Data import Computation_Data, Computation_Data_Interface
from python_files.interfaces.c_Gross_Neveu import Gross_Neveu_Interface


class Computation_Interface(Structure):
    def __init__(self, lib: CDLL, computation_data: POINTER(Computation_Data), return_data: POINTER(Return_Data)):
        self.compute = lib.compute
        self.compute.argtypes = [
            POINTER(Computation_Data), POINTER(Return_Data)]
        self.compute.restype = c_int

        result = self.compute(computation_data, return_data)

        self.computation_status = result


def main():

    if platform.system() == 'Linux':
        lib = CDLL("./build/src/libsim.so")
    elif platform.system() == 'Darwin':
        lib = CDLL("./build/src/libsim.dylib")
    else:
        raise RuntimeError('Platform not defined')

    sigma_max = 10
    N_Grid = 11
    samples = 10
    Lambda = 1000
    kir = 1e-4
    tolerances = 1e-10
    h = 1
    sigma_0 = 1
    d = 1
    d_gamma = 2
    mu = 0.1
    T = 0.1
    N_Flavor = np.infty
    # N_Flavor = 2

    grid_points = np.linspace(0, sigma_max, N_Grid)

    grid = Grid_Interface(lib, grid_points)

    return_data = Return_Data_Interface(lib, samples, grid.pointer)

    physics_data = Gross_Neveu_Interface(
        lib, Lambda, h, sigma_0, d, d_gamma, mu, T, N_Flavor)

    computation_data = Computation_Data_Interface(
        lib, Lambda, kir, grid.pointer, physics_data.pointer, tolerances)

    Computation_Interface(lib, computation_data.pointer, return_data.pointer)

    N = return_data.grid_size
    samples = return_data.samples
    print(return_data.samples)
    for i in range(N):
        print(return_data.pointer.contents.grid[i])


if __name__ == "__main__":
    main()
