from ctypes import CDLL, Structure, c_int, c_double, POINTER, c_void_p

from c_Grid import Grid


class Computation_Data(Structure):
    _fields_ = [("Lambda", c_double),
                ("kir", c_double),
                ("tir", c_double),
                ("tolerances", c_double),
                ("computation_grid", POINTER(Grid)),
                ("data", c_void_p)]


class Computation_Data_Interface():
    def __init__(self, lib: CDLL, Lambda: float, kir: float, grid: POINTER(Grid), physics_data: c_void_p, tolerances: float):
        self.initialize_computation_data = lib.initialize_computation_data
        self.initialize_computation_data.argtypes = [
            c_double, c_double, POINTER(Grid), c_void_p, c_double
        ]
        self.initialize_computation_data.restype = POINTER(Computation_Data)

        self.computation_data_pointer = self.initialize_computation_data(
            Lambda, kir, grid, physics_data, tolerances)

        self.destroy_computation_data = lib.destroy_computation_data
        self.destroy_computation_data.argtypes = [POINTER(Computation_Data)]
        self.destroy_computation_data.restype = None

    def __del__(self):
        self.destroy_computation_data(self.computation_data_pointer)

    @property
    def pointer(self):
        return self.computation_data_pointer
