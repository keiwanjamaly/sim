from ctypes import CDLL, Structure, c_int, c_double, POINTER
from numpy.typing import ArrayLike


class Grid(Structure):
    _fields_ = [("N", c_int),
                ("grid_points", POINTER(c_double)),
                ("dx", POINTER(c_double)),
                ("dx_midpoints", POINTER(c_double)),
                ("grid_left", c_double),
                ("grid_right", c_double)]


class Grid_Interface():
    def __init__(self, lib: CDLL, grid_points: ArrayLike = None):
        self.create_grid = lib.create_grid
        self.create_grid.argtypes = [c_int,
                                     POINTER(c_double)]
        self.create_grid.restype = POINTER(Grid)

        self.destroy_grid = lib.destroy_grid
        self.destroy_grid.argtypes = [POINTER(Grid)]
        self.destroy_grid.restype = None

        if grid_points is not None:
            self.set_grid(grid_points)

    def set_grid(self, grid: ArrayLike):
        grid_array = (c_double * len(grid))(*grid)
        self.grid_pointer = self.create_grid(len(grid_array), grid_array)

    def __del__(self):
        self.destroy_grid(self.grid_pointer)

    @property
    def pointer(self):
        return self.grid_pointer
