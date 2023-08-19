from ctypes import CDLL, Structure, c_int, c_double, POINTER
from python_files.interfaces.c_Grid import Grid


class Return_Data(Structure):
    _fields_ = [("grid_size", c_int),
                ("samples", c_int),
                ("grid", POINTER(c_double)),
                ("solution_y", POINTER(POINTER(c_double))),
                ("solution_time", POINTER(c_double)),
                ("diffusion", POINTER(POINTER(c_double))),
                ("source", POINTER(POINTER(c_double)))]


class Return_Data_Interface():
    def __init__(self, lib: CDLL, samples: int, grid: POINTER(Grid)) -> None:
        self.generated_grid = None
        self.generated_time = None
        self.generated_solution = None
        self.generated_source = None
        self.generated_diffusion = None

        self.create_return_data = lib.create_return_data
        self.create_return_data.argtypes = [
            c_int, POINTER(Grid)]
        self.create_return_data.restype = POINTER(Return_Data)

        self.return_data_pointer = self.create_return_data(samples, grid)

        self.destroy_return_data = lib.destroy_return_data
        self.destroy_return_data.argtypes = [POINTER(Return_Data)]
        self.destroy_return_data.restype = None

    def __del__(self):
        self.destroy_return_data(self.return_data_pointer)

    @property
    def pointer(self) -> POINTER(Return_Data):
        return self.return_data_pointer

    @property
    def grid_size(self) -> int:
        return self.return_data_pointer.contents.grid_size

    @property
    def samples(self) -> int:
        return self.return_data_pointer.contents.samples

    @property
    def grid(self):
        if self.generated_grid is None:
            self.generated_grid = [
                self.return_data_pointer.contents.grid[i] for i in range(self.grid_size)]
        return self.generated_grid

    @property
    def time(self):
        if self.generated_time is None:
            self.generated_time = [
                self.return_data_pointer.contents.solution_time[i] for i in range(self.samples)]
        return self.generated_time

    @property
    def solution(self):
        def generate_solution(i: int):
            return [self.return_data_pointer.contents.solution_y[i][j] for j in range(self.grid_size)]
        if self.generated_solution is None:
            self.generated_solution = [
                generate_solution(i) for i in range(self.samples)
            ]
        return self.generated_solution

    @property
    def source(self):
        def generate_source(i: int):
            return [self.return_data_pointer.contents.source[i][j] for j in range(self.grid_size)]
        if self.generated_source is None:
            self.generated_source = [
                generate_source(i) for i in range(self.samples)
            ]
        return self.generated_source

    @property
    def diffusion(self):
        def generate_diffusion(i: int):
            return [self.return_data_pointer.contents.diffusion[i][j] for j in range(self.grid_size)]
        if self.generated_diffusion is None:
            self.generated_diffusion = [
                generate_diffusion(i) for i in range(self.samples)
            ]
        return self.generated_diffusion
