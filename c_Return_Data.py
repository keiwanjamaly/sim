from ctypes import CDLL, Structure, c_int, c_double, POINTER
from c_Grid import Grid


class Return_Data(Structure):
    _fields_ = [("grid_size", c_int),
                ("samples", c_int),
                ("grid", POINTER(c_double)),
                ("solution_y", POINTER(POINTER(c_double))),
                ("solution_time", POINTER(c_double))]


class Return_Data_Interface():
    def __init__(self, lib: CDLL, samples: int, grid: POINTER(Grid)) -> None:
        self.create_return_data = lib.create_return_data
        self.create_return_data.argtypes = [
            c_int, POINTER(Grid)]
        self.create_return_data.restype = POINTER(Return_Data)

        self.return_data_pointer = self.create_return_data(samples, grid)

        self.destroy_return_data = lib.destroy_return_data
        self.destroy_return_data.argtypes = [POINTER(Return_Data)]
        self.destroy_return_data.restype = None

        self.grid_generated = False
        self.solution_generated = False

    # def get_data():
    #     self.return_data_pointer

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
        if not self.grid_generated:
            self.generated_grid = [
                self.return_data_pointer.contents.grid[i] for i in range(self.grid_size)]
            self.grid_generated = True
        return self.generated_grid

    @property
    def solution(self):
        def generate_solution(i: int):
            return [self.return_data_pointer.contents.solution_y[i][j] for j in range(self.grid_size)]
        if not self.solution_generated:
            self.generated_solution = [
                generate_solution(i) for i in range(self.samples)
            ]
            self.solution_generated
        return self.generated_solution
