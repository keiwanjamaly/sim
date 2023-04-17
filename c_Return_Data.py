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
        self.initialize_return_data = lib.initialize_return_data
        self.initialize_return_data.argtypes = [
            c_int, POINTER(Grid)]
        self.initialize_return_data.restype = POINTER(Return_Data)

        self.return_data_pointer = self.initialize_return_data(samples, grid)

        self.free_return_data = lib.free_return_data
        self.free_return_data.argtypes = [POINTER(Return_Data)]
        self.free_return_data.restype = None

    # def get_data():
    #     self.return_data_pointer

    def __del__(self):
        print("freeing return_data")
        self.free_return_data(self.return_data_pointer)

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
        return self.return_data_pointer.contents.grid
