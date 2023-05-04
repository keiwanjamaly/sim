from ctypes import CDLL, Structure, c_int, c_double, POINTER


class Gross_Neveu(Structure):
    _fields_ = [("h", c_double),
                ("one_over_g2", c_double),
                ("dimension", c_int),
                ("dimension_gamma", c_int),
                ("A_d", c_double),
                ("N", c_double),
                ("one_over_N", c_double),
                ("mu", c_double),
                ("beta", c_double)]


class Gross_Neveu_Vacuum(Structure):
    _fields_ = [("h", c_double),
                ("one_over_g2", c_double),
                ("dimension", c_int),
                ("dimension_gamma", c_int),
                ("A_d", c_double),
                ("N", c_double),
                ("one_over_N", c_double)]


class Gross_Neveu_Interface():
    def __init__(self, lib: CDLL, h: float, one_over_g2: float, dimension: int, dimension_gamma: int, mu: float, T: float, N: float) -> None:
        self.initialize_physics_data = lib.initialize_physics_data
        self.initialize_physics_data.argtypes = [
            c_double, c_double, c_int, c_int, c_double, c_double, c_double]
        self.initialize_physics_data.restype = POINTER(Gross_Neveu)

        self.physics_data_pointer = self.initialize_physics_data(
                h, one_over_g2, dimension, dimension_gamma, mu, T, N)

        self.free_physics_data = lib.free_physics_data
        self.free_physics_data.argtypes = [POINTER(Gross_Neveu)]
        self.free_physics_data.restype = None

    def __del__(self):
        self.free_physics_data(self.physics_data_pointer)

    @property
    def pointer(self):
        return self.physics_data_pointer


class Gross_Neveu_Interface_Vacuum():
    def __init__(self, lib: CDLL, h: float, one_over_g2: float, dimension: int, dimension_gamma: int, N: float) -> None:
        self.initialize_physics_data = lib.initialize_physics_data
        self.initialize_physics_data.argtypes = [
            c_double, c_double, c_int, c_int, c_double]
        self.initialize_physics_data.restype = POINTER(Gross_Neveu_Vacuum)

        self.physics_data_pointer = self.initialize_physics_data(
                h, one_over_g2, dimension, dimension_gamma, N)

        self.free_physics_data = lib.free_physics_data
        self.free_physics_data.argtypes = [POINTER(Gross_Neveu_Vacuum)]
        self.free_physics_data.restype = None

    def __del__(self):
        self.free_physics_data(self.physics_data_pointer)

    @property
    def pointer(self):
        return self.physics_data_pointer
