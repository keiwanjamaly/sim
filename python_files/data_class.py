import numpy as np
from scipy.integrate import cumulative_trapezoid
from typing import Tuple
from scipy.interpolate import interp1d


class Line:
    def __init__(self, m: float, c: float) -> None:
        self.m = m
        self.c = c

    def __call__(self, x: float) -> float:
        return self.m*x + self.c

    def get_root(self) -> float:
        return - self.c / self.m

    def get_area(self, x_start: float, x_end: float) -> float:
        return 0.5 * self.m * (x_end**2 - x_start**2) + self.c * (x_end - x_start)

    @classmethod
    def from_points(cls, f_1: float, x_1: float, f_2: float, x_2: float) -> "Line":
        m = (f_1 - f_2) / (x_1 - x_2)
        c = f_1 - m*x_1
        return cls(m, c)

    @staticmethod
    def get_intersection(line_1: "Line", line_2: "Line") -> Tuple[float, float]:
        x = - (line_1.c - line_2.c) / (line_1.m - line_2.m)
        f = line_1(x)
        return f, x


class Potential:
    def __init__(self, grid: np.ndarray, u: np.ndarray):
        self.grid = grid
        self.h = grid[1] - grid[0]
        self.u = u
        self.interval = self.__compute_sign_change_interval()
        self.U = self.__compute_potential()
        self.sigma = self.__sigma()
        self.mass_square = self.__compute_curvature_mass_square()
        self.second_div_at_zero = self.__compute_second_div_at_zero()
        self.forth_div_at_zero = self.__compute_forth_div_at_zero()

    def __compute_second_div_at_zero(self):
        result = self.u[1] / self.h
        return result

    def __compute_forth_div_at_zero(self):
        result = self.u[2] - 2 * self.u[1]
        result /= self.h**3
        return result

    def __compute_curvature_mass_square(self) -> float:
        lower_index, upper_index = self.interval

        if (lower_index, upper_index) != (0, 1):
            lower_index, upper_index = lower_index + 1, upper_index + 1
        else:
            return 0

        lower_point, upper_point = self.grid[lower_index], self.grid[upper_index]
        lower_u, upper_u = self.u[lower_index], self.u[upper_index]

        line = Line.from_points(upper_u, upper_point, lower_u, lower_point)

        return line.m

    def __calculate_position_of_root(self) -> float:
        lower_index, upper_index = self.interval

        if (lower_index, upper_index) == (0, 1):
            return 0

        # Use the next interval for interpolation to account for potential non-analytic behavior of u
        lower_index, upper_index = lower_index + 1, upper_index + 1
        lower_point, upper_point = self.grid[lower_index], self.grid[upper_index]
        lower_u, upper_u = self.u[lower_index], self.u[upper_index]

        line = Line.from_points(upper_u, upper_point, lower_u, lower_point)

        return line.get_root()

    def __compute_potential(self):
        lower_index, upper_index = self.interval

        if (lower_index, upper_index) == (0, 1):
            y_int = cumulative_trapezoid(self.u, self.grid, initial=0)
            points = self.grid
        else:
            # function, which returns the slope of two points
            lower_point, upper_point = self.grid[lower_index - 1], \
                self.grid[upper_index - 1]
            lower_u, upper_u = self.u[lower_index - 1], \
                self.u[upper_index - 1]

            line_l = Line.from_points(
                upper_u, upper_point, lower_u, lower_point)

            lower_point, upper_point = self.grid[lower_index + 1], \
                self.grid[upper_index + 1]
            lower_u, upper_u = self.u[lower_index + 1], \
                self.u[upper_index + 1]

            line_r = Line.from_points(
                upper_u, upper_point, lower_u, lower_point)

            (f_x_0, x_0) = Line.get_intersection(line_l, line_r)

            points = np.zeros(len(self.grid) + 1)
            values = np.zeros(len(self.u) + 1)

            points[:lower_index + 1] = self.grid[:lower_index + 1]
            values[:lower_index + 1] = self.u[:lower_index + 1]
            points[lower_index + 1] = x_0
            values[lower_index + 1] = f_x_0
            points[upper_index + 1:] = self.grid[upper_index:]
            values[upper_index + 1:] = self.u[upper_index:]

            y_int = cumulative_trapezoid(values, points, initial=0)

        interpolation = interp1d(points, y_int, kind='linear')
        return interpolation

    def __compute_sign_change_interval(self):
        interval = (0, 1)
        for i in range(1, len(self.u)):
            if self.u[i - 1] < 0 and self.u[i] > 0:
                interval = (i-1, i)
                break

        return interval

    def __sigma(self):
        potential_root = self.__calculate_position_of_root()
        potential = self.U(potential_root)
        if potential < 0:
            return potential_root

        return 0
