import numpy as np
from scipy.integrate import trapezoid
from typing import Tuple


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


def get_sign_change_interval_u(u: np.ndarray) -> Tuple[int, int]:
    interval = (0, 1)
    for i in range(1, len(u)):
        if u[i - 1] < 0 and u[i] > 0:
            interval = (i-1, i)
            break

    return interval


def calculate_position_of_root(grid_points: np.ndarray, u: np.ndarray, interval: Tuple[int, int]) -> float:
    lower_index, upper_index = interval

    if (lower_index, upper_index) == (0, 1):
        return 0

    # Use the next interval for interpolation to account for potential non-analytic behavior of u
    lower_index, upper_index = lower_index + 1, upper_index + 1
    lower_point, upper_point = grid_points[lower_index], grid_points[upper_index]
    lower_u, upper_u = u[lower_index], u[upper_index]

    line = Line.from_points(upper_u, upper_point, lower_u, lower_point)

    return line.get_root()


def get_curvature_mass_square(grid_points: np.ndarray, u: np.ndarray, interval: Tuple[int, int]) -> float:
    lower_index, upper_index = interval

    if (lower_index, upper_index) != (0, 1):
        lower_index, upper_index = lower_index + 1, upper_index + 1

    lower_point, upper_point = grid_points[lower_index], grid_points[upper_index]
    lower_u, upper_u = u[lower_index], u[upper_index]

    line = Line.from_points(upper_u, upper_point, lower_u, lower_point)

    return line.m


def get_potential_at_interval_from_u(grid_points: np.ndarray, u: np.ndarray, interval: Tuple[int, int]) -> float:
    root = calculate_position_of_root(grid_points, u, interval)

    lower_index, upper_index = interval

    if (lower_index, upper_index) == (0, 1):
        return 0

    # function, which returns the slope of two points
    lower_point, upper_point = grid_points[lower_index - 1], \
        grid_points[upper_index - 1]
    lower_u, upper_u = u[lower_index - 1], \
        u[upper_index - 1]

    line_l = Line.from_points(upper_u, upper_point, lower_u, lower_point)

    lower_point, upper_point = grid_points[lower_index + 1], \
        grid_points[upper_index + 1]
    lower_u, upper_u = u[lower_index + 1], \
        u[upper_index + 1]

    line_r = Line.from_points(upper_u, upper_point, lower_u, lower_point)

    (_, x_0_prime) = Line.get_intersection(line_l, line_r)

    last_triangle = line_r.get_area(x_0_prime, root)
    first_triangle = line_l.get_area(grid_points[lower_index], x_0_prime)

    area_till_bracketed = trapezoid(u[:lower_index], grid_points[:lower_index])
    return area_till_bracketed + first_triangle + last_triangle

    # result = trapezoid(u[:lower_index], grid_points[:lower_index])
    # print(result)
    # return result


def get_sigma_0(grid_points: np.ndarray, u: np.ndarray, interval: Tuple[float, float]) -> float:
    potential = get_potential_at_interval_from_u(grid_points, u, interval)
    if potential < 0:
        return calculate_position_of_root(grid_points, u, interval)

    return calculate_position_of_root(grid_points, u, (0, 1))


class DataClass:
    def __init__(self, grid: np.ndarray, time: np.ndarray, u: np.ndarray[2]) -> None:
        self.grid = grid
        self.time = time
        self.u = u

        self.sigma_0_generated = False
        self.interval_generated = False
        self.mass_square_generated = False

    @property
    def sigma_0(self) -> np.ndarray:
        if not self.sigma_0_generated:
            self.sigma_0_computed = [get_sigma_0(
                self.grid, u_points, interval) for u_points, interval in zip(self.u, self.interval)]
            self.sigma_0_generated = True

        return self.sigma_0_computed

    @property
    def mass_square(self) -> np.ndarray:
        if not self.mass_square_generated:
            self.mass_square_computed = [get_curvature_mass_square(
                self.grid, u_points, interval) for u_points, interval in zip(self.u, self.interval)]
            self.mass_square_generated = True

        return self.mass_square_computed

    @property
    def interval(self) -> np.ndarray[1, Tuple[float, float]]:
        if not self.interval_generated:
            self.interval_computed = [get_sign_change_interval_u(
                u_points) for u_points in self.u]
            self.sigma_0_generated = False

        return self.interval_computed


def main():
    grid = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    u = np.array([[0.0, -1.0, -1.0, 1.0, 3.0]])
    time = np.array([0.0])

    # print(get_sign_change_interval_u(u[0]))

    data = DataClass(grid, time, u)
    print(data.interval)
    print(data.sigma_0)
    print(data.mass_square)


if __name__ == "__main__":
    main()
