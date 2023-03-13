import numpy as np


class Grid():
    def __init__(self, sigma, extrapolation) -> None:
        self.ex_l = -(sigma[1] - sigma[0])
        self.ex_r = sigma[-3] - 3 * sigma[-2] + 3*sigma[-1]
        self.dx_direct = np.ediff1d(
            sigma, to_begin=(sigma[0] - self.ex_l), to_end=(self.ex_r - sigma[-1]))
        midpoints = np.empty(len(sigma) + 1)
        midpoints[0] = (self.ex_l + sigma[0])/1
        midpoints[1:-1] = (sigma[1:] + sigma[:-1])/2
        midpoints[-1] = (self.ex_r + sigma[-1])/1
        self.dx_half = np.ediff1d(midpoints)

        self.extrapolation = extrapolation

    def compute_extrapolation_factors(self):
        match self.extrapolation:
            case 1:
                c_0 = 0

                c_1 = (self.ex_r - self.sigma[-1])
                c_1 /= (self.sigma[-2] - self.sigma[-1])

                c_2 = (self.ex_r - self.sigma[-2])
                c_2 /= (self.sigma[-1] - self.sigma[-2])

                return [c_0, c_1, c_2]
            case 2:
                c_0 = (
                    self.ex_r - self.sigma[-2]) * (self.ex_r - self.sigma[-1])
                c_0 /= (self.sigma[-3] - self.sigma[-2]) * \
                    (self.sigma[-3] - self.sigma[-1])

                c_1 = (
                    self.ex_r - self.sigma[-3]) * (self.ex_r - self.sigma[-1])
                c_1 /= (self.sigma[-2] - self.sigma[-3]) * \
                    (self.sigma[-2] - self.sigma[-1])

                c_2 = (
                    self.ex_r - self.sigma[-3]) * (self.ex_r - self.sigma[-2])
                c_2 /= (self.sigma[-1] - self.sigma[-3]) * \
                    (self.sigma[-1] - self.sigma[-2])

                return [c_0, c_1, c_2]

    def __str__(self):
        sigma_smaller_one = np.ediff1d(self.sigma[self.sigma <= 1])
        return f'Grid from 0 to {self.sigma[-1]:.2f}, with {len(self.sigma)} grid points. Using extrapolation on right boundary of order {self.extrapolation}.\n' \
            f'In the intervall [0, 1] the smalles (largest) spacing between gridpoints is {np.min(sigma_smaller_one):.3f} ({np.max(sigma_smaller_one):.3f})'


class UniformGrid(Grid):
    def __init__(self, sigma_max, N, extrapolation) -> None:
        self.sigma = np.linspace(0, sigma_max, N)
        super().__init__(self.sigma, extrapolation)


class RescaledGeomspace(Grid):
    def __init__(self, sigma_max, N, extrapolation) -> None:
        upper = np.log10(sigma_max + 1)
        h = np.linspace(0, upper, N)
        self.sigma = 10**h - 1
        super().__init__(self.sigma, extrapolation)


class SinhGrid(Grid):
    def __init__(self, sigma_max, N, extrapolation) -> None:
        # create nonuniform grid
        hmin = np.arcsinh(0)
        hmax = np.arcsinh(sigma_max)
        h = np.linspace(hmin, hmax, N)
        self.sigma = np.sinh(h)

        # generate all other shit
        super().__init__(self.sigma, extrapolation)
