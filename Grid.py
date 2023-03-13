import numpy as np


class Grid():
    def __init__(self, sigma) -> None:
        ex_l = -(sigma[1] - sigma[0])
        ex_r = sigma[-3] - 3 * sigma[-2] + 3*sigma[-1]
        self.dx_direct = np.ediff1d(
            sigma, to_begin=(sigma[0] - ex_l), to_end=(ex_r - sigma[-1]))
        midpoints = np.empty(len(sigma) + 1)
        midpoints[0] = (ex_l + sigma[0])/1
        midpoints[1:-1] = (sigma[1:] + sigma[:-1])/2
        midpoints[-1] = (ex_r + sigma[-1])/1
        self.dx_half = np.ediff1d(midpoints)

        # compute extrapolation factors
        self.c_0 = (ex_r - sigma[-2]) * (ex_r - sigma[-1])
        self.c_0 /= (sigma[-3] - sigma[-2]) * \
            (sigma[-3] - sigma[-1])

        self.c_1 = (ex_r - sigma[-3]) * (ex_r - sigma[-1])
        self.c_1 /= (sigma[-2] - sigma[-3]) * \
            (sigma[-2] - sigma[-1])

        self.c_2 = (ex_r - sigma[-3]) * (ex_r - sigma[-2])
        self.c_2 /= (sigma[-1] - sigma[-3]) * \
            (sigma[-1] - sigma[-2])

    def __str__(self):
        sigma_smaller_one = np.ediff1d(self.sigma[self.sigma <= 1])
        return f'Grid from 0 to {self.sigma[-1]:.2f}, with {len(self.sigma)} grid points.\n' \
            f'In the intervall [0, 1] the smalles (largest) spacing between gridpoints is {np.min(sigma_smaller_one):.3f} ({np.max(sigma_smaller_one):.3f})'


class UniformGrid(Grid):
    def __init__(self, sigma_max, N) -> None:
        self.sigma = np.linspace(0, sigma_max, N)
        super().__init__(self.sigma)


class RescaledGeomspace(Grid):
    def __init__(self, sigma_max, N) -> None:
        upper = np.log10(sigma_max + 1)
        h = np.linspace(0, upper, N)
        self.sigma = 10**h - 1
        super().__init__(self.sigma)


class SinhGrid(Grid):
    def __init__(self, sigma_max, N) -> None:
        # create nonuniform grid
        hmin = np.arcsinh(0)
        hmax = np.arcsinh(sigma_max)
        h = np.linspace(hmin, hmax, N)
        self.sigma = np.sinh(h)

        # generate all other shit
        super().__init__(self.sigma)
