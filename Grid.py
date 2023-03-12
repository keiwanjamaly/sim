import numpy as np
from scipy.interpolate import RegularGridInterpolator


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


class SinhGrid(Grid):
    def __init__(self, sigma_max, N) -> None:
        # create nonuniform grid
        hmin = np.arcsinh(0)
        hmax = np.arcsinh(sigma_max)
        h = np.linspace(hmin, hmax, N)
        self.sigma = np.sinh(h)

        # generate all other shit
        super().__init__(self.sigma)


class Adaptive(Grid):
    def __init__(self, sigma_max, tolerance=1e-4) -> None:
        self.init = False
        self.N_start = 11
        self.sigma_max = sigma_max
        self.sigma = np.linspace(0, self.sigma_max, self.N_start)
        self.tolerance = tolerance
        super().__init__(self.sigma)

    def save_old_parameters_and_set_new_simga(self, sigma_old, u, t, sigma_new):
        self.sigma_old = sigma_old
        self.u_old = u
        self.t_old = t
        # self.simga = np.array(sigma_new)
        super().__init__(np.array(sigma_new))
        return np.array(sigma_new)

    def remesh(self, u, t):
        # refine every interval at the beginning
        if not self.init:
            new_sigma = []
            for interval in zip(self.sigma[:-1], self.sigma[1:]):
                new_sigma.append(interval[0])
                new_sigma.append((interval[0] + interval[-1])/2)
            new_sigma.append(self.sigma[-1])
            self.sigma = self.save_old_parameters_and_set_new_simga(
                self.sigma, u, t, new_sigma)
            self.init = True
            return True

        old_interpolation = RegularGridInterpolator(
            (self.sigma_old, self.t_old), self.u_old)
        new_interpolation = RegularGridInterpolator(
            (self.sigma, t), u)

        t_test = np.sort(np.unique(np.concatenate(
            (t, self.t_old))))
        sigma_test = np.sort(np.unique(np.concatenate(
            (self.sigma, self.sigma_old))))

        X, T = np.meshgrid(sigma_test, t_test, indexing='ij')

        error = np.abs(old_interpolation((X, T)) - new_interpolation((X, T)))

        refine_points = error > self.tolerance
        # corse_points = error < self.tolerance
        final_refinement = np.empty(refine_points.shape[0], dtype=bool)
        for i in range(len(final_refinement)):
            final_refinement[i] = np.any(refine_points[i, :])

        intervals_to_be_refined = np.logical_or(
            final_refinement[:-1], final_refinement[1:])

        if not intervals_to_be_refined.any():
            return False

        new_sigma = []
        for i, refine in enumerate(intervals_to_be_refined):
            new_sigma.append(self.sigma[i])
            if refine:
                new_sigma.append((self.sigma[i] + self.sigma[i+1])/2)
        new_sigma.append(self.sigma[-1])

        self.sigma = self.save_old_parameters_and_set_new_simga(
            self.sigma, u, t, new_sigma)
        print(len(self.sigma))

        return True


def main():
    def cal_u(x, t):
        return x ** 2 * np.sin(t**2)
    grid = Adaptive(10, tolerance=1e-4)
    t = np.linspace(0, 20, 10000)

    X, T = np.meshgrid(grid.sigma, t, indexing='ij')
    u = cal_u(X, T)
    while grid.remesh(u, t):
        X, T = np.meshgrid(grid.sigma, t, indexing='ij')
        u = cal_u(X, T)

    print(np.ediff1d(grid.sigma))


if __name__ == "__main__":
    main()
