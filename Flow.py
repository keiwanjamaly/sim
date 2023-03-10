import numpy as np
import numpy.typing as npt
from scipy.integrate import solve_ivp, trapezoid
from scipy.optimize import newton
import pandas as pd
import h5py
import os
import time


@np.errstate(over="ignore")
def sech(x):
    return 1/np.cosh(x)


@np.errstate(over="ignore")
def csch(x):
    return 1/np.sinh(x)


def solve_A(x_max, x_min):
    def f(A):
        return A*np.exp(x_min/A) - x_max

    return newton(f, 1.1, maxiter=1000)


class Flow():
    def __init__(self, Lambda, kir, sigma_max_1, sigma_max_2, dx_first, dx_last, mu, T, N_Flavor=np.Inf, save_flow_flag=False, console_logging=False, number_of_observables=None, tolerance=1e-12) -> None:
        self.Lambda = Lambda
        self.kir = kir

        # build non-uniform-grid
        grid_linear = np.arange(0, sigma_max_1, dx_first)
        logarithmic_start_point = 2*grid_linear[-1] - grid_linear[-2]
        dh = 10**(grid_linear[-1] - grid_linear[-2] +
                  logarithmic_start_point) - 10**logarithmic_start_point
        grid_logarithmic = [logarithmic_start_point]
        i = 0
        b = solve_b(sigma_max_2, dx_last, dh)

        N = int(np.ceil(np.log(sigma_max_2/logarithmic_start_point)/(dh*np.log(b)))) - 1
        grid_logarithmic = np.logspace(np.log(
            logarithmic_start_point)/np.log(b), np.log(sigma_max_2)/np.log(b), N, base=b)
        # factor = np.log(logarithmic_start_point)/np.log(b)
        # while grid_logarithmic[-1] < sigma_max_2:
        #     i += 1
        #     grid_logarithmic.append(
        #         b**(factor + i * dh))
        ex_l = -(grid_linear[1] - grid_linear[0])
        i += 1
        ex_r = b**(factor + i * dh)
        self.grid = np.concatenate((grid_linear, grid_logarithmic))
        self.dx_direct = np.ediff1d(
            self.grid, to_begin=(self.grid[0] - ex_l), to_end=(ex_r - self.grid[-1]))
        midpoints = np.empty(len(self.grid) + 1)
        midpoints[0] = (ex_l + self.grid[0])/1
        midpoints[1:-1] = (self.grid[1:] + self.grid[:-1])/2
        midpoints[-1] = (ex_r + self.grid[-1])/1
        self.dx_half = np.ediff1d(midpoints)

        # compute extrapolation factors
        self.c_0 = (ex_r - self.grid[-2]) * (ex_r - self.grid[-1])
        self.c_0 /= (self.grid[-3] - self.grid[-2]) * \
            (self.grid[-3] - self.grid[-1])

        self.c_1 = (ex_r - self.grid[-3]) * (ex_r - self.grid[-1])
        self.c_1 /= (self.grid[-2] - self.grid[-3]) * \
            (self.grid[-2] - self.grid[-1])

        self.c_2 = (ex_r - self.grid[-3]) * (ex_r - self.grid[-2])
        self.c_2 /= (self.grid[-1] - self.grid[-3]) * \
            (self.grid[-1] - self.grid[-2])

        self.spatial_dimension = 2
        self.u_init = self.initial_condition()
        self.mu = mu
        self.T = T
        self.beta = 1/T
        self.N_Flavor = N_Flavor
        self.mean_field_flag = np.isinf(N_Flavor)

        self.save_flow_flag = save_flow_flag
        self.console_logging_flag = console_logging
        self.print_counter = 0
        self.number_of_observables = number_of_observables
        self.tolerance = tolerance

    def __str__(self):
        return f'using {len(self.grid)} grid points'

    def initial_condition(self) -> npt.NDArray[np.float64]:
        match self.spatial_dimension:
            case 1:
                # for (1+1)
                intermediate = 1/np.sqrt(1+(1/self.Lambda)**2)
                return (self.grid/np.pi)*(np.arctanh(intermediate) - intermediate)
            case 2:
                # for (2+1)
                intermediate = (2+self.Lambda**2 - 2*np.sqrt(1+self.Lambda**2)
                                ) / (2*np.pi*np.sqrt(1+self.Lambda**2))
                return intermediate*self.grid

    def S(self, k, sigma):
        e = self.e_f(k, sigma)
        mu = self.mu
        beta = self.beta
        minus = beta * (e-mu) / 2
        plus = beta * (e+mu) / 2

        match self.spatial_dimension:
            case 1:
                # for (1+1)
                return sigma*k**3/(4*e**3*np.pi) * (e * beta * (sech(minus)**2 + sech(plus)**2)
                                                    - 2*(np.tanh(minus) + np.tanh(plus)))
            case 2:
                # for (2+1)
                return sigma*k**4/(8*e**3*np.pi) * (e * beta * (sech(minus)**2 + sech(plus)**2)
                                                    - 2*(np.tanh(minus) + np.tanh(plus)))

    def Q(self, k, ux):
        beta = self.beta
        N = self.N_Flavor
        e = self.e_b(k, ux)

        match self.spatial_dimension:
            case 1:
                # for (1+1)
                return - k**3 / (2*np.pi*e*N) * (1+2*self.n_b(beta * e))
            case 2:
                # for (2+1)
                return - k**4 / (4*np.pi*2*e*N) * (1+2*self.n_b(beta * e))

    def compute_diffusion(self, k, u):
        if self.mean_field_flag:
            return np.zeros_like(self.grid)
        else:
            extrapolation_u_right = self.c_0 * \
                u[-3] + self.c_1 * u[-2] + self.c_2 * u[-1]
            ux = np.ediff1d(
                u, to_begin=u[1], to_end=(extrapolation_u_right - u[-1])) / self.dx_direct
            Q_cal = self.Q(k, ux)
            return np.ediff1d(Q_cal)/self.dx_half

    def f(self, t, u):
        k = self.k(t)

        diffusion = self.compute_diffusion(k, u)

        if self.console_logging_flag:
            self.print_counter += 1
            if self.print_counter > 0:
                print_string = f'{t:.7f} ({self.k(t):.5e})/{self.t(self.kir):.2f}'
                self.print_counter -= 1000
                print(print_string)

        return diffusion + self.S(k, self.grid)

    def compute(self):
        tir = self.t(self.kir)
        if self.number_of_observables is None:
            t_eval_points = None
        elif self.number_of_observables == 0:
            raise RuntimeWarning(
                f'numer of observables must be larger then 1, it is currently {self.number_of_observables}')
        elif self.number_of_observables == 1:
            t_eval_points = [tir]
        else:
            t_eval_points = np.linspace(0, tir, self.number_of_observables)

        time_start = time.time()
        self.solution = solve_ivp(
            self.f, [0, tir], self.u_init, lband=1, uband=2, method="LSODA", rtol=self.tolerance, atol=self.tolerance, t_eval=t_eval_points)
        self.time_elapsed = (time.time() - time_start)

        if self.solution.status != 0:
            raise RuntimeError(
                f'Incomplete flow for mu={self.mu}, T={self.T}, sigmaMax={self.grid[-1]}, Lambda={self.Lambda}, N={len(self.grid)}\nSolution broke at t={self.solution.t[-1]} ({self.k(self.solution.t[-1])})')

    def k(self, t):
        return self.Lambda * np.exp(-t)

    def t(self, k):
        return - np.log(k/self.Lambda)

    def e_f(self, k, sigma):
        return np.sqrt(k**2 + sigma**2)

    def e_b(self, k, u_x):
        return np.sqrt(k**2 + u_x)

    @np.errstate(over="ignore")
    def n_f(self, x):
        return 1/(np.exp(x) + 1)

    @np.errstate(over="ignore")
    def n_b(self, x):
        return 1/(np.exp(x) - 1)

    def get_grid(self):
        return self.grid

    def get_observables_at_pos(self, j):
        t = self.solution.t[j]
        current_root = 0
        current_root_value = 0
        y = self.solution.y[:, j]
        grid = self.grid
        dx_0 = self.grid[1] - self.grid[0]
        massSquare = (y[1] - y[0])/dx_0

        # calculate third derivative at sigma = 0
        third_div = (y[2] - 2*y[1])/dx_0**3
        first_div = massSquare

        for i in range(1, len(y)-1):
            # find bracketed root
            if y[i] < 0 and y[i+1] > 0:
                # linear interpolate the bracketed interval and find the root
                test_root = grid[i] - y[i] * \
                    (grid[i+1] - grid[i])/(y[i+1] - y[i])

                # calculate the potential at the minimum
                rest_potential_interval = test_root * \
                    y[i] - grid[i] * y[i] + \
                    ((test_root-grid[i])**2*(y[i+1] - y[i])) / \
                    (2*(grid[i+1] - grid[i]))
                testPotential = trapezoid(
                    y[:i+1], grid[:i+1]) + rest_potential_interval
                if current_root_value > testPotential:
                    current_root_value = testPotential
                    current_root = test_root
                    # NOTE: remember that you have coosen mass square from the next interval
                    massSquare = (y[i+2] - y[i+1])/(grid[i+2] - grid[i+1])

        return {"sigma": current_root, "massSquare": massSquare, "first_div": first_div, "third_div": third_div, "pressure": - current_root_value, "t": t, "k": self.k(t), "y": y}

    def get_observables_for_all_positions(self):
        self.observable_array = pd.DataFrame([self.get_observables_at_pos(
            i) for i in range(len(self.solution.t))])

        return self.observable_array

    def save(self, path):
        # create folder, if it doesn't exist
        if not os.path.exists(path):
            os.makedirs(path)

        # generate filenames
        filename = f'mu={self.mu}_T={self.T}_sigmaMax={self.grid[-1]}_Lambda={self.Lambda}_kir={self.kir}_nGrid={len(self.grid)}_nFlavor={self.N_Flavor}_tolerance={self.tolerance:e}_d={self.spatial_dimension}.hdf5'
        path_and_filename = os.path.join(path, filename)

        # compute observables
        observables = self.observable_array

        with h5py.File(path_and_filename, "w") as f:
            f.attrs["sigmaMax"] = self.grid[-1]
            f.attrs["NFlavor"] = self.N_Flavor
            f.attrs["Lambda"] = self.Lambda
            f.attrs["kir"] = self.kir
            f.attrs["NGrid"] = len(self.grid)
            f.attrs["mu"] = self.mu
            f.attrs["T"] = self.T
            f.attrs["computation_time"] = self.time_elapsed
            f.attrs["tolerance"] = self.tolerance
            f.create_dataset("t", data=observables["t"])
            f.create_dataset("sigma", data=observables["sigma"])
            f.create_dataset("massSquare", data=observables["massSquare"])
            f.create_dataset("first_div", data=observables["first_div"])
            f.create_dataset("third_div", data=observables["third_div"])
            f.create_dataset("pressure", data=observables["pressure"])
            f.create_dataset("k", data=observables["k"])

            if self.save_flow_flag:
                f.create_dataset("grid", data=self.grid)
                y_dset = f.create_dataset("y", (observables["y"].to_numpy().shape[0], len(
                    self.grid)))
                Q_dset = f.create_dataset("Q", (observables["y"].to_numpy().shape[0], len(
                    self.grid)))
                S_dset = f.create_dataset("S", (observables["y"].to_numpy().shape[0], len(
                    self.grid)))

                for i, elem in enumerate(observables["y"].to_numpy()):
                    y_dset[i, :] = elem[:]
                    k = observables["k"][i]
                    Q_dset[i, :] = self.compute_diffusion(k, elem[:])
                    S_dset[i, :] = self.S(k, self.grid)


def main():
    Lambda = 1e3
    kir = 1e-4
    n_flavor = 2
    # n_flavor = np.Inf
    dx_first = 0.006
    sigma_max_1 = 1.2
    dx_last = 10
    sigma_max_2 = 1000
    mu = 0.0
    T = 0.3
    path = './'

    # kir, sigma_max_1, sigma_max_2, dx_first, dx_last, mu

    flow = Flow(Lambda, kir, sigma_max_1, sigma_max_2, dx_first, dx_last, mu, T,
                n_flavor, save_flow_flag=True, console_logging=True, number_of_observables=1000, tolerance=1e-12)
    print(flow)

    flow.compute()
    flow.get_observables_for_all_positions()
    if not (path is None):
        flow.save(path)


if __name__ == "__main__":
    main()
