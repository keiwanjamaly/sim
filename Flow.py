import numpy as np
import numpy.typing as npt
from scipy.integrate import solve_ivp, trapezoid
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


class Flow():
    def __init__(self, Lambda, kir, sigma_max, n_grid, mu, T, N_Flavor=np.Inf, save_flow_flag=False, console_logging=False, number_of_observables=None, tolerance=1e-12) -> None:
        self.Lambda = Lambda
        self.kir = kir
        self.grid_u, self.dx = np.linspace(0,
                                           sigma_max, n_grid, retstep=True)
        self.grid_log_mass = self.grid_u[:-1] + self.dx/2
        self.spatial_dimension = 1
        self.log_mass_init = np.full(
            self.grid_log_mass.shape, np.log(self.Lambda**2 + self.initial_condition()))
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

        self.k_pow = self.spatial_dimension + 2
        self.diff_factor = 4 if self.spatial_dimension == 2 else 1

    def initial_condition(self) -> npt.NDArray[np.float64]:
        match self.spatial_dimension:
            case 1:
                # for (1+1)
                intermediate = 1/np.sqrt(1+(1/self.Lambda)**2)
                return (1/np.pi)*(np.arctanh(intermediate) - intermediate)
            case 2:
                # for (2+1)
                intermediate = (2+self.Lambda**2 - 2*np.sqrt(1+self.Lambda**2)
                                ) / (2*np.pi*np.sqrt(1+self.Lambda**2))
                return intermediate

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

    def Q(self, k, log_mass):
        beta = self.beta
        N = self.N_Flavor
        e = self.e_b(k, log_mass)

        match self.spatial_dimension:
            case 1:
                # for (1+1)
                return - k**3 / (2*np.pi*e*N) * (1+2*self.n_b(beta * e))
            case 2:
                # for (2+1)
                return - k**4 / (4*np.pi*2*e*N) * (1+2*self.n_b(beta * e))

    def f(self, t, log_mass):
        k = self.k(t)
        if self.mean_field_flag:
            diffusion = np.zeros_like(self.grid_u)
        else:
            Q_cal = self.Q(k, log_mass)
            diffusion = np.ediff1d(Q_cal, to_begin=0, to_end=0)/self.dx

        W = self.S(k, self.grid_u) + diffusion

        if self.console_logging_flag:
            self.print_counter += 1
            if self.print_counter > 0:
                print_string = f'{t:.7f} ({self.k(t):.5e})/{self.t(self.kir):.2f}'
                self.print_counter -= 1000
                print(print_string)

        result = np.exp(-log_mass)*(np.ediff1d(W)/self.dx - 2*k**2)
        return result

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
            self.f, [0, tir], self.log_mass_init, lband=1, uband=1, method="LSODA", t_eval=t_eval_points, atol=self.tolerance, rtol=self.tolerance)
        self.time_elapsed = (time.time() - time_start)

        if self.solution.status != 0:
            raise RuntimeError(
                f'Incomplete flow for mu={self.mu}, T={self.T}, sigmaMax={self.grid_u[-1]}, Lambda={self.Lambda}, N={len(self.grid_u)}')

    def k(self, t):
        return self.Lambda * np.exp(-t)

    def t(self, k):
        return - np.log(k/self.Lambda)

    def e_f(self, k, sigma):
        return np.sqrt(k**2 + sigma**2)

    def e_b(self, k, log_mass):
        return np.exp(log_mass/2)

    @np.errstate(over="ignore")
    def n_f(self, x):
        return 1/(np.exp(x) + 1)

    @np.errstate(over="ignore")
    def n_b(self, x):
        return 1/(np.exp(x) - 1)

    def get_observables_at_pos(self, j):
        t = self.solution.t[j]
        current_root = 0
        current_root_value = 0
        log_mass = self.solution.y[:, j]

        # calculate u
        enlarged_log_mass = np.empty(len(self.grid_u)+len(self.grid_log_mass))
        enlarged_log_mass[1:-1:2] = log_mass[:]
        enlarged_log_mass[2:-2:2] = (log_mass[:-1] + log_mass[1:])/2
        enlarged_log_mass[0] = 2*log_mass[0] - log_mass[1]
        enlarged_log_mass[-1] = 2*log_mass[-1] - log_mass[-2]
        y = np.empty(enlarged_log_mass.shape[0]+1)
        y[0] = 0
        y[1:] = np.cumsum(np.exp(enlarged_log_mass) - self.k(t)**2)*self.dx/2
        y = y[::2]

        grid = self.grid_u
        massSquare = (y[1] - y[0])/self.dx

        # calculate third derivative at sigma = 0
        third_div = (y[2] - 2*y[1])/self.dx**3
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
        filename = f'mu={self.mu}_T={self.T}_sigmaMax={self.grid_u[-1]}_Lambda={self.Lambda}_kir={self.kir}_nGrid={len(self.grid_u)}_nFlavor={self.N_Flavor}_tolerance={self.tolerance:e}.hdf5'
        path_and_filename = os.path.join(path, filename)

        # compute observables
        observables = self.observable_array

        with h5py.File(path_and_filename, "w") as f:
            f.attrs["sigmaMax"] = self.grid_u[-1]
            f.attrs["NFlavor"] = self.N_Flavor
            f.attrs["Lambda"] = self.Lambda
            f.attrs["kir"] = self.kir
            f.attrs["NGrid"] = len(self.grid_u)
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
                f.create_dataset("grid", data=self.grid_u)
                y_dset = f.create_dataset("y", (observables["y"].to_numpy().shape[0], len(
                    self.grid_u)))

                for i, elem in enumerate(observables["y"].to_numpy()):
                    y_dset[i, :] = elem[:]


def main():
    Lambda = 1e5
    kir = 1e-4
    n_flavor = 2
    # n_flavor = np.Inf
    sigma_max = 6.0
    mu = 0.1
    T = 0.1
    n_grid = 1000
    path = './'

    flow = Flow(Lambda, kir, sigma_max, n_grid, mu, T,
                n_flavor, save_flow_flag=True, console_logging=True, number_of_observables=1000, tolerance=5e-14)

    flow.compute()
    flow.get_observables_for_all_positions()
    if not (path is None):
        flow.save(path)


if __name__ == "__main__":
    main()
