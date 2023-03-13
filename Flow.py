import numpy as np
import numpy.typing as npt
from scipy.integrate import solve_ivp, trapezoid
import pandas as pd
import h5py
import os
import timeit
import Grid
from numba import njit


class Flow():
    def __init__(self, Lambda, kir, grid, mu, T, N_Flavor=np.Inf, save_flow_flag=False, console_logging=False, number_of_observables=None, tolerance=1e-12) -> None:
        self.Lambda = Lambda
        self.kir = kir

        self.grid = grid

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

        # compute extrapolation coefficients
        self.c_0, self.c_1, self.c_2 = self.grid.compute_extrapolation_factors()

    def __str__(self):
        return str(self.grid)

    def initial_condition(self) -> npt.NDArray[np.float64]:
        match self.spatial_dimension:
            case 1:
                # for (1+1)
                intermediate = 1/np.sqrt(1+(1/self.Lambda)**2)
                return (self.grid.sigma/np.pi)*(np.arctanh(intermediate) - intermediate)
            case 2:
                # for (2+1)
                intermediate = (2+self.Lambda**2 - 2*np.sqrt(1+self.Lambda**2)
                                ) / (2*np.pi*np.sqrt(1+self.Lambda**2))
                return intermediate*self.grid.sigma

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

        self.time_start = timeit.default_timer()

        args_for_integration = (self.Lambda, self.mean_field_flag, self.console_logging_flag, self.grid.sigma, self.grid.dx_direct,
                                self.grid.dx_half, self.c_0, self.c_1, self.c_2, self.mu, self.beta, self.N_Flavor, self.spatial_dimension, self.time_start, self.kir,)
        # using the extrapolation order, to define the uband of the Jacobi matrix
        self.solution = solve_ivp(
            f, [0, tir], self.u_init, lband=1, uband=self.grid.extrapolation, method="LSODA", rtol=self.tolerance, atol=self.tolerance, t_eval=t_eval_points, args=args_for_integration)
        self.time_elapsed = (timeit.default_timer() - self.time_start)

        if self.solution.status != 0:
            raise RuntimeError(
                f'Incomplete flow for mu={self.mu}, T={self.T}, sigmaMax={self.grid[-1]}, Lambda={self.Lambda}, N={len(self.grid)}\nSolution broke at t={self.solution.t[-1]} ({self.k(self.solution.t[-1])})')
        else:
            print(
                f'integration done from 0 ({self.Lambda:.1e}) to {self.solution.t[-1]:.3f} ({self.k(self.solution.t[-1]):.1e}); time elapsed = {self.time_elapsed:.2f} seconds')

    def k(self, t):
        return self.Lambda * np.exp(-t)

    def t(self, k):
        return - np.log(k/self.Lambda)

    def get_observables_at_pos(self, j):
        t = self.solution.t[j]
        current_root = 0
        current_root_value = 0
        y = self.solution.y[:, j]
        grid = self.grid.sigma
        dx_0 = grid[1] - grid[0]
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
        filename = f'mu={self.mu}_T={self.T}_sigmaMax={self.grid.sigma[-1]}_Lambda={self.Lambda}_kir={self.kir}_nGrid={len(self.grid.sigma)}_nFlavor={self.N_Flavor}_tolerance={self.tolerance:e}_d={self.spatial_dimension}.hdf5'
        path_and_filename = os.path.join(path, filename)

        # compute observables
        observables = self.observable_array

        with h5py.File(path_and_filename, "w") as f:
            f.attrs["sigmaMax"] = self.grid.sigma[-1]
            f.attrs["NFlavor"] = self.N_Flavor
            f.attrs["Lambda"] = self.Lambda
            f.attrs["kir"] = self.kir
            f.attrs["NGrid"] = len(self.grid.sigma)
            f.attrs["mu"] = self.mu
            f.attrs["T"] = self.T
            f.attrs["computation_time"] = self.time_elapsed
            f.attrs["tolerance"] = self.tolerance
            f.attrs["spatial_dimension"] = self.spatial_dimension
            f.attrs["grid_style"] = type(self.grid).__name__
            f.attrs["extrapolation_order"] = self.grid.extrapolation
            f.create_dataset("t", data=observables["t"])
            f.create_dataset("sigma", data=observables["sigma"])
            f.create_dataset("massSquare", data=observables["massSquare"])
            f.create_dataset("first_div", data=observables["first_div"])
            f.create_dataset("third_div", data=observables["third_div"])
            f.create_dataset("pressure", data=observables["pressure"])
            f.create_dataset("k", data=observables["k"])

            if self.save_flow_flag:
                f.create_dataset("grid", data=self.grid.sigma)
                y_dset = f.create_dataset("y", (observables["y"].to_numpy().shape[0], len(
                    self.grid.sigma)))
                Q_dset = f.create_dataset("Q", (observables["y"].to_numpy().shape[0], len(
                    self.grid.sigma)))
                S_dset = f.create_dataset("S", (observables["y"].to_numpy().shape[0], len(
                    self.grid.sigma)))

                for i, elem in enumerate(observables["y"].to_numpy()):
                    y_dset[i, :] = elem[:]
                    k = observables["k"][i]
                    # Q_dset[i, :] = self.compute_diffusion(k, elem[:])
                    # S_dset[i, :] = self.S(k, self.grid.sigma)


def main():
    Lambda = 5e2
    kir = 1e-4
    n_flavor = 2
    # n_flavor = np.Inf
    mu = 0.0
    T = 0.3
    path = './test'

    # configure spatial domain
    n_grid = 1000
    sigma_max = 1000
    extrapolation_oder = 1
    grid = Grid.RescaledGeomspace(sigma_max, n_grid, extrapolation_oder)

    flow = Flow(Lambda, kir, grid, mu, T,
                n_flavor, save_flow_flag=True, console_logging=True, number_of_observables=1000, tolerance=1e-12)
    print(flow)

    for i in range(10):
        flow.compute()
    flow.get_observables_for_all_positions()
    if not (path is None):
        flow.save(path)


if __name__ == "__main__":
    main()
