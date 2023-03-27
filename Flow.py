import numpy as np
import numpy.typing as npt
from scipy.integrate import solve_ivp, trapezoid
import pandas as pd
import h5py
import os
import timeit
import Grid
import integration_functions as int_fun
from numba import njit, objmode
import timeit


@njit(cache=True)
def compute_diffusion(k, u, mean_field_flag, Q, grid, dx_direct, dx_half, c_0, c_1, c_2, *args):
    if mean_field_flag:
        return np.zeros_like(grid)
    else:
        extrapolation_u_right = c_0 * \
            u[-3] + c_1 * u[-2] + c_2 * u[-1]
        endpoint = extrapolation_u_right - u[-1]
        # perform u_x derivative
        ux = np.empty_like(dx_direct)
        ux[0] = u[1]
        ux[1:-1] = u[1:] - u[:-1]
        ux[-1] = endpoint
        ux = ux / dx_direct
        Q_cal = Q(k, ux, *args)

        return (Q_cal[1:] - Q_cal[:-1])/dx_half


@njit(cache=True)
def f(t, u, mean_field_flag, Q, S, Lambda, grid, dx_direct, dx_half, c_0, c_1, c_2, console_logging_flag, time_start, tir, iterator, *args):
    """
    The *args conventions are as follows
    *args = (Lambda, mean_field_flag, console_logging_flag, grid, dx_direct,
      dx_half, c_0, c_1, c_2, mu, beta, N_Flavor, spatial_dimension, time_start, kir, iterations,)
    """
    k = Lambda * np.exp(-t)

    diffusion = compute_diffusion(
        k, u, mean_field_flag, Q, grid, dx_direct, dx_half, c_0, c_1, c_2, *args)

    source = S(k, grid, *args)

    if console_logging_flag:
        iterator[0] += 1
        if iterator[0] > 0:
            with objmode():
                time_elapsed = timeit.default_timer() - time_start
                print_string = '{time:.7f} ({kval:.5e})/{tirval:.2f}; time elapsed = {te:.2f} seconds'.format(
                    time=t, kval=k, tirval=tir, te=time_elapsed)
                iterator[0] -= 1000
                print(print_string, end="\r")

    return diffusion + source


class Flow():
    def __init__(self, Lambda, kir, grid, mu, T, filename=None,
                 initial_condition=None, Q=None, S=None, args=(),
                 save_flow_flag=False, console_logging=False,
                 number_of_observables=None, tolerance=1e-12, file_attributes={}):
        self.Lambda = Lambda
        self.kir = kir

        self.function_args = args

        self.grid = grid

        self.u_init = initial_condition(
            self.grid.sigma, *self.function_args)
        self.mu = mu
        self.T = T
        self.beta = 1/T
        self.S = S
        self.Q = Q
        self.mean_field_flag = Q is None

        self.save_flow_flag = save_flow_flag
        self.console_logging_flag = console_logging
        self.print_counter = 0
        self.number_of_observables = number_of_observables
        self.tolerance = tolerance

        # compute extrapolation coefficients
        self.c_0, self.c_1, self.c_2 = self.grid.compute_extrapolation_factors()

        self.filename = filename
        self.file_attributes = file_attributes

    def __str__(self):
        return str(self.grid)

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

        # args_for_integration = (self.Lambda, self.mean_field_flag, self.console_logging_flag, self.grid.sigma, self.grid.dx_direct,
        #                         self.grid.dx_half, self.c_0, self.c_1, self.c_2, self.mu, self.beta, self.N_Flavor, self.spatial_dimension, self.time_start, tir, np.array([0]),)
        # using the extrapolation order, to define the uband of the Jacobi matrix
        # Q, S, Lambda, grid, dx_direct, dx_half, c_0, c_1, c_2, console_logging_flag, time_start, tir, iterator, *args
        self.solution = solve_ivp(
            f, [0, tir], self.u_init, lband=1, uband=self.grid.extrapolation, method="LSODA", rtol=self.tolerance, atol=self.tolerance, t_eval=t_eval_points,
            args=(self.mean_field_flag, self.Q, self.S, self.Lambda, self.grid.sigma, self.grid.dx_direct, self.grid.dx_half, self.c_0, self.c_1, self.c_2, self.console_logging_flag, self.time_start,
                  tir, np.array([0]), *self.function_args, ))
        self.time_elapsed = (timeit.default_timer() - self.time_start)

        if self.solution.status != 0:
            raise RuntimeError(
                f'Incomplete flow for mu={self.mu}, T={self.T}, sigmaMax={self.grid[-1]}, Lambda={self.Lambda}, N={len(self.grid.sigma)}\nSolution broke at t={self.solution.t[-1]} ({self.k(self.solution.t[-1])})')
        else:
            print(
                f'    Integration done from 0 ({self.Lambda:.1e}) to {self.solution.t[-1]:.3f} ({self.k(self.solution.t[-1]):.1e}); time elapsed = {self.time_elapsed:.2f} seconds')

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

        k = self.k(t)
        return {"sigma": current_root, "massSquare": massSquare, "first_div": first_div, "third_div": third_div,
                "pressure": - current_root_value, "t": t, "k": k, "y": y,
                "diffusion": compute_diffusion(k, y, self.mean_field_flag, self.Q, self.grid.sigma,
                                               self.grid.dx_direct, self.grid.dx_half, self.c_0, self.c_1,
                                               self.c_2, *self.function_args),
                "source": self.S(
                    k, self.grid.sigma, *self.function_args)}

    def get_observables_for_all_positions(self):
        self.observable_array = pd.DataFrame([self.get_observables_at_pos(
            i) for i in range(len(self.solution.t))])

        return self.observable_array

    def save(self, path):
        # create folder, if it doesn't exist
        if not os.path.exists(path):
            os.makedirs(path)

        # generate filenames
        path_and_filename = os.path.join(path, self.filename)

        # compute observables
        observables = self.observable_array

        with h5py.File(path_and_filename, "w") as f:
            f.attrs["sigmaMax"] = self.grid.sigma[-1]
            f.attrs["Lambda"] = self.Lambda
            f.attrs["kir"] = self.kir
            f.attrs["NGrid"] = len(self.grid.sigma)
            f.attrs["computation_time"] = self.time_elapsed
            f.attrs["tolerance"] = self.tolerance
            f.attrs["grid_style"] = type(self.grid).__name__
            f.attrs["extrapolation_order"] = self.grid.extrapolation
            # TODO: include manual storage of additional attrs
            for key, value in self.file_attributes.items():
                f.attrs[key] = value

            # create datasets and store them
            f.create_dataset("t", data=observables["t"])
            f.create_dataset("sigma", data=observables["sigma"])
            f.create_dataset("massSquare", data=observables["massSquare"])
            f.create_dataset("first_div", data=observables["first_div"])
            f.create_dataset("third_div", data=observables["third_div"])
            f.create_dataset("pressure", data=observables["pressure"])
            f.create_dataset("k", data=observables["k"])

            # create additional datasets and store them
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
                    Q_dset[i, :] = observables["diffusion"][i]
                    S_dset[i, :] = observables["source"][i]

        return path_and_filename


def main():
    spatial_dimension = 2
    Lambda = 100
    kir = 1e-4
    n_flavor = 2
    # n_flavor = np.Inf
    mu = 0.0
    T = 0.1
    path = './'
    tolerance = 1e-12

    # configure spatial domain
    n_grid = 1000
    sigma_max = 6
    extrapolation_oder = 1
    grid = Grid.UniformGrid(sigma_max, n_grid, extrapolation_oder)

    args = int_fun.generate_args(spatial_dimension, Lambda, mu, T, n_flavor)

    storage_dict = int_fun.storage_dictionary(*args)

    filename = int_fun.generate_filename(
        mu, T, sigma_max, n_grid, kir, tolerance, *args)

    flow = Flow(Lambda, kir, grid, mu, T, filename=filename, initial_condition=int_fun.initial_condition, Q=int_fun.Q, S=int_fun.S,
                args=args, save_flow_flag=True, console_logging=True, number_of_observables=100, tolerance=tolerance, file_attributes=storage_dict)
    print(flow)

    flow.compute()
    flow.get_observables_for_all_positions()
    if not (path is None):
        flow.save(path)


if __name__ == "__main__":
    main()
