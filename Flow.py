import numpy as np
import numpy.typing as npt
from scipy.integrate import solve_ivp, trapezoid
import pandas as pd
import h5py
import os
import argparse


@np.errstate(over="ignore")
def sech(x):
    return 1/np.cosh(x)


@np.errstate(over="ignore")
def csch(x):
    return 1/np.sinh(x)


class Flow():
    def __init__(self, Lambda, kir, sigma_max, n_grid, mu, T, N_Flavor=np.Inf) -> None:
        self.__Lambda = Lambda
        self.__kir = kir
        self.__grid = np.linspace(0, sigma_max, n_grid)
        self.__dx = self.__grid[1] - self.__grid[0]
        self.__u_init = self.initial_condition()
        self.__mu = mu
        self.__T = T
        self.__beta = 1/T
        self.__N_Flavor = N_Flavor
        self.__mean_field_flag = np.isinf(N_Flavor)

        self.__print_counter = 0

    def initial_condition(self) -> npt.NDArray[np.float64]:
        intermediate = 1/np.sqrt(1+(1/self.__Lambda)**2)
        return (self.__grid/np.pi)*(np.arctanh(intermediate) - intermediate)

    def S(self, t, sigma):
        k = self.k(t)
        e = self.e_f(k, sigma)
        mu = self.__mu
        beta = self.__beta
        minus = beta * (e-mu) / 2
        plus = beta * (e+mu) / 2
        return sigma*k**3/(4*e**3*np.pi) * (e * beta * (sech(minus)**2 + sech(plus)**2)
                                            - 2*(np.tanh(minus) + np.tanh(plus)))

    def Q(self, t, ux):
        k = self.k(t)
        beta = self.__beta
        N = self.__N_Flavor
        return - k**3 / (2*np.pi*self.e_b(k, ux)*N) * (1+2*self.n_b(beta * self.e_b(k, ux)))

    def f(self, t, u):
        self.__print_counter += 1
        if self.__print_counter % 300 == 0:
            print(t, "/", self.t(self.__kir))

        if self.__mean_field_flag:
            return self.S(t, self.__grid)
        else:
            u_x = np.ediff1d(u, to_begin=u[1], to_end=u[-1]-u[-2]) / self.__dx
            return (self.Q(t, u_x[1:]) - self.Q(t, u_x[:-1]))/self.__dx + self.S(t, self.__grid)

    def compute(self):
        tir = self.t(self.__kir)

        self.__solution = solve_ivp(
            self.f, [0, tir], self.__u_init, lband=1, uband=1, method="LSODA", rtol=1e-10, atol=1e-10)

    def k(self, t):
        return self.__Lambda * np.exp(-t)

    def t(self, k):
        return - np.log(k/self.__Lambda)

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
        return self.__grid

    def get_observables_at_pos(self, j):
        t = self.__solution.t[j]
        current_root = 0
        current_root_value = 0
        y = self.__solution.y[:, j]
        grid = self.__grid
        massSquare = (y[1] - y[0])/self.__dx

        # calculate third derivative at sigma = 0
        third_div = (y[2] - 2*y[1])/self.__dx**3
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
                    massSquare = (y[i+1] - y[i])/(grid[i+1] - grid[i])

        return {"sigma": current_root, "massSquare": massSquare, "first_div": first_div, "third_div": third_div, "pressure": - current_root_value, "t": t, "k": self.k(t), "y": y}

    def get_observables_for_all_positions(self):
        self.observable_array = pd.DataFrame([self.get_observables_at_pos(
            i) for i in range(len(self.__solution.t))])

        return self.observable_array

    def save(self, path):
        # create folder, if it doesn't exist
        if not os.path.exists(path):
            os.makedirs(path)

        # generate filenames
        filename = f'mu={self.__mu}_T={self.__T}_sigmaMax={self.__grid[-1]}_Lambda={self.__Lambda}_kir={self.__kir}_nGrid={len(self.__grid)}_nFlavor={self.__N_Flavor}.hdf5'
        path_and_filename = os.path.join(path, filename)

        # compute observables
        observables = self.observable_array

        with h5py.File(path_and_filename, "w") as f:
            f.attrs["sigmaMax"] = self.__grid[-1]
            f.attrs["NFlavor"] = self.__N_Flavor
            f.attrs["Lambda"] = self.__Lambda
            f.attrs["kir"] = self.__kir
            f.attrs["NGrid"] = len(self.__grid)
            f.attrs["mu"] = self.__mu
            f.attrs["T"] = self.__T
            f.create_dataset("t", data=observables["t"])
            f.create_dataset("sigma", data=observables["sigma"])
            f.create_dataset("massSquare", data=observables["massSquare"])
            f.create_dataset("first_div", data=observables["first_div"])
            f.create_dataset("third_div", data=observables["third_div"])
            f.create_dataset("pressure", data=observables["pressure"])
            f.create_dataset("k", data=observables["k"])
            # f.create_dataset("grid", data=self.__grid)
            # y_dset = f.create_dataset("y", (observables["y"].to_numpy().shape[0], len(
            #     self.__grid)))

            # for i, elem in enumerate(observables["y"].to_numpy()):
            #     y_dset[i, :] = elem[:]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-N", type=int,
                        help="flavor number")
    parser.add_argument("-L", type=float,
                        help="ultra-violet cutoff")
    parser.add_argument("-kir", type=float,
                        help="infrared cutoff")
    parser.add_argument("-sigmaMax", type=float,
                        help="maximum of spatial computation domain")
    parser.add_argument("-grid", type=int,
                        help="number of grid points")
    parser.add_argument("-mu", type=float,
                        help="chemical potential")
    parser.add_argument("-T", type=float,
                        help="Temperature")
    parser.add_argument("-o", type=str,
                        help="output file path")
    args = parser.parse_args()
    Lambda = args.L
    kir = args.kir
    n_flavor = args.N if args.N != -1 else np.Inf
    sigma_max = args.sigmaMax
    mu = args.mu
    T = args.T
    n_grid = args.grid

    path = args.o

    flow = Flow(Lambda, kir, sigma_max, n_grid, mu, T, n_flavor)
    # flow = Flow(1e5, 1e-4, 6, 1000, 0.4, 0.0125, 2)

    flow.compute()
    flow.get_observables_for_all_positions()
    if not (path is None):
        flow.save(path)


if __name__ == "__main__":
    main()
