import numpy as np
import numpy.typing as npt
from scipy.integrate import solve_ivp, trapezoid
import pandas as pd
import h5py
import os


@np.errstate(over="ignore")
def sech(x):
    return 1/np.cosh(x)


@np.errstate(over="ignore")
def csch(x):
    return 1/np.sinh(x)


class Flow():
    def __init__(self, Lambda, kir, sigma_max, n_grid, mu, T, N_Flavor=np.Inf, save_flow_flag=False) -> None:
        self.Lambda = Lambda
        self.kir = kir
        self.grid = np.linspace(0, sigma_max, n_grid)
        self.dx = self.grid[1] - self.grid[0]
        self.spatial_dimension = 2
        self.u_init = self.initial_condition()
        self.mu = mu
        self.T = T
        self.beta = 1/T
        self.N_Flavor = N_Flavor
        self.mean_field_flag = np.isinf(N_Flavor)

        self.save_flow_flag = save_flow_flag
        self.print_counter = 0

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

    def S(self, t, sigma):
        k = self.k(t)
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

    def Q(self, t, ux):
        k = self.k(t)
        beta = self.beta
        N = self.N_Flavor

        match self.spatial_dimension:
            case 1:
                # for (1+1)
                return - k**3 / (2*np.pi*self.e_b(k, ux)*N) * (1+2*self.n_b(beta * self.e_b(k, ux)))
            case 2:
                # for (2+1)
                return - k**4 / (4*np.pi*2*self.e_b(k, ux)*N) * (1+2*self.n_b(beta * self.e_b(k, ux)))

    def f(self, t, u):
        if self.mean_field_flag:
            diffusion = np.zeros_like(self.grid)
        else:
            u_x = np.ediff1d(u, to_begin=u[1], to_end=u[-1]-u[-2]) / self.dx
            Q_cal = self.Q(t, u_x)
            diffusion = np.ediff1d(Q_cal)/self.dx

        self.print_counter += 1
        if self.print_counter > 0:
            print_string = f'{t:.7f} ({self.k(t):.5e})/{self.t(self.kir):.2f}'
            self.print_counter -= 1000
            print(print_string)

        return diffusion + self.S(t, self.grid)

    def compute(self):
        tir = self.t(self.kir)
        t_eval_points = np.linspace(0, tir, 1000)

        self.solution = solve_ivp(
            self.f, [0, tir], self.u_init, lband=1, uband=1, method="LSODA", rtol=1e-12, atol=1e-12, t_eval=t_eval_points)

        if self.solution.status != 0:
            raise RuntimeError(
                f'Incomplete flow for mu={self.mu}, T={self.T}, sigmaMax={self.grid[-1]}, Lambda={self.Lambda}, N={len(self.grid)}')

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
                    # TODO: remember that you have coosen a different mass square
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
        filename = f'mu={self.mu}_T={self.T}_sigmaMax={self.grid[-1]}_Lambda={self.Lambda}_kir={self.kir}_nGrid={len(self.grid)}_nFlavor={self.N_Flavor}.hdf5'
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

                for i, elem in enumerate(observables["y"].to_numpy()):
                    y_dset[i, :] = elem[:]


def main():
    Lambda = 1000
    kir = 1e-3
    n_flavor = 20
    # n_flavor = np.Inf
    sigma_max = 6.0
    mu = 0.4
    T = 0.0125
    n_grid = 1000
    path = './'

    flow = Flow(Lambda, kir, sigma_max, n_grid, mu, T, n_flavor, True)

    flow.compute()
    flow.get_observables_for_all_positions()
    if not (path is None):
        flow.save(path)


if __name__ == "__main__":
    main()
