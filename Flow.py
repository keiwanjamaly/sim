import numpy as np
import numpy.typing as npt
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from numba import njit
import sys
import pickle


@njit
def sech(x):
    return 1/np.cosh(x)


@njit
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
        return - k**3 / (2*np.pi*self.e_b(k, ux)*N) * (1+self.n_b(beta * self.e_b(k, ux)))

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
        self.__solution = solve_ivp(
            self.f, [0, self.t(self.__kir)], self.__u_init, lband=1, uband=1, method="LSODA", rtol=1e-10, atol=1e-10)
        print(self.__solution)
        pickle.dump(self.__solution, open('test.pkl', 'wb'))
        print(sys.getsizeof(self.__solution.y))

        return self.__solution

    def k(self, t):
        return self.__Lambda * np.exp(-t)

    def t(self, k):
        return - np.log(k/self.__Lambda)

    def e_f(self, k, sigma):
        return np.sqrt(k**2 + sigma**2)

    def e_b(self, k, u_x):
        return np.sqrt(k**2 + u_x)

    def n_f(self, x):
        return 1/(np.exp(x) + 1)

    def n_b(self, x):
        return 1/(np.exp(x) - 1)

    def get_grid(self):
        return self.__grid


def main():
    simple_heat = Flow(1e5, 1e-4, 6, 1000, 0.1, 0.01, 2)
    solution = simple_heat.compute()
    grid = simple_heat.get_grid()
    # plt.plot(grid, solution.y[:, -1])
    # plt.show()


if __name__ == "__main__":
    main()
