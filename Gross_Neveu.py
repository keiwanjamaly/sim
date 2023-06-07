from c_Grid import Grid_Interface
from compute_observables import DataClass
from c_Gross_Neveu import Gross_Neveu_Interface, Gross_Neveu_Interface_Vacuum
from c_Computation import Computation_Interface
from c_Computation_Data import Computation_Data_Interface
from c_Return_Data import Return_Data_Interface
from ctypes import CDLL
import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple
from math import sqrt, atanh
# from mpmath import polylog
import scipy


def expit(x):
    return scipy.special.expit(-x)


def log_expit(x):
    return scipy.special.log1p(np.exp(x))
    # return -scipy.special.log_expit(-x)


class GN_1p1():
    d = 1
    d_gamma = 2
    tolerances = 1e-10

    def __init__(self, sigma_max, Lambda, kir, N_Grid, samples, mu, T,
                 N_Flavor=np.infty, h=1, one_over_g2=None, sigma_0=1) -> None:

        if T == 0:
            if mu != 0:
                raise RuntimeError(
                    f'There is no implementation for T=0 and mu != 0, mu is {mu}')
            self.lib = CDLL("./build/src/libgross_neveu_vacuum.dylib")
        else:
            self.lib = CDLL("./build/src/libgross_neveu.dylib")

        grid_points = np.linspace(0, sigma_max, N_Grid)
        grid = Grid_Interface(self.lib, grid_points)
        self.return_data = Return_Data_Interface(
            self.lib, samples, grid.pointer)

        # calculate 1/g**2 if not set by user
        if (one_over_g2 is None):
            self.one_over_g2 = self.calculate_one_g2(h, sigma_0, Lambda)
        else:
            self.one_over_g2 = one_over_g2

        if T == 0:
            physics_data = Gross_Neveu_Interface_Vacuum(
                self.lib, h, self.one_over_g2, self.dimension(), self.dimension_gamma(), N_Flavor)
        else:
            physics_data = Gross_Neveu_Interface(
                self.lib, h, self.one_over_g2, self.dimension(),
                self.dimension_gamma(), mu, T, N_Flavor)

        computation_data = Computation_Data_Interface(
            self.lib, Lambda, kir, grid.pointer, physics_data.pointer,
            self.tolerances)

        Computation_Interface(
            self.lib, computation_data.pointer, self.return_data.pointer)

    @classmethod
    def calculate_one_g2(cls, h, sigma_0, Lambda) -> float:
        tmp = 1/sqrt(1.0 + (h * sigma_0 / Lambda)**2)
        return cls.d_gamma * (1 / (2 * np.pi)) * (atanh(tmp) - tmp)

    @classmethod
    def dimension(cls):
        return cls.d

    @classmethod
    def dimension_gamma(cls):
        return cls.d_gamma


class GN_2p1(GN_1p1):
    d = 2
    d_gamma = 4
    tolerances = 1e-12

    @classmethod
    def calculate_one_g2(cls, h, sigma_0, Lambda) -> float:
        tmp = sqrt((h * sigma_0)**2 + Lambda**2)
        tmp = cls.d_gamma * (Lambda**2 + 2 * h * sigma_0 * (h * sigma_0 - tmp)) \
            / (8 * np.pi * tmp)
        return tmp


def run_1_1_tests():
    import json
    fig, ((ax1_1, ax1_2, ax1_3), (ax2_1, ax2_2, ax2_3)) = plt.subplots(2, 3)
    fig.tight_layout()

    sigma_max = 6.0
    Lambda = 1e5
    kir = 1e-4
    N_Grid = 1000
    samples = 100

    def compute_1_1_sol(mu: float, T: float, N_Flavor: int) -> Tuple[float, float]:
        sol = GN_1p1(sigma_max, Lambda, kir, N_Grid,
                     samples, mu, T, N_Flavor=N_Flavor)
        y = sol.return_data.solution
        x = sol.return_data.grid
        time = sol.return_data.time
        return x, y, time

    # compare MF results
    mu = 0.1
    T = 0.1
    N_Flavor = np.infty
    x, y, time = compute_1_1_sol(mu, T, N_Flavor)
    data_class_MF = DataClass(x, time, y)

    # test mean field computation for 1+1 dimension
    with open('1+1_test_files/flow_MF_T=0.1,mu=0.1.json') as user_file:
        parsed_json = json.load(user_file)
        x_ref = np.array(parsed_json["graphs"][-1]["x"])
        y_ref = np.array(parsed_json["graphs"][-1]["y"])

    if (x_ref != x).all():
        raise RuntimeWarning("For MF the grids do not match")
    ax1_1.plot(x, y[-1], label="computation")
    ax1_1.plot(x_ref, y_ref, label="reference", linestyle="dashed")
    ax1_1.legend()
    ax1_1.grid()
    ax1_1.set_xlabel(r'$\sigma$')
    difference_y = np.abs(y[-1] - y_ref)
    rel_error_y = difference_y / np.abs(y[-1] + y_ref)

    ax1_2.plot(x[:-1], difference_y[:-1], label="abs_error")
    ax1_2.plot(x[:-1], rel_error_y[:-1], label="rel_error")
    ax1_2.legend()
    ax1_2.grid()
    ax1_2.set_xlabel(r'$\sigma$')

    sigma_array = data_class_MF.sigma_0
    ax1_3.plot(time, sigma_array)
    ax1_3.grid()
    ax1_3.set_xlabel(r'$t$')
    ax1_3.set_ylabel(r'$\sigma_0$')

    # compare N_Flavor = 2 results
    mu = 0.1
    T = 0.1
    N_Flavor = 2
    x, y, time = compute_1_1_sol(mu, T, N_Flavor)
    data_class_MF = DataClass(x, time, y)

    # test N_Flavor = 2 for 1+1 dimension
    with open('1+1_test_files/flow_N=2,T=0.1,mu=0.1.json') as user_file:
        parsed_json = json.load(user_file)
        x_ref = np.array(parsed_json["graphs"][-1]["x"])
        y_ref = np.array(parsed_json["graphs"][-1]["y"])

    if (x_ref != x).all():
        raise RuntimeWarning("For N_flavor=2 the grids do not match")
    ax2_1.plot(x, y[-1], label="computation")
    ax2_1.plot(x_ref, y_ref, label="reference", linestyle="dashed")
    ax2_1.legend()
    ax2_1.grid()
    ax2_1.set_xlabel(r'$\sigma$')
    difference_y = np.abs(y[-1] - y_ref)
    rel_error_y = difference_y / np.abs(y[-1] + y_ref)

    ax2_2.plot(x[:-1], difference_y[:-1], label="abs_error")
    ax2_2.plot(x[:-1], rel_error_y[:-1], label="rel_error")
    ax2_2.legend()
    ax2_2.grid()
    ax2_2.set_xlabel(r'$\sigma$')

    sigma_array = data_class_MF.sigma_0
    ax2_3.plot(time, sigma_array)
    ax2_3.grid()
    ax2_3.set_xlabel(r'$t$')
    ax2_3.set_ylabel(r'$\sigma_0$')

    plt.show()


def run_2_1_tests():
    def analytic_vacuum_solution(sigma: float, Lambda: float, h: float, sigma_0: float):
        tmp = np.sqrt(Lambda**2 + (h*sigma_0)**2)
        d_gamma = GN_2p1.dimension_gamma()
        e_vac = np.sqrt(Lambda**2 + h**2 * sigma_0**2)
        e = np.sqrt(Lambda**2 + h**2 * sigma**2)
        tmp = -6*e + 6*e_vac + 6*h * \
            (sigma - sigma_0) + 3*Lambda**2*(1/e - 1/e_vac)
        return (d_gamma*h**2*sigma*tmp)/(24*np.pi)

    # def expression2(h, sigma, beta, mu, sigma0, capital_lambda):
    # def analythic_mean_field_solution(sigma: float, Lambda: float, h: float, sigma_0: float, mu: float, T: float):
    #     beta = 1/T
    #     term1 = (1/(2 * np.pi)) * h
    #     term2 = -((h * sigma * (Lambda**2 + 2 * h * sigma * (h * sigma - np.sqrt(
    #         Lambda**2 + h**2 * sigma**2)))) / np.sqrt(Lambda**2 + h**2 * sigma**2))
    #     term3 = (h * sigma * (Lambda**2 + 2 * h * sigma_0 * (h * sigma_0 - np.sqrt(
    #         Lambda**2 + h**2 * sigma_0**2)))) / np.sqrt(Lambda**2 + h**2 * sigma_0**2)
    #     term4 = (2 * h * sigma * np.log(1 +
    #              np.exp(beta * (mu - h * sigma)))) / beta
    #     term5 = (2 * h * sigma * np.log(1 +
    #              np.exp(-beta * (mu + h * sigma)))) / beta
    #     term6 = ((Lambda * (Lambda + 2 * h * beta * sigma)) / (1 + np.exp(Lambda - beta * mu + h *
    #              beta * sigma)) + 2 * h * beta * sigma * np.log(1 + np.exp(-Lambda + beta * mu - h * beta * sigma))) / beta**2
    #     term7 = (Lambda * (Lambda + 2 * h * beta * sigma) + 2 * (1 + np.exp(Lambda + beta * mu + h * beta * sigma)) * h * beta *
    #              sigma * np.log(1 + np.exp(-Lambda - beta * (mu + h * sigma)))) / ((1 + np.exp(Lambda + beta * (mu + h * sigma))) * beta**2)
    #
    #     result = term1 * (term2 + term3 + term4 + term5 - term6 + term7)
    #     return result

    # def analythic_mean_field_solution(sigma: float, Lambda: float, h: float, sigma_0: float, mu: float, T: float):
    #     beta=1 / T
    #     term1=h**2 * sigma
    #     term2=h * beta * (sigma - sigma_0)
    #     term3=np.log(1 + np.exp(beta * (mu - h * sigma)))
    #     term4=np.log(1 + np.exp(-beta * (mu + h * sigma)))
    #
    #     result=(term1 * (term2 + term3 + term4)) / (np.pi * beta)
    #     return result

    def analythic_mean_field_solution(sigma: float, Lambda: float, h: float,
                                      sigma_0: float, mu: float, T: float):
        beta = 1/T
        e = np.sqrt(Lambda**2 + (h*sigma)**2)
        e_vac = np.sqrt(Lambda**2 + (h*sigma_0)**2)

        exp_plus = beta * (h*sigma + mu)
        exp_minus = beta * (h*sigma - mu)

        vacuum = - (h*sigma*(Lambda**2 + 2*h*sigma * (h*sigma - e)))/e \
            + (h*sigma*(Lambda**2 + 2*h*sigma_0 * (h*sigma_0 - e_vac)))/e_vac

        thermal = 2*h*sigma*log_expit(-exp_plus) + \
            2*h*sigma*log_expit(-exp_minus)
        thermal /= beta

        regularization_1 = Lambda * \
            (Lambda + 2 * h * beta * sigma) * expit(Lambda + exp_minus) \
            + 2*h * beta*sigma * log_expit(-Lambda-exp_minus)
        regularization_2 = Lambda * \
            (Lambda + 2 * h * beta * sigma) * expit(Lambda + exp_plus) \
            + 2*h * beta*sigma * log_expit(-Lambda-exp_plus)
        regularization = (regularization_1 + regularization_2) / beta**2

        return h/(2*np.pi) * (vacuum + thermal - regularization)

    sigma_max = 6.0
    Lambda = 1e5
    kir = 1e-4
    N_Grid = 1000
    samples = 100
    h = 1
    sigma_0 = 1

    # test vacuum solution
    def compute_2_1_sol(mu: float, T: float, N_Flavor: int) -> Tuple[float, float]:
        sol = GN_2p1(sigma_max, Lambda, kir, N_Grid,
                     samples, mu, T, N_Flavor=N_Flavor)
        y = sol.return_data.solution
        x = sol.return_data.grid
        time = sol.return_data.time
        return x, y, time

    fig, ((ax1_1, ax1_2), (ax2_1, ax2_2)) = plt.subplots(2, 2)

    mu = 0.0
    T = 0.0
    N_Flavor = np.infty
    x, y, time = compute_2_1_sol(mu, T, N_Flavor)

    y_ref = [analytic_vacuum_solution(x_i, Lambda, h, sigma_0) for x_i in x]
    ax1_1.plot(x, y[-1], label="final")
    ax1_1.plot(x, y_ref, label="reference", linestyle="dashed")
    ax1_1.set_xlim([0, 1.2])
    ax1_1.set_ylim([-.1, .1])
    ax1_1.legend()
    ax1_1.grid()

    ax1_2.plot(x, y[-1], label="final")
    ax1_2.plot(x, y_ref, label="reference", linestyle="dashed")
    ax1_2.legend()
    ax1_2.grid()

    mu = 0.1
    T = 0.4
    N_Flavor = np.infty
    x, y, time = compute_2_1_sol(mu, T, N_Flavor)

    y_ref = [analythic_mean_field_solution(
        x_i, Lambda, h, sigma_0, mu, T) for x_i in x]
    ax2_1.plot(x, y[-1], label="final")
    ax2_1.plot(x, y_ref, label="reference", linestyle="dashed")
    ax2_1.set_xlim([0, 1.2])
    ax2_1.set_ylim([-.1, .1])
    ax2_1.legend()
    ax2_1.grid()

    ax2_2.plot(x, y[-1], label="final")
    ax2_2.plot(x, y_ref, label="reference", linestyle="dashed")
    ax2_2.legend()
    ax2_2.grid()
    plt.show()


if __name__ == "__main__":
    run_1_1_tests()
    run_2_1_tests()
