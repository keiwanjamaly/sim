from python_files.gross_neveu.Gross_Neveu import GN_1p1, GN_2p1
from python_files.grid_creator import create_homogenious_grid
import numpy as np
import scipy
import json


def test_1_plus_1_MF():
    sigma_max = 6.0
    Lambda = 1e5
    kir = 1e-4
    N_Grid = 1000
    samples = 100
    mu = 0.1
    T = 0.1
    N_Flavor = np.infty
    grid_points = create_homogenious_grid(sigma_max, N_Grid)
    sol = GN_1p1(grid_points, Lambda, kir,
                 samples, mu, T, N_Flavor=N_Flavor)
    y = sol.return_data.solution[-1]
    x = sol.return_data.grid
    time = sol.return_data.time
    with open('./1+1_test_files/flow_MF_T=0.1,mu=0.1.json') as user_file:
        parsed_json = json.load(user_file)
        x_ref = np.array(parsed_json["graphs"][-1]["x"])
        y_ref = np.array(parsed_json["graphs"][-1]["y"])
    assert (x_ref == x).all()
    np.testing.assert_allclose(y, y_ref, atol=1e-5, rtol=0)


def test_1_plus_1_N_2():
    sigma_max = 6.0
    Lambda = 1e5
    kir = 1e-4
    N_Grid = 1000
    samples = 100
    mu = 0.1
    T = 0.1
    N_Flavor = 2
    grid_points = create_homogenious_grid(sigma_max, N_Grid)
    sol = GN_1p1(grid_points, Lambda, kir,
                 samples, mu, T, N_Flavor=N_Flavor)
    y = sol.return_data.solution[-1]
    x = sol.return_data.grid
    time = sol.return_data.time
    with open('1+1_test_files/flow_N=2,T=0.1,mu=0.1.json') as user_file:
        parsed_json = json.load(user_file)
        x_ref = np.array(parsed_json["graphs"][-1]["x"])
        y_ref = np.array(parsed_json["graphs"][-1]["y"])
    assert (x_ref == x).all()
    np.testing.assert_allclose(y, y_ref, atol=1e-3, rtol=0)


def analytic_vacuum_solution(sigma: float, Lambda: float, h: float, sigma_0: float):
    tmp = np.sqrt(Lambda**2 + (h*sigma_0)**2)
    d_gamma = GN_2p1.dimension_gamma()
    e_vac = np.sqrt(Lambda**2 + h**2 * sigma_0**2)
    e = np.sqrt(Lambda**2 + h**2 * sigma**2)
    tmp = -6*e + 6*e_vac + 6*h * \
        (sigma - sigma_0) + 3*Lambda**2*(1/e - 1/e_vac)
    return (d_gamma*h**2*sigma*tmp)/(24*np.pi)


def analythic_mean_field_solution(sigma: float, Lambda: float, h: float,
                                  sigma_0: float, mu: float, T: float):
    def expit(x):
        return scipy.special.expit(-x)

    def log_expit(x):
        return scipy.special.log1p(np.exp(x))
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


def test_2_plus_1_vacuum():
    sigma_max = 6.0
    Lambda = 1e5
    kir = 1e-4
    N_Grid = 1000
    samples = 100
    mu = 0.0
    T = 0.0
    h = 1
    sigma_0 = 1
    N_Flavor = np.infty
    grid_points = create_homogenious_grid(sigma_max, N_Grid)
    sol = GN_2p1(grid_points, Lambda, kir,
                 samples, mu, T, N_Flavor=N_Flavor)
    y = sol.return_data.solution[-1]
    x = sol.return_data.grid
    time = sol.return_data.time
    y_ref = [analytic_vacuum_solution(x_i, Lambda, h, sigma_0) for x_i in x]
    np.testing.assert_allclose(y, y_ref, atol=1e-4, rtol=0)


def test_2_plus_1_N_2():
    sigma_max = 6.0
    Lambda = 1e5
    kir = 1e-4
    N_Grid = 1000
    samples = 100
    mu = 0.1
    T = 0.4
    h = 1
    sigma_0 = 1
    N_Flavor = np.infty
    grid_points = create_homogenious_grid(sigma_max, N_Grid)
    sol = GN_2p1(grid_points, Lambda, kir,
                 samples, mu, T, N_Flavor=N_Flavor)
    y = sol.return_data.solution[-1]
    x = sol.return_data.grid
    time = sol.return_data.time

    y_ref = [analythic_mean_field_solution(
        x_i, Lambda, h, sigma_0, mu, T) for x_i in x]
    np.testing.assert_allclose(y, y_ref, atol=1e-4, rtol=0)
