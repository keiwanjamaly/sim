import numpy as np
from numba import njit, objmode
import timeit
import sys


@njit(cache=True)
def sech(x):
    return 1/np.cosh(x)


@njit(cache=True)
def csch(x):
    return 1/np.sinh(x)


@njit(cache=True)
def e_f(k, sigma):
    return np.sqrt(k**2 + sigma**2)


@njit(cache=True)
def e_b(k, u_x):
    return np.sqrt(k**2 + u_x)


@njit(cache=True)
def n_f(x):
    return 1/(np.exp(x) + 1)


@njit(cache=True)
def n_b(x):
    return 1/(np.exp(x) - 1)


@njit(cache=True)
def S(k, sigma, *args):
    spatial_dimension = args[12]
    e = e_f(k, sigma)
    mu = args[9]
    beta = args[10]
    minus = beta * (e-mu) / 2
    plus = beta * (e+mu) / 2

    match spatial_dimension:
        case 1:
            # for (1+1)
            return sigma*k**3/(4*e**3*np.pi) * (e * beta * (sech(minus)**2 + sech(plus)**2)
                                                - 2*(np.tanh(minus) + np.tanh(plus)))
        case 2:
            # for (2+1)
            return sigma*k**4/(8*e**3*np.pi) * (e * beta * (sech(minus)**2 + sech(plus)**2)
                                                - 2*(np.tanh(minus) + np.tanh(plus)))


@njit(cache=True)
def Q(k, ux, *args):
    spatial_dimension = args[12]
    beta = args[10]
    N_Flavor = args[11]

    e = e_b(k, ux)

    match spatial_dimension:
        case 1:
            # for (1+1)
            return - k**3 / (2*np.pi*e*N_Flavor) * (1+2*n_b(beta * e))
        case 2:
            # for (2+1)
            return - k**4 / (4*np.pi*2*e*N_Flavor) * (1+2*n_b(beta * e))


@njit(cache=True)
def compute_diffusion(k, u, *args):
    mean_field_flag = args[1]
    grid = args[3]
    dx_direct = args[4]
    dx_half = args[5]
    c_0 = args[6]
    c_1 = args[7]
    c_2 = args[8]
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
def f(t, u, *args):
    """
    The *args conventions are as follows
    *args = (Lambda, mean_field_flag, console_logging_flag, grid, dx_direct,
      dx_half, c_0, c_1, c_2, mu, beta, N_Flavor, spatial_dimension, time_start, kir, iterations)
    """

    Lambda = args[0]
    k = Lambda * np.exp(-t)
    grid = args[3]

    diffusion = compute_diffusion(k, u, *args)

    source = S(k, grid, *args)

    console_logging_flag = args[2]
    time_start = args[13]
    tir = args[14]

    if console_logging_flag:
        args[-1][0] += 1
        if args[-1][0] > 0:
            with objmode():
                time_elapsed = timeit.default_timer() - time_start
                print_string = '{time:.7f} ({kval:.5e})/{tirval:.2f}; time elapsed = {te:.2f} seconds'.format(
                    time=t, kval=k, tirval=tir, te=time_elapsed)
                args[-1][0] -= 1000
                print(print_string, end="\r")

    return diffusion + source
