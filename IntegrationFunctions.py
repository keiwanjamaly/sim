import numpy as np
from numba import njit


@njit
def sech(x):
    return 1/np.cosh(x)


@njit
def csch(x):
    return 1/np.sinh(x)


@njit
def e_f(k, sigma):
    return np.sqrt(k**2 + sigma**2)


@njit
def e_b(k, u_x):
    return np.sqrt(k**2 + u_x)


@njit
def n_f(x):
    return 1/(np.exp(x) + 1)


@njit
def n_b(x):
    return 1/(np.exp(x) - 1)


@njit
def S(k, sigma, mu, beta, spatial_dimension):
    e = e_f(k, sigma)
    mu = mu
    beta = beta
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


@njit
def Q(k, ux, beta, N_Flavor, spatial_dimension):
    beta = beta
    N = N_Flavor
    e = e_b(k, ux)

    match spatial_dimension:
        case 1:
            # for (1+1)
            return - k**3 / (2*np.pi*e*N) * (1+2*n_b(beta * e))
        case 2:
            # for (2+1)
            return - k**4 / (4*np.pi*2*e*N) * (1+2*n_b(beta * e))


@njit
def compute_diffusion(k, u, mean_field_flag, grid, dx_direct, dx_half, c_0, c_1, c_2, beta, N_Flavor, spatial_dimension):
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
        Q_cal = Q(k, ux, beta, N_Flavor, spatial_dimension)

        return (Q_cal[1:] - Q_cal[:-1])/dx_half


@njit
def f(t, u, Lambda, mean_field_flag, console_logging_flag, grid, dx_direct,
      dx_half, c_0, c_1, c_2, mu, beta, N_Flavor, spatial_dimension, time_start, kir):
    k = Lambda * np.exp(-t)

    diffusion = compute_diffusion(k, u, mean_field_flag, grid, dx_direct,
                                  dx_half, c_0, c_1, c_2, beta, N_Flavor, spatial_dimension)

    # if console_logging_flag:
    #     time_elapsed = timeit.default_timer() - time_start
    #     print_string = '{time:.7f} ({k:.5e})/{kirval:.2f}; time elapsed = {te:.2f} seconds'.format(
    #         time=t, kval=k, kirval=kir, te=time_elapsed)
    #     print(print_string)

    return diffusion + S(k, grid, mu, beta, spatial_dimension)
