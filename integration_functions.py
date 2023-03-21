import numpy as np
from numba import njit


def generate_filename(mu, T, sigma_max, n_grid, kir, tolerances, *args):
    Lambda = get_lambda(*args)
    n_Flavor = get_n_flavor(*args)
    spatial_dimension = get_spatial_dimension(*args)
    return f'mu={mu}_T={T}_sigmaMax={sigma_max}_Lambda={Lambda}_kir={kir}_nGrid={n_grid}_nFlavor={n_Flavor}_tolerance={tolerances:e}_d={spatial_dimension}.hdf5'


def generate_args(spatial_dimension, Lambda, mu, T, n_flavor):
    return (spatial_dimension, Lambda, mu, 1/T, n_flavor)


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


@njit
def get_lambda(*args):
    return args[1]


@njit
def get_spatial_dimension(*args):
    return args[0]


@njit
def get_spatial_mu(*args):
    return args[2]


@njit
def get_beta(*args):
    return args[3]


@njit
def get_n_flavor(*args):
    return args[4]


@njit(cache=True)
def initial_condition(grid, *args):
    print(args)
    spatial_dimension = get_spatial_dimension(*args)
    Lambda = get_lambda(*args)

    match spatial_dimension:
        case 1:
            # for (1+1)
            intermediate = 1/np.sqrt(1+(1/Lambda)**2)
            return (grid/np.pi)*(np.arctanh(intermediate) - intermediate)
        case 2:
            # for (2+1)
            intermediate = (2+Lambda**2 - 2*np.sqrt(1+Lambda**2)
                            ) / (2*np.pi*np.sqrt(1+Lambda**2))
            return intermediate*grid


@njit(cache=True)
def S(k, sigma, *args):
    spatial_dimension = get_spatial_dimension(*args)
    e = e_f(k, sigma)
    mu = get_spatial_mu(*args)
    beta = get_beta(*args)
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
    spatial_dimension = get_spatial_dimension(*args)
    beta = get_beta(*args)
    N_Flavor = get_n_flavor(*args)

    e = e_b(k, ux)

    match spatial_dimension:
        case 1:
            # for (1+1)
            return - k**3 / (2*np.pi*e*N_Flavor) * (1+2*n_b(beta * e))
        case 2:
            # for (2+1)
            return - k**4 / (4*np.pi*2*e*N_Flavor) * (1+2*n_b(beta * e))
