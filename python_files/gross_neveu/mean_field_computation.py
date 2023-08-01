import numpy as np


def thermal_part(sigma, mu, T, h, sigma0):
    return T * np.log(1+np.exp(-(mu + h * sigma)/T))


def compute_potential_2p1(sigma, mu, T, h=1.0, sigma0=1.0):
    result = h * (sigma - sigma0)
    result += thermal_part(sigma, mu, T, h, sigma0)
    result += thermal_part(sigma, -mu, T, h, sigma0)
    result *= h**2 * sigma / np.pi
    return result
