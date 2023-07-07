import numpy as np


def analytic_2p1_MF(mu: float, T: float, h: float = 1, sigma_0: float = 1) -> float:
    intermediate = np.exp(h*sigma_0/T)/2 - 1
    if intermediate >= 1:
        mu_crit = T*np.arccosh(intermediate)
    else:
        mu_crit = -np.infty

    if mu <= mu_crit:
        return T/(h*sigma_0) * np.arccosh(np.exp(h*sigma_0/T)/2 - np.cosh(mu/T))
    else:
        return 0
