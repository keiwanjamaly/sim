from python_files.gross_neveu.compute_observable import sigma
# from multiprocessing import Lock


def compute_sigma_spread(x):
    return compute_sigma(*x)


def compute_sigma(one_over_g2, dimension, mu, T, sigma_max,
                  Lambda, kir, delta_sigma, N_Flavor, h, sigma_0):
    result = sigma(one_over_g2, dimension, mu, T, sigma_max,
                   Lambda, kir, delta_sigma, N_Flavor, h, sigma_0)

    print(f'mu = {mu}, T = {T}, result = {result} - done')

    return mu, T, result
