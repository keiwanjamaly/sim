from python_files.gross_neveu.compute_observable import sigma


def compute_sigma_spread(x):
    return compute_sigma(*x)


def compute_sigma(one_over_g2, dimension, mu, T, sigma_max,
                  Lambda, kir, delta_sigma, N_Flavor, h, sigma_0):
    try:
        result = sigma(one_over_g2, dimension, mu, T, sigma_max,
                       Lambda, kir, delta_sigma, N_Flavor, h, sigma_0)
        print(f'mu = {mu}, T = {T}, result = {result} - done')
    except:
        result = -1
        print(f'mu = {mu}, T = {T} result = {result}- error')

    return mu, T, result
