from python_files.gross_neveu.compute_observable import sigma


def compute_simga(one_over_g2, dimension, mu, T, sigma_max,
                  Lambda, kir, delta_sigma, N_Flavor, h, sigma_0):
    try:
        result = sigma(one_over_g2, dimension, mu, T, sigma_max,
                       Lambda, kir, delta_sigma, N_Flavor, h, sigma_0)
        print(f'T = {T}, mu = {mu} - done')
        return result
    except:
        print(f'T = {T}, mu = {mu} - error')
        return -1
