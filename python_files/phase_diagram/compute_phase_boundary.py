import numpy as np
from python_files.gross_neveu.Gross_Neveu import get_model
from python_files.gross_neveu.couplings.couplings_io import get_exact_coupling_from_file
from python_files.gross_neveu.couplings.couplings_io import generate_filename
from python_files.phase_diagram.computation_function import compute_u, compute_sigma
import h5py
from joblib import Parallel, delayed
from python_files.data_class import Potential, Line
import matplotlib.pyplot as plt
import io


def get_curve_point(t, angle, r_min, r_max, origin):
    r = t * r_max + (1-t) * r_min
    x = r * np.cos(angle)
    y = r * np.sin(angle)
    x += origin[0]
    y += origin[1]
    return x, y


def compute_point(Lambda, N_Flavor, mu, T):
    mu_max = 1.2
    T_max = 1.0
    if np.isinf(N_Flavor):
        one_over_g2 = get_model(2).calculate_one_g2(
            h=1.0, sigma_0=1.0, Lambda=Lambda)
    else:
        filename = generate_filename(Lambda, N_Flavor, "./data")
        one_over_g2 = get_exact_coupling_from_file(filename)
    dimension = 2
    sigma_max = 2000
    kir = 1e-2
    delta_sigma = 0.006
    h = 1.0
    sigma_0 = 1.0
    mu, T, x, result = compute_u(one_over_g2, dimension, mu, T, sigma_max,
                                 Lambda, kir, delta_sigma, N_Flavor, h, sigma_0)
    potential = Potential(x, result)
    sigma = potential.sigma
    return sigma, x, result


def compute_ray(Lambda, N_Flavor, angle, r_min, r_max, origin):
    a = 0
    mu, T = get_curve_point(a, angle, r_min, r_max, origin)
    f_a, x, result_a = compute_point(Lambda, N_Flavor, mu, T)
    b = 1
    mu, T = get_curve_point(b, angle, r_min, r_max, origin)
    f_b, _, result_b = compute_point(Lambda, N_Flavor, mu, T)
    while np.abs(a - b) > 0.0001:
        new = (a + b)/2
        mu, T = get_curve_point(new, angle, r_min, r_max, origin)
        f_new, x, result_new = compute_point(Lambda, N_Flavor, mu, T)
        if f_new <= 0:
            b, f_b, result_b = new, f_new, result_new
        else:
            a, f_a, result_a = new, f_new, result_new

    return get_curve_point(b, angle, r_min, r_max, origin), x, result_b


def compute_boundary(Lambda, N_Flavor, save: io.TextIOWrapper, save_LP: io.TextIOWrapper):
    angles = np.linspace(0, np.pi/2, 1000)
    # results = []
    job_list = []
    for i, angle in enumerate(angles):
        # results.append(compute_ray(Lambda, N_Flavor,
        #                angle, 0.5, 1.2, (0,  0.01)))
        job_list.append([Lambda, N_Flavor, angle, 0.5, 1.2, (0, 0.01)])
        # print(f'done - {i}/{len(angles)-1}')

    results = Parallel(n_jobs=-1, verbose=20, batch_size=1)(
        delayed(compute_ray)(*x) for x in job_list)
    mus = []
    Ts = []
    for result in results:
        mus.append(result[0][0])
        Ts.append(result[0][1])

    forth_divs = [
        Potential(result[1], result[2]).forth_div_at_zero for result in results]
    # compute sign change
    sign_change_interval = None
    for i in range(len(forth_divs) - 1):
        if (forth_divs[i] > 0 and forth_divs[i+1] < 0) or (forth_divs[i] < 0 and forth_divs[i+1] > 0):
            sign_change_interval = (i, i+1)
            print('sign change:', sign_change_interval)

    plt.plot(mus, Ts, 'bx')

    if sign_change_interval is not None:
        previous_forth_div = forth_divs[sign_change_interval[0]]
        next_forth_div = forth_divs[sign_change_interval[1]]
        previous_mu, previous_T = results[sign_change_interval[0]][0]
        next_mu, next_T = results[sign_change_interval[1]][0]
        print(previous_forth_div, next_forth_div)
        print(previous_mu, previous_T)
        print(next_mu, next_T)
        root = Line.from_points(
            next_forth_div, 0, previous_forth_div, 1).get_root()
        mu_LP = next_mu * (1-root) + previous_mu * root
        T_LP = next_T * (1-root) + previous_T * root
        print(mu_LP, T_LP)

        plt.plot([mu_LP], [T_LP], 'ro')
        if save_LP is not None:
            print(f'{mu_LP}\t{T_LP}', file=save_LP)

    if save is None:
        plt.show()
    else:
        for mu, T in zip(mus, Ts):
            print(f'{mu}\t{T}', file=save)

    # print(result[0][0], result[0][1])
