import numpy as np
from scipy.optimize import newton
import h5py
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import timeit
from python_files.gross_neveu.couplings.couplings_io import get_closest_coupling_approximation_or_MF_coupling
from python_files.gross_neveu.couplings.couplings_io import save_coupling
from python_files.gross_neveu.couplings.couplings_io import get_computed_couplings_from_file
from python_files.gross_neveu.couplings.couplings_io import coupling_fit_function
from python_files.gross_neveu.couplings.couplings_io import get_all_cutoffs_with_more_than_two_couplings
from python_files.gross_neveu.compute_observable import sigma as calculate_sigma
from python_files.gross_neveu.Gross_Neveu import get_model

def calculate_sigma_difference(one_over_g2, dimension, sigma_max, Lambda, kir, delta_sigma, N_Flavor, h, sigma_0):
    sigma_0_ir = calculate_sigma(one_over_g2 ,dimension, sigma_max, Lambda, kir, delta_sigma, N_Flavor, h, sigma_0)
    sigma_0 = 1.0
    result = sigma_0_ir - sigma_0
    print(f'with 1/g^2 = {one_over_g2} the deviation of sigma_0 = {result}')
    return result

def compute_couping(Lambda, N_Flavor):
    kir = 1e-2
    h = 1
    sigma_0 = 1.0
    sigma_max = 2000.0
    delta_sigma = 0.006

    dimension = 2
    dir = "./data"

    # use mean field coupling or an approximation of the fit for the initial parameter
    one_over_g2 = get_closest_coupling_approximation_or_MF_coupling(
        Lambda, N_Flavor, dir, dimension)

    start = timeit.default_timer()
    result = newton(calculate_sigma_difference, one_over_g2, args=(
        dimension, sigma_max, Lambda, kir, delta_sigma, N_Flavor, h, sigma_0))
    time = timeit.default_timer() - start

    save_coupling(Lambda, result, N_Flavor, time, dir)


def compute_couping_fit(Lambda, plot=False):
    flavours, couplings = get_computed_couplings_from_file(Lambda, './data')

    coupling_mf = get_model(dimension=2).calculate_one_g2(1.0, 1.0, Lambda)

    popt, pcov = curve_fit(lambda x, a: coupling_fit_function(x, a, coupling_mf), flavours, couplings)

    if plot:
        plt.axhline(y=coupling_mf, color='r', linestyle='-',
                    label="mean field couplings")
        xdata = np.linspace(2, 16)
        plt.plot(xdata, coupling_fit_function(xdata, popt[0], coupling_mf))
        plt.title(f'$\\Lambda = {Lambda:.0f}$')
        plt.scatter(flavours, couplings,
                    label=f'$2/(\\pi * g^2_{{mf}}) \\cdot \\arctan({popt[0]:.1f} \\cdot N_f)$')
        plt.xlabel(r'$N_f$')
        plt.ylabel(r'$\frac{1}{g^2}$', rotation=0)

        plt.legend()
        plt.show()
    return popt[0], coupling_mf


def compute_all_coupling_fits():
    Lambda_Array = get_all_cutoffs_with_more_than_two_couplings("./data")
    data = []
    for Lambda in Lambda_Array:
        beta, mean_field_coupling = compute_couping_fit(Lambda)
        data.append([Lambda, beta, mean_field_coupling])
    data = np.array(data)
    with h5py.File('./data/couping_pre_computations.hdf5', "w") as f:
        f.create_dataset("couplings", data=data)


