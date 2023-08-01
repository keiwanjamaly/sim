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
from python_files.gross_neveu.couplings.couplings_io import get_all_pre_computed_alphas_and_lambdas
from python_files.gross_neveu.compute_observable import sigma as calculate_sigma
from python_files.gross_neveu.Gross_Neveu import get_model
from lmfit import Model


def calculate_sigma_difference(one_over_g2, dimension, sigma_max, Lambda, kir, delta_sigma, N_Flavor, h, sigma_0):
    mu = 0.0
    T = 0.01
    sigma_0_ir = calculate_sigma(one_over_g2, dimension, mu, T,
                                 sigma_max, Lambda, kir, delta_sigma, N_Flavor, h, sigma_0)
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


def compute_couping_fit(Lambda, save=None, save_fit=None, plot=True):
    flavours, couplings = get_computed_couplings_from_file(Lambda, './data')

    coupling_mf = get_model(dimension=2).calculate_one_g2(1.0, 1.0, Lambda)

    popt, pcov = curve_fit(lambda x, a: coupling_fit_function(
        x, a, coupling_mf), flavours, couplings)

    print(Lambda, popt[0], coupling_mf)
    gmodel = Model(lambda x, a: coupling_fit_function(
        x, a, coupling_mf))
    result = gmodel.fit(couplings, x=flavours, a=200)
    alpha_val = result.best_values['a']
    alpha_err = np.sqrt(result.covar)[0][0]

    skip_plotting = False

    if save is not None:
        for flavour, coupling in zip(flavours, couplings):
            print(f'{flavour}\t{coupling}', file=save)
        skip_plotting = True

    xdata = np.linspace(2, 16)
    fdata = coupling_fit_function(xdata, alpha_val, coupling_mf)
    if save_fit is not None:
        for x, y in zip(xdata, fdata):
            print(f'{x}\t{y}\t{coupling_mf}', file=save_fit)
        skip_plotting = True

    if not skip_plotting and plot:
        plt.axhline(y=coupling_mf, color='r', linestyle='-',
                    label="mean field couplings")
        plt.plot(xdata, fdata)
        plt.title(f'$\\Lambda = {Lambda:.0f}$')
        plt.scatter(flavours, couplings,
                    label=f'$2/(\\pi * g^2_{{mf}}) \\cdot \\arctan({popt[0]:.1f} \\cdot N_f)$')
        plt.xlabel(r'$N_f$')
        plt.ylabel(r'$\frac{1}{g^2}$', rotation=0)

        plt.legend()
        plt.show()
    return alpha_val, alpha_err, coupling_mf


def compute_all_coupling_fits():
    Lambda_Array = get_all_cutoffs_with_more_than_two_couplings("./data")
    data = []
    for Lambda in Lambda_Array:
        beta, beta_error, mean_field_coupling = compute_couping_fit(
            Lambda, plot=False)
        data.append([Lambda, beta, mean_field_coupling, beta_error])
    data = np.array(data)
    with h5py.File('./data/couping_pre_computations.hdf5', "w") as f:
        f.create_dataset("couplings", data=data)


def plot_alpha_vs_lambda(dir: str, save=None):
    # def f(x, a, b, c, d, e):
    #     return b * x**a + c * x ** d
    def f(x, a, b, d):
        return b * x**a + x/4 + d
    # def f(x, a, b):
    #     return b * x**a + x/4
    Lambdas, alphas, alpha_errors = get_all_pre_computed_alphas_and_lambdas(
        dir)
    print(Lambdas, alphas, alpha_errors)
    if save is not None:
        for Lambda, alpha, alpha_error in zip(Lambdas, alphas, alpha_errors):
            print(Lambda, alpha, alpha_error, file=save)
    gmodel = Model(f)
    result = gmodel.fit(alphas, x=Lambdas,
                        weights=alpha_errors, a=0.26, b=0.66, d=11.0)
    print(result.fit_report())
    fig = result.plot()
    axes = fig.gca()
    axes.set_xscale('log')
    plt.show()
