import numpy as np
import os.path
import argparse
from scipy.optimize import newton
import Gross_Neveu
from compute_observables import DataClass
from grid_creator import create_inhomogenious_grid_from_cell_spacing
import h5py
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import timeit


def calculate_sigma(one_over_g2: float, Model, sigma_max, Lambda, kir,
                    delta_sigma, N_Flavor, h, sigma_0):
    mu = 0.0
    T = 0.01
    samples = 3
    N_Grid = 1000
    grid_points = create_inhomogenious_grid_from_cell_spacing(
        sigma_max, delta_sigma)
    model = Model(grid_points, Lambda, kir, samples,
                  mu, T, N_Flavor, h, one_over_g2, sigma_0)

    y = model.return_data.solution
    x = model.return_data.grid
    time = model.return_data.time

    data_class = DataClass(x, time, y)
    sigma_0_ir = data_class.sigma_0[-1]

    result = sigma_0_ir - sigma_0
    # plt.plot(x, y[-1])
    # plt.show()
    print(f'with 1/g^2 = {one_over_g2} the deviation of sigma_0 = {result}')
    return result


def get_coupling_from_file(Lambda: float, N_Flavor: float, model, sigma_0: float = 1.0, h: float = 1.0):
    filename = f'./data/couplings_Lambda={Lambda}_N={N_Flavor}.hdf5'
    if np.isinf(N_Flavor):
        return model.calculate_one_g2(h, sigma_0, Lambda)
    elif os.path.isfile(filename):
        with h5py.File(filename, "r") as f:
            return f.attrs["coupling"]
    else:
        with h5py.File('./data/couping_pre_computations.hdf5', "r") as f:
            if Lambda in f["couplings"][:, 0]:
                pos = list(f["couplings"][:, 0]).index(Lambda)
                beta, alpha = f["couplings"][pos][1:]
                return alpha * np.arctan(beta * N_Flavor) * 2 / np.pi
            else:
                print(
                    "There is no reference data for this Cutoff, using mean field coupling instead")
                return model.calculate_one_g2(h, sigma_0, Lambda)


def compute_couping(Lambda, N_Flavor):
    kir = 1e-2
    h = 1
    sigma_0 = 1.0
    sigma_max = 2000.0
    delta_sigma = 0.006

    model = Gross_Neveu.GN_2p1

    # use mean field coupling or an approximation of the fit for the initial parameter
    one_over_g2 = get_coupling_from_file(Lambda, N_Flavor, model)

    start = timeit.default_timer()
    result = newton(calculate_sigma, one_over_g2, args=(
        model, sigma_max, Lambda, kir, delta_sigma, N_Flavor, h, sigma_0))
    time = timeit.default_timer() - start
    print(f'N = {N_Flavor}, 1/g^2 = {result} it took {time:.2f} seconds')

    with h5py.File(f'./data/couplings_Lambda={Lambda}_N={N_Flavor}.hdf5', "w") as f:
        f.attrs["sigma_max"] = sigma_max
        f.attrs["Lambda"] = Lambda
        f.attrs["kir"] = kir
        f.attrs["coupling"] = result
        f.attrs["N_Flavor"] = N_Flavor
        f.attrs["time"] = time


def compute_couping_fit(Lambda, plot=False):
    couplings = []
    N_Flavor_List = list(range(2, 17))
    for N_Flavor in N_Flavor_List:
        with h5py.File(f'./data/couplings_Lambda={Lambda}_N={N_Flavor}.hdf5', "r") as f:
            couplings.append(f.attrs["coupling"])

    coupling_mf = Gross_Neveu.GN_2p1.calculate_one_g2(1.0, 1.0, Lambda)

    def func(x, a):
        return coupling_mf * np.arctan(a*x) * 2 / np.pi

    popt, pcov = curve_fit(func, N_Flavor_List, couplings)

    if plot:
        plt.axhline(y=coupling_mf, color='r', linestyle='-',
                    label="mean field couplings")
        xdata = np.linspace(2, 16)
        plt.plot(xdata, func(xdata, *popt))
        plt.title(f'$\\Lambda = {Lambda:.0f}$')
        plt.scatter(N_Flavor_List, couplings,
                    label=f'$2/(\\pi * g^2_{{mf}}) \\cdot \\arctan({popt[0]:.1f} \\cdot N_f)$')
        plt.xlabel(r'$N_f$')
        plt.ylabel(r'$\frac{1}{g^2}$', rotation=0)

        plt.legend()
        plt.show()
    return popt[0], coupling_mf


def compute_all_coupling_fits(Lambda):
    Lambda_Array = [10.0, 200.0, 500.0, 1000.0]
    data = []
    for Lambda in Lambda_Array:
        beta, mean_field_coupling = compute_couping_fit(Lambda)
        data.append([Lambda, beta, mean_field_coupling])
    print(data)
    data = np.array(data)
    with h5py.File('./data/couping_pre_computations.hdf5', "w") as f:
        f.create_dataset("couplings", data=data)


def compute_time_fits(Lambda):
    def func(x, a, b):
        return a * np.exp(b * x)

    N_Flavor_List = list(range(2, 17))

    time = []

    for N_Flavor in N_Flavor_List:
        with h5py.File(f'./data/couplings_Lambda={Lambda}_N={N_Flavor}.hdf5', "r") as f:
            time.append(f.attrs["time"])

    time_log = np.log(time)
    popt, pcov = curve_fit(lambda x, a, b: np.log(
        func(x, a, b)), N_Flavor_List, time_log)

    print(popt)
    xdata = np.linspace(2, 16)
    plt.plot(xdata, func(xdata, *popt))
    plt.scatter(N_Flavor_List, time)
    plt.title(f'exponential is {popt[0]:.1f}exp({popt[1]:.2f}*x)')
    plt.yscale('log')

    plt.show()


def main():
    parser = argparse.ArgumentParser(prog='couplings calculator')
    parser.add_argument(
        '-N', type=int, help='Set the number of flavors', required=True, default=None)
    parser.add_argument('-L', type=float,
                        help='Set the UV Cutoff', required=True, default=None)
    parser.add_argument('--time',
                        help='extrapolate computation time data', action='store_true')
    parser.add_argument('--fit',
                        help='extraplate coupling data', action='store_true')
    parser.add_argument('--all',
                        help='extraplate coupling data', action='store_true')

    args = parser.parse_args()

    if args.time:
        compute_time_fits(args.L)
    elif args.fit:
        if args.all:
            compute_all_coupling_fits(args.L)
        else:
            compute_couping_fit(args.L, True)
    else:
        compute_couping(args.L, args.N)


if __name__ == "__main__":
    # compute_time_fits()
    main()
    # data = []
    # for Lambda in [10.0, 200.0, 500.0, 1000.0]:
    #     beta = compute_time_fits(Lambda)[0]
    #     data.append([Lambda, beta])
    #     print(data)
    # with h5py.File('./data/couping_pre_computations.hdf5', "w") as f:
