import numpy as np
import argparse
from scipy.optimize import newton
import Gross_Neveu
from compute_observables import DataClass
from grid_creator import create_inhomogenious_grid_from_cell_spacing
import h5py


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


def compute_couping(Lambda, N_Flavor):
    kir = 1e-2
    h = 1
    sigma_0 = 1.0
    sigma_max = 2000.0
    delta_sigma = 0.006

    model = Gross_Neveu.GN_2p1
    one_over_g2 = model.calculate_one_g2(h, sigma_0, Lambda)
    result = newton(calculate_sigma, one_over_g2, args=(
        model, sigma_max, Lambda, kir, delta_sigma, N_Flavor, h, sigma_0))
    print(f'N = {N_Flavor}, 1/g^2 = {result}')

    with h5py.File(f'./data/couplings_Lambda={Lambda}_N={N_Flavor}.hdf5', "w") as f:
        f.attrs["sigma_max"] = sigma_max
        f.attrs["Lambda"] = Lambda
        f.attrs["kir"] = kir
        f.attrs["coupling"] = result
        f.attrs["N_Flavor"] = N_Flavor


def compute_time_fits(Lambda):
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit

    def func(x, a, b):
        return a * np.exp(b * x)

    runs = {
        10: ([2, 3, 4, 5, 6, 7], [37.5, 58.48, 152.44, 342.30, 1234.52, 1919.02]),
        100: ([2, 3, 4, 5, 6, 7], [34.08, 54.0, 89.87, 468.4, 824.78, 1439.66]),
        200: ([2, 3, 4, 5, 6, 7], [34.17, 60.78, 84.88, 613.28, 681.30, 1398.38]),
        500: ([2, 3, 4, 5, 6, 7], [31.91, 58.23, 77.6, 136.26, 974.79, 1106.92]),
        1000: ([2, 3, 4, 5, 6, 7], [32.04, 64.55, 124.19, 417.64, 718.57, 1127.69]),
    }

    N = runs[Lambda][0]
    time = runs[Lambda][1]
    time_log = np.log(runs[Lambda][1])
    popt, pcov = curve_fit(lambda x, a, b: np.log(func(x, a, b)), N, time_log)

    print(popt)
    xdata = np.linspace(2, 16)
    plt.plot(xdata, func(xdata, *popt))
    plt.scatter(N, time)
    plt.title(f'exponential is {popt[0]:.1f}exp({popt[1]:.2f}*x)')
    plt.yscale('log')

    plt.show()

    couplings = []
    for N_Flavor in N:
        with h5py.File(f'./data/couplings_Lambda={Lambda}_N={N_Flavor}.hdf5', "r") as f:
            couplings.append(f.attrs["coupling"])

    coupling_mf = Gross_Neveu.GN_2p1.calculate_one_g2(1.0, 1.0, Lambda)

    # def func(x, a, b):
    #     return coupling_mf - a*np.exp(-b*x)
    def func(x, a):
        return coupling_mf * np.arctan(a*x) * 2 / np.pi

    plt.axhline(y=coupling_mf, color='r', linestyle='-')
    popt, pcov = curve_fit(func, N, couplings)
    xdata = np.linspace(2, 16)
    plt.plot(xdata, func(xdata, *popt))
    plt.title(f'the factor inside arctan tangent is {popt[0]:.2f}')
    plt.scatter(N, couplings)

    plt.show()


def main():
    parser = argparse.ArgumentParser(
        prog='Program Name',
        description='What the program does',
        epilog='Text at the bottom of help')
    parser.add_argument(
        '-N', type=int, help='Set the number of flavors', required=True, default=None)
    parser.add_argument('-L', type=float,
                        help='Set the UV Cutoff', required=True, default=None)
    parser.add_argument('--time',
                        help='Set the UV Cutoff', action='store_false')

    args = parser.parse_args()

    if not args.time:
        # if args.L is None and args.N is None:
        #     raise RuntimeError("can't have args.N and ")
        compute_time_fits(args.L)
    else:
        compute_couping(args.L, args.N)


if __name__ == "__main__":
    # compute_time_fits()
    main()
