import h5py
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from scipy.interpolate import LinearNDInterpolator
from scipy.optimize import root
import numpy as np
from python_files.phase_diagram.file_io import get_filename
from mpl_toolkits.axes_grid1 import make_axes_locatable


def get_result(Lambda, N_Flavor):
    filename = f'./data/phase_diagram/phase_diagram_Lambda_{Lambda}_N_Flavor_{N_Flavor}.hdf5'
    with h5py.File(filename, "r") as f:
        result = f["phase_diagram"][:]
    return result


def density_plot(result):
    mu_max = np.max(result[:, 0])
    T_max = np.max(result[:, 1])
    # plot the results
    points = list(zip(result[:, 0], result[:, 1]))
    values = list(result[:, 2])
    interpolation = LinearNDInterpolator(points, values)
    X, Y = np.meshgrid(np.linspace(0.0, mu_max, 1000),
                       np.linspace(0.01, T_max, 1000))
    Z = interpolation(X, Y)
    plt.pcolormesh(X, Y, Z, shading='auto')
    plt.colorbar()


def get_mesh(mu_array, T_array):
    mu_max = np.max(mu_array)
    mu_min = np.min(mu_array)
    T_max = np.max(T_array)
    T_min = np.min(T_array)
    X, Y = np.meshgrid(np.arange(mu_min, mu_max, 0.001),
                       np.arange(T_min, T_max, 0.001))
    return X, Y


def get_phase_diagram(file: h5py.File):
    mu_array = file["phase_diagram"][:, 0]
    T_array = file["phase_diagram"][:, 1]
    return mu_array, T_array


def create_interpolation(mu_array, T_array, values):
    interpolation = LinearNDInterpolator(
        list(zip(mu_array, T_array)), values)
    return interpolation


def create_Z_data(mu_array, T_array, values, X, Y):
    interpolation = create_interpolation(mu_array, T_array, values)
    return interpolation(X, Y)


def compute_and_plot_value(mu_array, T_array, value, X, Y, ax: Axes, fig: plt.Figure, title: str = ""):
    Z = create_Z_data(mu_array, T_array, value, X, Y)
    im = ax.pcolormesh(X, Y, Z, shading='auto')
    ax.set_xlabel(r'$\mu$')
    ax.set_ylabel(r'$T$')
    ax.set_title(title)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')


def compute_liftschitz(mu_array, T_array, second_div_at_zero, forth_div_at_zero):
    interpolation_second_div_at_zero = create_interpolation(
        mu_array, T_array, second_div_at_zero)
    interpolation_forth_div_at_zero = create_interpolation(
        mu_array, T_array, forth_div_at_zero)

    def f(arg):
        x, y = arg
        return interpolation_second_div_at_zero(x, y), interpolation_forth_div_at_zero(x, y)

    tol = None
    x0 = [9.921e-01, 3.157e-02]
    result = root(f, x0=x0, tol=tol)
    if not result.success:
        x0 = [0.9, 0.2]
        result = root(f, x0=x0, tol=tol)

    print(result)
    return result.x


def plot_observables(Lambda, N_Flavor):
    filename = get_filename(Lambda, N_Flavor)
    fig, ((sigma_0_plot, mass_square_plot), (second_div_plot, forth_div_plot)) = plt.subplots(
        2, 2, tight_layout=True)
    with h5py.File(filename, "r") as file:
        mu_array, T_array = get_phase_diagram(file)
        X, Y = get_mesh(mu_array, T_array)

        mu_l, T_l = compute_liftschitz(
            mu_array, T_array, file["second_div_at_zero"][:], file["forth_div_at_zero"][:])

        sigmas = file["sigma_0"][:]
        compute_and_plot_value(
            mu_array, T_array, sigmas, X, Y, sigma_0_plot, fig, title=r'Condensate $\sigma_0$')
        sigma_0_plot.scatter([mu_l], [T_l])

        mass_squares = file["mass_square"][:]
        compute_and_plot_value(
            mu_array, T_array, mass_squares, X, Y, mass_square_plot, fig, title=r'$m_0^2$')

        second_div = 2*np.heaviside(file["second_div_at_zero"][:], 0) - 1
        compute_and_plot_value(
            mu_array, T_array, second_div, X, Y, second_div_plot, fig)
        second_div_plot.scatter([mu_l], [T_l])

        forth_div = 2*np.heaviside(file["forth_div_at_zero"][:], 0) - 1
        compute_and_plot_value(
            mu_array, T_array, forth_div, X, Y, forth_div_plot, fig)
        forth_div_plot.scatter([mu_l], [T_l])

        plt.show()
