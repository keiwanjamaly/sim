import h5py
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import UnivariateSpline
from scipy import signal
from scipy.optimize import root
import numpy as np
from python_files.phase_diagram.file_io import get_filename
from mpl_toolkits.axes_grid1 import make_axes_locatable
from shapely import geometry
from skimage import measure


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
    X, Y = np.meshgrid(np.arange(mu_min, mu_max, 0.01),
                       np.arange(T_min, T_max, 0.01))
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


def compute_and_density_plot_value(mu_array, T_array, value, X, Y, ax: Axes, fig: plt.Figure, title: str = ""):
    Z = create_Z_data(mu_array, T_array, value, X, Y)
    im = ax.pcolormesh(X, Y, Z, shading='auto')
    ax.set_xlabel(r'$\mu$')
    ax.set_ylabel(r'$T$')
    ax.set_title(title)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(im, cax=cax, orientation='vertical')


def compute_and_contour_plot_value(mu_array, T_array, value, X, Y, ax: Axes, fig: plt.Figure, title: str = "", color=None):
    Z = create_Z_data(mu_array, T_array, value, X, Y)
    im = ax.contour(X, Y, Z, 0, colors=[color] if color is not None else None)
    ax.set_xlabel(r'$\mu$')
    ax.set_ylabel(r'$T$', rotation=0)
    ax.set_title(title)
    return im


def get_path(contour):
    paths = [path for collection in contour.collections
             for path in collection.get_paths()]
    path_index = np.argmax([len(path) for path in paths])
    return paths[path_index]


def findIntersection(contour1, contour2):
    p1 = get_path(contour1)
    p2 = get_path(contour2)
    v1 = p1.vertices
    v2 = p2.vertices

    poly1 = geometry.LineString(v1)
    poly2 = geometry.LineString(v2)

    intersection = poly1.intersection(poly2)

    print(intersection)
    return intersection.x, intersection.y


def plot_observables(Lambda, N_Flavor, save):
    filename = get_filename(Lambda, N_Flavor)
    fig, ((sigma_0_plot, mass_square_plot), (div_plot, phase_boundary_plot)) = plt.subplots(
        2, 2, tight_layout=True)
    with h5py.File(filename, "r") as file:
        mu_array, T_array = get_phase_diagram(file)
        X, Y = get_mesh(mu_array, T_array)

        second_div = file["second_div_at_zero"][:]
        im_1 = compute_and_contour_plot_value(
            mu_array, T_array, second_div, X, Y, div_plot, fig, color="red")
        forth_div = file["forth_div_at_zero"][:]
        im_2 = compute_and_contour_plot_value(
            mu_array, T_array, forth_div, X, Y, div_plot, fig, color="blue")
        mu_l, T_l = findIntersection(im_1, im_2)
        div_plot.plot([mu_l], [T_l], 'ro')

        sigmas = file["sigma_0"][:]
        compute_and_density_plot_value(
            mu_array, T_array, sigmas, X, Y, sigma_0_plot, fig, title=r'Condensate $\sigma_0$')
        sigma_0_plot.plot([mu_l], [T_l], 'rx')

        im_c = compute_and_contour_plot_value(
            mu_array, T_array, sigmas, X, Y, phase_boundary_plot, fig)
        phase_boundary_plot.plot([mu_l], [T_l], 'rx')

        mass_squares = file["mass_square"][:]
        compute_and_density_plot_value(
            mu_array, T_array, mass_squares, X, Y, mass_square_plot, fig, title=r'$m_0^2$')

        if save is not None:
            points = get_path(im_c).vertices
            for (mu, T) in points:
                print(f'{mu}\t{T}', file=save)
        else:
            plt.show()
