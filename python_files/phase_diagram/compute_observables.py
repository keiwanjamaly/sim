from python_files.phase_diagram.file_io import get_filename
from scipy.optimize import root
import h5py
from python_files.data_class import Potential
from python_files.phase_diagram.plot_phase_diagram import create_interpolation


def compute_sigma(file: h5py.File, potentials: list[Potential]):
    sigmas = [potential.sigma for potential in potentials]
    if 'sigma_0' not in file.keys():
        file.create_dataset('sigma_0', (len(sigmas)), dtype=float)
    file['sigma_0'][:] = sigmas[:]


def compute_mass_square(file: h5py.File, potentials: list[Potential]):
    masses = [potential.mass_square for potential in potentials]
    if 'mass_square' not in file.keys():
        file.create_dataset('mass_square', (len(masses)), dtype=float)
    file['mass_square'][:] = masses[:]


def compute_second_div_at_zero(file: h5py.File, potentials: list[Potential]):
    second_divs_at_zero = [
        potential.second_div_at_zero for potential in potentials]
    if 'second_div_at_zero' not in file.keys():
        file.create_dataset('second_div_at_zero',
                            (len(second_divs_at_zero)), dtype=float)
    file['second_div_at_zero'][:] = second_divs_at_zero


def compute_forth_div_at_zero(file: h5py.File, potentials: list[Potential]):
    forth_divs_at_zero = [
        potential.forth_div_at_zero for potential in potentials]
    if 'forth_div_at_zero' not in file.keys():
        file.create_dataset('forth_div_at_zero',
                            (len(forth_divs_at_zero)), dtype=float)
    file['forth_div_at_zero'][:] = forth_divs_at_zero


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


def compute_observables(Lambda, N_Flavor):
    filename = get_filename(Lambda, N_Flavor)
    with h5py.File(filename, 'r+') as file:
        phase_diagram = file["phase_diagram"][:]
        us = file["u"][:]
        grid = file["grid"]
        potentials: list[Potential] = [
            Potential(grid, u) for u in us]
        compute_sigma(file, potentials)
        compute_mass_square(file, potentials)
        compute_second_div_at_zero(file, potentials)
        compute_forth_div_at_zero(file, potentials)
