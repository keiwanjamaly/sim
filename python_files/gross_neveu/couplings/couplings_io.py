import h5py
import os
import numpy as np
import re
from python_files.gross_neveu.Gross_Neveu import get_model

def get_all_cutoffs_with_more_than_two_couplings(dir: str):
    all_files = os.listdir(dir)
    filtered_files = filter(lambda x: x.startswith(f'couplings_Lambda='), all_files)
    Lambdas = map(lambda x: float(re.search("[0-9]*\\.[0-9]*", x)[0]), filtered_files)
    return list(set(Lambdas))

def get_computed_couplings_from_file(Lambda: float, dir: str):
    all_files = os.listdir(dir)
    filtered_files = filter(lambda x: x.startswith(f'couplings_Lambda={Lambda}_'), all_files)
    couplings = []
    flavors = []
    for filename in filtered_files:
        with h5py.File(os.path.join(dir, filename)) as f:
            couplings.append(f.attrs["coupling"])
            flavors.append(f.attrs["N_Flavor"])

    return flavors, couplings


def save_coupling(Lambda, coupling, N_Flavor, time, dir):
    filename = generate_filename(Lambda, N_Flavor, dir)
    with h5py.File(filename, "w") as f:
        f.attrs["Lambda"] = Lambda
        f.attrs["coupling"] = coupling
        f.attrs["N_Flavor"] = N_Flavor
        f.attrs["time"] = time


def get_closest_coupling_approximation_or_MF_coupling(Lambda: float, N_Flavor: int, dir, dimension: int = 2):
    try:
        return get_closest_coupling_approximation(Lambda, N_Flavor, dir)
    except (NoCouplingComputed, CouplingNotInFile):
        model = get_model(dimension)
        print("using MF coupling instead")
        return model.calculate_one_g2(h=1, sigma_0=1, Lambda=Lambda)


def get_closest_coupling_approximation(Lambda: float, N_Flavor: int, dir: str):
    filename = generate_filename(Lambda, N_Flavor, dir)
    if file_exists(filename):
        return get_exact_coupling_from_file(filename)
    elif file_exists(os.path.join(dir, 'couping_pre_computations.hdf5')):
        try:
            return get_pre_computed_coupling_from_file(Lambda, N_Flavor, dir)
        except CouplingNotInFile as e:
            raise e
    else:
        raise NoCouplingComputed()


class CouplingNotInFile(Exception):
    pass


class NoCouplingComputed(Exception):
    pass


def coupling_fit_function(x, a, coupling_mf):
    return coupling_mf * np.arctan(a*x) * 2 / np.pi


def get_pre_computed_coupling_from_file(Lambda: float, N_Flavor: float, dir: str):
    with h5py.File(os.path.join(dir, 'couping_pre_computations.hdf5'), "r") as f:
        if Lambda in f["couplings"][:, 0]:
            pos = list(f["couplings"][:, 0]).index(Lambda)
            beta, alpha = f["couplings"][pos][1:]
            return alpha * np.arctan(beta * N_Flavor) * 2 / np.pi
        else:
            raise CouplingNotInFile()


def generate_filename(Lambda: float, N_Flavor: float, dir: str):
    filename = os.path.join(dir, f'couplings_Lambda={Lambda}_N={N_Flavor}.hdf5')
    return filename


def file_exists(filename: str):
    return os.path.isfile(filename)


def get_exact_coupling_from_file(filename: str):
    with h5py.File(filename, "r") as f:
        return f.attrs["coupling"]
