import numpy as np
from joblib import Parallel, delayed
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
    print(result)
    return result


def calculate_parameter(N_Flavor, sigma_0, Lambda, kir, sigma_max, delta_sigma, h):
    model = Gross_Neveu.GN_2p1
    one_over_g2 = model.calculate_one_g2(h, sigma_0, Lambda)
    result = newton(calculate_sigma, one_over_g2, args=(
        model, sigma_max, Lambda, kir, delta_sigma, N_Flavor, h, sigma_0))
    print(f'N = {N_Flavor}, 1/g^2 = {result}')
    return (N_Flavor, result)


def main():
    Lambda = 1000
    kir = 1e-2
    h = 1
    sigma_0 = 1.0
    sigma_max = 2000.0
    delta_sigma = 0.006

    N_Flavor_list = range(2, 3)

    job_list = []
    for N_Flavor in N_Flavor_list:
        job_list.append(delayed(calculate_parameter)(
            N_Flavor, sigma_0, Lambda, kir, sigma_max, delta_sigma, h))

    result = Parallel(n_jobs=1)(job_list)
    result = np.array(result)
    print(result)

    for N_Flavor in N_Flavor_list:
        with h5py.File(f'./data/couplings_Lambda={Lambda}_N={N_Flavor}.hdf5', "w") as f:
            f.attrs["sigma_max"] = sigma_max
            f.attrs["Lambda"] = Lambda
            f.attrs["kir"] = kir

    # with h5py.File(path_and_filename, "w") as f:
    #     f.attrs["sigmaMax"] = self.grid.sigma[-1]
    #     f.attrs["Lambda"] = self.Lambda
    #     f.attrs["kir"] = self.kir
    #     f.attrs["NGrid"] = len(self.grid.sigma)
    #     f.attrs["computation_time"] = self.time_elapsed
    #     f.attrs["tolerance"] = self.tolerance
    #     f.attrs["grid_style"] = type(self.grid).__name__
    #     f.attrs["extrapolation_order"] = self.grid.extrapolation
    #     # TODO: include manual storage of additional attrs
    #     for key, value in self.file_attributes.items():
    #         f.attrs[key] = value
    #
    #     # create datasets and store them
    #     f.create_dataset("t", data=observables["t"])
    #     f.create_dataset("sigma", data=observables["sigma"])
    #     f.create_dataset("massSquare", data=observables["massSquare"])
    #     f.create_dataset("first_div", data=observables["first_div"])
    #     f.create_dataset("third_div", data=observables["third_div"])
    #     f.create_dataset("pressure", data=observables["pressure"])
    #     f.create_dataset("k", data=observables["k"])
    #
    #     # create additional datasets and store them
    #     if self.save_flow_flag:
    #         f.create_dataset("grid", data=self.grid.sigma)
    #         y_dset = f.create_dataset("y", (observables["y"].to_numpy().shape[0], len(
    #             self.grid.sigma)))
    #         Q_dset = f.create_dataset("Q", (observables["y"].to_numpy().shape[0], len(
    #             self.grid.sigma)))
    #         S_dset = f.create_dataset("S", (observables["y"].to_numpy().shape[0], len(
    #             self.grid.sigma)))
    #
    #         for i, elem in enumerate(observables["y"].to_numpy()):
    #             y_dset[i, :] = elem[:]
    #             k = observables["k"][i]
    #             Q_dset[i, :] = observables["diffusion"][i]
    #             S_dset[i, :] = observables["source"][i]
    #
    # return path_and_filename


if __name__ == "__main__":
    main()
