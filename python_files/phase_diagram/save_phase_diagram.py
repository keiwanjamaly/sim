import h5py
from python_files.data_class import Potential
from swinging_door import swinging_door


def save_final_phase_diagram_to_file(filename, outname):
    with h5py.File(filename, "r") as file:
        grid = file['grid'][:]
        mu_array = file['phase_diagram'][:, 0]
        T_array = file['phase_diagram'][:, 1]
        u_array = file['u'][:]
        potentials = [Potential(grid, u) for u in u_array]
        sigma_array = [potential.sigma for potential in potentials]
        mass_square_array = [potential.mass_square for potential in potentials]
        pressure_array = [
            potential.pressure for potential in potentials]
        mu_back = mu_array[0]
        points = []
        tmp = []
        for mu, T, *arguments in zip(mu_array, T_array, sigma_array, mass_square_array, pressure_array):
            if mu_back != mu:
                # result = list(swinging_door(tmp, deviation=0.01))
                result = tmp[::5]
                for T_back, *arguments_back in result:
                    print(mu_back, T_back, *arguments_back,
                          sep='\t', file=outname)
                print(sep='\t', file=outname)
                mu_back = mu
                tmp = []

            tmp.append((T, *arguments))

        # for mu, T, sigma, mass_square, pressure in zip(mu_array, T_array, sigma_array, mass_square_array, pressure_array):
        #     if mu_back != mu:
        #         mu_back = mu
        #         print(sep='\t', file=outname)
        #     print(mu, T, sigma, mass_square, pressure,
        #           sep='\t', file=outname)
