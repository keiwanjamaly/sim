from python_files.gross_neveu.Gross_Neveu import get_model
from python_files.data_class import DataClass
from python_files.grid_creator import create_inhomogenious_grid_from_cell_spacing


def sigma(one_over_g2: float, dimension: int, mu: float, T:float, sigma_max, Lambda, kir,
          delta_sigma, N_Flavor, h, sigma_0):
    samples = 3
    grid_points = create_inhomogenious_grid_from_cell_spacing(
        sigma_max, delta_sigma)
    model = get_model(dimension)
    model = model(grid_points, Lambda, kir, samples,
                  mu, T, N_Flavor, h, one_over_g2, sigma_0)

    y = model.return_data.solution
    x = model.return_data.grid
    time = model.return_data.time

    data_class = DataClass(x, time, y)
    sigma_0_ir = data_class.sigma_0[-1]

    return sigma_0_ir
