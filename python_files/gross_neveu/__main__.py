from python_files.gross_neveu.Gross_Neveu import GN_2p1, GN_1p1
from python_files.grid_creator import create_homogenious_grid


def main():
    sigma_max = 6.0
    N_Grid = 1000
    Lambda = 1e5
    kir = 1e-4
    samples = 10
    mu = 0.1
    T = 0.1
    N_Flavor = 2
    grid_points = create_homogenious_grid(sigma_max, N_Grid)
    sol = GN_1p1(grid_points, Lambda, kir,
                 samples, mu, T, N_Flavor=N_Flavor)
    pass


if __name__ == "__main__":
    main()
