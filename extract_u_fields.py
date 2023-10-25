import h5py
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-f', help='filename of the h5py file', type=str, required=True)
    parser.add_argument(
        '-T', help='Temperature to select from', type=float, required=True)
    parser.add_argument(
        '-mu', help='chemical potential so select from', type=float, required=True)
    args = parser.parse_args()

    with h5py.File(args.f, 'r') as f:
        mu_array = f['phase_diagram'][:, 0]
        T_array = f['phase_diagram'][:, 1]
        for i, (mu, T) in enumerate(zip(mu_array, T_array)):
            if mu == args.mu and T == args.T:
                u_array = f['u'][i, :]
                grid = f['grid'][:]
                for sigma, u in zip(grid, u_array):
                    print(sigma, u, sep='\t')


if __name__ == "__main__":
    main()
