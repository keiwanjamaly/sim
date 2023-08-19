import argparse
from python_files.numeric_tests.grid_resolution_test import compute_grid_spacing_tests
from python_files.numeric_tests.Lambda_max_test import Lambda_max_test
from python_files.numeric_tests.sigma_max_test import sigma_max_test
from python_files.numeric_tests.infrared_cutoff_test import kir_max_test
from python_files.numeric_tests.inhomogenious_grid_test import test_inhomogenious_grid


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--save', help='set the output file type', type=argparse.FileType('w', encoding='latin-1'))
    parser.add_argument(
        '--spacing', help='compute grid spacing', action='store_true')
    parser.add_argument(
        '--max', help='compute sigma max spacing', action='store_true')
    parser.add_argument(
        '--Lambda', help='compute Lambda max', action='store_true')
    parser.add_argument(
        '--kir', help='test infrared cutoff', action='store_true')
    parser.add_argument(
        '--grid', help='test an inhomogenious grid', action='store_true')
    parser.add_argument(
        '--T', help='set the temperature', required=True, type=float)
    parser.add_argument(
        '--mu', help='set the chemical potential', required=True, type=float)
    args = parser.parse_args()

    mu = args.mu
    T = args.T

    if args.spacing:
        compute_grid_spacing_tests(mu, T, args.save)
    if args.Lambda:
        Lambda_max_test(mu, T, args.save)
    if args.max:
        sigma_max_test(mu, T, args.save)
    if args.kir:
        kir_max_test(mu, T, args.save)
    if args.grid:
        test_inhomogenious_grid(mu, T, args.save)


    # print(args.save)
if __name__ == "__main__":
    main()
