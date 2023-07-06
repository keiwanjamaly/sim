from python_files.gross_neveu.couplings.compute_couplings import compute_all_coupling_fits
from python_files.gross_neveu.couplings.compute_couplings import compute_couping_fit
from python_files.gross_neveu.couplings.compute_couplings import compute_couping
import argparse

def main():
    parser = argparse.ArgumentParser(prog='couplings calculator')
    parser.add_argument(
        '-N', type=int, help='Set the number of flavors', required=True, default=None)
    parser.add_argument('-L', type=float,
                        help='Set the UV Cutoff', required=True, default=None)
    parser.add_argument('--plot', action='store_true')
    parser.add_argument('--fit',
                        help='extraplate coupling data', action='store_true')
    parser.add_argument('--all',
                        help='extraplate coupling data', action='store_true')

    args = parser.parse_args()

    if args.fit:
        if args.all:
            compute_all_coupling_fits()
        else:
            compute_couping_fit(args.L, True)
    else:
        compute_couping(args.L, args.N)

if __name__ == "__main__":
    main()
