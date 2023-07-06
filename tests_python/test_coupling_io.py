from python_files.gross_neveu.couplings.couplings_io import get_computed_couplings_from_file
from python_files.gross_neveu.couplings.couplings_io import get_all_cutoffs_with_more_than_two_couplings


def test_computed_couplings_from_file():
    flavours, couplings = get_computed_couplings_from_file(10.0, "./tests_python/data")
    assert flavours == [9, 10]
    assert couplings == [1.2937570902210294, 1.2940825471463189]

def test_get_all_couplings_from_files():
    Lambdas = get_all_cutoffs_with_more_than_two_couplings("./tests_python/data")
    assert Lambdas == [10.0, 20.0]

