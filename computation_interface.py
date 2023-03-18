import ctypes
import os


def main():
    os.system(
        'gcc -fPIC -shared -o computation_library.so computation_library.c')
    lib = ctypes.CDLL("./computation_library.so")

    initialize = lib.initialize
    initialize.argtypes = [ctypes.c_double,
                           ctypes.c_double, ctypes.c_double, ctypes.c_double]
    initialize(1e5, 1e-4, 0.0, 0.1)
    lib.compute()
    lib.cleanup()


if __name__ == "__main__":
    main()
