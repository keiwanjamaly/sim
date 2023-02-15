import glob
import h5py
import numpy as np

files = glob.glob("./mean_field_test/Lambda=10/*")
outfile = "./mean_field_test/Lambda=10.hdf5"

mu = np.empty(len(files))
T = np.empty(len(files))


for i, filename in enumerate(files):
    with h5py.File(filename, "r") as f:
        mu[i] = f.attrs["mu"]
        T[i] = f.attrs["T"]
        print(i, "/", len(files))

mu_unique = np.sort(np.unique(mu))
T_unique = np.sort(np.unique(T))
mu_and_T = np.column_stack((mu, T))
X, Y = np.meshgrid(mu_unique, T_unique)
sigma_grid = np.empty((len(mu_unique), len(T_unique)))

print(mu_unique.shape, T_unique.shape)
print(mu_unique)
print(T_unique)

for mu_i in range(len(mu_unique)):
    for T_i in range(len(T_unique)):
        exists = ([X[mu_i, T_i], Y[mu_i, T_i]] == mu_and_T).all(1).any()
        if not exists:
            print((X[mu_i, T_i], Y[mu_i, T_i]))
