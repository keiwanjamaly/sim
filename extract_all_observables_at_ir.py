import glob
import h5py
import numpy as np

files = glob.glob("./mean_field_test/Lambda=10/*")
outfile = "./mean_field_test/Lambda=10.hdf5"
mu = np.empty(len(files))
T = np.empty(len(files))
sigma = np.empty(len(files))
firstDiv = np.empty(len(files))
thirdDiv = np.empty(len(files))

for i, filename in enumerate(files):
    with h5py.File(filename, "r") as f:
        mu[i] = f.attrs["mu"]
        T[i] = f.attrs["T"]
        sigma[i] = f["sigma"][-1]
        firstDiv[i] = f["first_div"][-1]
        thirdDiv[i] = f["third_div"][-1]
        print(i, "/", len(files))

mu_unique = np.sort(np.unique(mu))
T_unique = np.sort(np.unique(T))
X, Y = np.meshgrid(mu_unique, T_unique)
sigma_grid = np.empty((len(mu_unique), len(T_unique)))

for mu_i in range(len(mu_unique)):
    for T_i in range(len(T_unique)):
        mu_filter = mu == X[mu_i, T_i]
        T_filter = T == Y[mu_i, T_i]
        filter_combined = np.logical_and(mu_filter, T_filter)
        sigma_grid[mu_i, T_i] = sigma[filter_combined]

boundary_x = []
boundary_y = []
for mu_i in range(len(mu_unique)-1):
    for T_i in range(len(T_unique)-1):
        current_element = sigma_grid[mu_i, T_i]
        ce_x = X[mu_i, T_i]
        ce_y = Y[mu_i, T_i]
        # print current point
        # print(
        #     f'currently working on mu={X[mu_i,T_i]}, T={Y[mu_i,T_i]}, with sigma={sigma_grid[mu_i,T_i]}')
        # handle north direction
        if current_element != 0 and sigma_grid[mu_i, T_i+1] == 0:
            boundary_x.append((X[mu_i, T_i+1] + ce_x)/2)
            boundary_y.append((Y[mu_i, T_i+1] + ce_y)/2)
        # print(
        #     f'north mu={X[mu_i,T_i+1]}, T={Y[mu_i,T_i+1]}, with sigma={sigma_grid[mu_i,T_i+1]}'
        # )
        # print east
        if current_element != 0 and sigma_grid[mu_i+1, T_i] == 0:
            boundary_x.append((X[mu_i+1, T_i] + ce_x)/2)
            boundary_y.append((Y[mu_i+1, T_i] + ce_y)/2)
        # print(
        #     f'north mu={X[mu_i+1,T_i]}, T={Y[mu_i+1,T_i]}, with sigma={sigma_grid[mu_i+1,T_i]}'
        # )
        # print north-west
        if current_element != 0 and sigma_grid[mu_i+1, T_i+1] == 0:
            boundary_x.append((X[mu_i+1, T_i+1] + ce_x)/2)
            boundary_y.append((Y[mu_i+1, T_i+1] + ce_y)/2)
        # print(
        #     f'north mu={X[mu_i+1,T_i+1]}, T={Y[mu_i+1,T_i+1]}, with sigma={sigma_grid[mu_i+1,T_i+1]}'
        # )

# print(list(zip(boundary_x, boundary_y)))

with h5py.File(outfile, "w") as f:
    f.create_dataset("sigma", data=sigma)
    f.create_dataset("third_div", data=thirdDiv)
    f.create_dataset("first_div", data=firstDiv)
    f.create_dataset("T", data=T)
    f.create_dataset("mu", data=mu)
    f.create_dataset("boundary_mu", data=boundary_x)
    f.create_dataset("boundary_T", data=boundary_y)
