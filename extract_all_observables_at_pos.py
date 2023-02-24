import glob
import h5py
import numpy as np

Lambda = 1e5
kir = 1e-4

tir = - np.log(kir/Lambda)

# n = 10

positions_k = np.array([0.4, 0.3, 0.2, 0.1, 0.05, 0.03, 0.01, 0.003])
positions = - np.log(positions_k/Lambda)
# positions = np.linspace(0, tir, n)
# positions_k = Lambda*np.exp(-positions)

files = glob.glob("./phase_diagram/*")
outfile = "./phase_diagram=2.hdf5"
mu = np.empty(len(files))
T = np.empty(len(files))
sigma = np.empty((len(positions), len(files)))
firstDiv = np.empty((len(positions), len(files)))
thirdDiv = np.empty((len(positions), len(files)))


def find_interval(t_test, f, starting_index=0):
    t_array = f["t"][starting_index:]
    for i in range(len(t_array)-1):
        if t_test >= t_array[i] and t_test <= t_array[i+1]:
            return (i+starting_index, i+1+starting_index, t_array[i], t_array[i+1])


def interpolate_observable(min, t_min, max, t_max, file, key, t):
    return np.interp(t, [t_min, t_max], [f[key][min], f[key][max]])


for i, filename in enumerate(files):
    with h5py.File(filename, "r") as f:
        mu[i] = f.attrs["mu"]
        T[i] = f.attrs["T"]
        backup_starting_index = 0
        for t_pos, t in enumerate(positions):
            t_arg_min, t_arg_max, t_min, t_max = find_interval(
                t, f, backup_starting_index)
            backup_starting_index = t_arg_max
            sigma[t_pos, i] = interpolate_observable(
                t_arg_min, t_min, t_arg_max, t_max, f, "sigma", t)
            firstDiv[t_pos, i] = interpolate_observable(
                t_arg_min, t_min, t_arg_max, t_max, f, "first_div", t)
            thirdDiv[t_pos, i] = interpolate_observable(
                t_arg_min, t_min, t_arg_max, t_max, f, "third_div", t)
    print(i+1, "/", len(files))


def calculate_boundary(mu_b, T_b, sigma_b):
    mu_unique = np.sort(np.unique(mu_b))
    T_unique = np.sort(np.unique(T_b))
    X, Y = np.meshgrid(mu_unique, T_unique, indexing="ij")
    sigma_grid = np.empty((len(mu_unique), len(T_unique)))

    for mu_i in range(len(mu_unique)):
        for T_i in range(len(T_unique)):
            mu_filter = mu_b == X[mu_i, T_i]
            T_filter = T_b == Y[mu_i, T_i]
            filter_combined = np.logical_and(mu_filter, T_filter)
            sigma_grid[mu_i, T_i] = sigma_b[filter_combined]

    boundary_x = []
    boundary_y = []
    for mu_i in range(len(mu_unique)-1):
        for T_i in range(len(T_unique)-1):
            # handle current point
            current_element = sigma_grid[mu_i, T_i]
            ce_x = X[mu_i, T_i]
            ce_y = Y[mu_i, T_i]
            # handle north
            if current_element != 0 and sigma_grid[mu_i, T_i+1] == 0:
                boundary_x.append((X[mu_i, T_i+1] + ce_x)/2)
                boundary_y.append((Y[mu_i, T_i+1] + ce_y)/2)
            # handle east
            if current_element != 0 and sigma_grid[mu_i+1, T_i] == 0:
                boundary_x.append((X[mu_i+1, T_i] + ce_x)/2)
                boundary_y.append((Y[mu_i+1, T_i] + ce_y)/2)
            # handle north-west
            if current_element != 0 and sigma_grid[mu_i+1, T_i+1] == 0:
                boundary_x.append((X[mu_i+1, T_i+1] + ce_x)/2)
                boundary_y.append((Y[mu_i+1, T_i+1] + ce_y)/2)

    return boundary_x, boundary_y


boundaries_x = []
boundaries_y = []

for i, t in enumerate(positions):
    b_x, b_y = calculate_boundary(mu, T, sigma[i])
    boundaries_x.append(b_x)
    boundaries_y.append(b_y)

with h5py.File(outfile, "w") as f:
    f.create_dataset("T", data=T)
    f.create_dataset("mu", data=mu)
    f.create_dataset("t", data=positions)
    f.create_dataset("k", data=positions_k)
    t_grp = f.create_group('t_grp')
    k_grp = f.create_group('k_grp')
    for i, (t, k) in enumerate(zip(positions, positions_k)):
        grp = t_grp.create_group(str(t))
        k_grp[str(k)] = grp
        grp.create_dataset("sigma", data=sigma[i])
        grp.create_dataset("third_div", data=thirdDiv[i])
        grp.create_dataset("first_div", data=firstDiv[i])
        grp.create_dataset("boundary_mu", data=boundaries_x[i])
        grp.create_dataset("boundary_T", data=boundaries_y[i])
