import numpy as np
import matplotlib.pyplot as plt


def plot_phase_boundary(boundary_first_order, boundary_second_order):
    if len(boundary_first_order) != 0:
        plt.scatter(boundary_first_order["mu"],
                    boundary_first_order["T"])
    if len(boundary_second_order) != 0:
        plt.scatter(boundary_second_order["mu"],
                    boundary_second_order["T"])


def liftschitz_point(result, tolerance):
    boundary_first_order, boundary_second_order = compute_phase_boundary(
        result, tolerance)
    first_mu, first_T = boundary_first_order[-1]
    for mu, T in reversed(boundary_first_order):
        if first_mu != mu:
            break
        else:
            first_T = T
    second_mu, second_T = boundary_second_order[0]
    # print(boundary_first_order[-1], boundary_first_order[-2])
    print(boundary_second_order)

    return (first_mu + second_mu)/2, (first_T + second_T)/2


def plot_tolerance(result, tolerance_array):
    liftschitz_point_array = []
    for tolerance in tolerance_array:
        mu, T = liftschitz_point(result, tolerance)
        print(
            f'for the tolerance of {tolerance:.2f}, to liftschitz point is ({mu}, {T})')
        liftschitz_point_array.append([mu, T])

    liftschitz_point_array = np.array(liftschitz_point_array)
    mu_to_be = 0.999233
    T_to_be = 0.105838
    lfp_mf = np.array([mu_to_be, T_to_be])

    print('done')

    fig, (ax1, ax2) = plt.subplots(1, 2)
    ax1.plot(tolerance_array, liftschitz_point_array[:, 0])
    ax1.set_ylabel("mu")
    ax1.axhline(y=mu_to_be, color='r', linestyle='-')
    ax2.plot(tolerance_array, liftschitz_point_array[:, 1])
    ax2.axhline(y=T_to_be, color='r', linestyle='-')
    ax2.set_ylabel("T")
    plt.show()

    # liftschitz_point_array = np.array(liftschitz_point_array)
    #
    # error = [np.linalg.norm(lfp_mf - point)
    #          for point in liftschitz_point_array]
    # print('done')
    # plt.plot(tolerance_array, error)
    # plt.show()


def compute_phase_boundary(result, tolerance):
    mu_array = np.sort(np.unique(result[:, 0]))
    T_array = np.sort(np.unique(result[:, 1]))
    xv, yv = np.meshgrid(mu_array, T_array, indexing='ij')
    # create result array
    mesh = np.empty((len(mu_array), len(T_array)))

    len_T = len(T_array)
    for index, (mu, T, value) in enumerate(result):
        first_index = index // len_T
        second_index = index % len_T

        mesh[first_index, second_index] = value

    boundary_first_order = []
    boundary_second_order = []
    derivative_array = []

    def direction(i, j, dir: int):
        ii = 1 if dir == 1 else 0
        jj = 1 if dir == 0 else 0
        if (mesh[i, j] == 0 and mesh[i+ii, j+jj] != 0) or (mesh[i, j] != 0 and mesh[i+ii, j+jj] == 0):
            spatial_distance = (
                T_array[j] - T_array[j+jj]) if dir == 0 else (mu_array[i] - mu_array[i+ii])
            derivative = np.abs(
                (mesh[i+ii, j+jj] - mesh[i, j])/spatial_distance)

            derivative_array.append(derivative)

            if derivative > tolerance:
                boundary_second_order.append(((mu_array[i] + mu_array[i+ii]) /
                                             2, (T_array[j] + T_array[j+jj])/2))
            else:
                boundary_first_order.append(((mu_array[i] + mu_array[i+ii]) /
                                             2, (T_array[j] + T_array[j+jj])/2))

    for i, mu in enumerate(mu_array[:-1]):
        for j, T in enumerate(T_array[:-1]):
            direction(i, j, 1)
            direction(i, j, 0)

    dtype = [('mu', float), ('T', float)]

    boundary_first_order = np.array(boundary_first_order, dtype=dtype)
    boundary_second_order = np.array(boundary_second_order, dtype=dtype)

    boundary_first_order = np.sort(boundary_first_order, order=['mu', 'T'])
    boundary_second_order = np.flip(
        np.sort(boundary_second_order, order=['T', 'mu']))

    return boundary_first_order, boundary_second_order
