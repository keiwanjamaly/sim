import fileinput
import matplotlib.pyplot as plt
import numpy as np
from python_files.data_class import Potential


def main():
    sigma_list = []
    u_list = []

    for line in fileinput.input():
        field = line.replace('\n', '')
        sigma, u = field.split('\t')
        sigma = float(sigma)
        u = float(u)
        sigma_list.append(sigma)
        u_list.append(u)

    sigma_array = np.array(sigma_list)
    u_array = np.array(u_list)
    potential = Potential(sigma_array, u_array)

    x = sigma_array
    xmax_show = 1.2
    y_show = np.array(u_array)[x <= xmax_show]
    x_show = x[x <= xmax_show]

    plt.plot(x_show, y_show)
    plt.title(f'$\\sigma = {potential.sigma}$')
    plt.show()


if __name__ == "__main__":
    main()
