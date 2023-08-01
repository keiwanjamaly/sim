from python_files.phase_diagram.file_io import get_filename
from matplotlib.lines import Line2D
from python_files.phase_diagram.plot_phase_diagram import compute_and_contour_plot_value
from python_files.phase_diagram.plot_phase_diagram import get_phase_diagram
from python_files.phase_diagram.plot_phase_diagram import get_mesh
import numpy as np
import matplotlib.pyplot as plt
import h5py
import matplotlib
plt.rcParams.update({'font.size': 15})


def convert_float_to_tex(number: float):
    exp = np.log10(number)
    before = number/10**exp
    if before == 1:
        return f'10^{{{exp:.0f}}}'
    return f'{before:.1f} \\cdot 10^{{{exp:.0f}}}'


def parse_flavour_number(N_Flavor):
    return N_Flavor if not np.isinf(N_Flavor) and N_Flavor != -1 else r'\infty'


def plot_boundaries(Lambda, ax, fig):
    N_Flavor_List = [2, 4, 8, np.inf]
    cmap = matplotlib.cm.get_cmap('Spectral')
    colors = cmap(np.linspace(0, 1, len(N_Flavor_List)))
    legend_elements = []
    for N_Flavor, color in zip(reversed(N_Flavor_List), colors):
        filename = get_filename(Lambda, N_Flavor)
        with h5py.File(filename, "r") as file:
            mu_array, T_array = get_phase_diagram(file)
            X, Y = get_mesh(mu_array, T_array)
            sigmas = file["sigma_0"][:]
            compute_and_contour_plot_value(
                mu_array, T_array, sigmas, X, Y, ax, fig, color=[tuple(color)])
            label_number = parse_flavour_number(N_Flavor)
            legend_elements.append(
                Line2D([0], [0], color=tuple(color), ls='-', lw=2, label=f'$N = {label_number}$'))
    ax.legend(handles=legend_elements, loc='upper right')
    ax.set_ylim([0.01, 0.8])
    ax.set_xlim([0.00, 1.1])
    ax.set_title(
        f'$\\Lambda = {convert_float_to_tex(Lambda)} \\rightarrow k_{{ir}} = {convert_float_to_tex(1e-2)}$')
    ax.grid()


def plot_liftschitz_points():
    lf_points = {
        2: {10: (0.9202967208651033, 0.1634152881800319),
            20: (0.9517595970576876, 0.1309504047777122),
            200: (0.9904523549485027, 0.0293433431158443),
            500: (0.9906491646819718, 0.0293041785351767),
            1000: (0.9907597005949752, 0.0292814439061708),
            1500: (0.9908057454792941, 0.0292698988639711),
            2000: (0.9908307820448028, 0.0292668647755322)},
        4: {10: (0.9202505088215505, 0.1903632034844181),
            20: (0.9601522725665624, 0.152512077371942),
            200: (0.9905494129304002, 0.0304874666972015),
            500: (0.9908982202237707, 0.029781653717461),
            1000: (0.9911350354917215, 0.0296489931760655),
            1500: (0.9912462861892417, 0.0295403312802001),
            2000: (0.9913104883478119, 0.0294774908575736)},
        8: {10: (0.9201228293767697, 0.2012607484244602),
            20: (0.9602379619695697, 0.1708635144519281),
            200: (0.9906244014679839, 0.0803209451971805),
            500: (0.9911598598648622, 0.0307844989479824),
            1000: (0.9916049095517139, 0.0299461404861025),
            1500: (0.9918385234428402, 0.0297179153084838),
            2000: (0.9919828264981188, 0.029577165140295)},
        16: {10: (0.9116381175473001, 0.21146272063925),
             20: (0.9602054490502934, 0.1807033356790888),
             200: (0.9906770443170007, 0.110307274394645),
             500: (0.9913661133413946, 0.0693887300753366),
             1000: (0.9920528932253558, 0.0315724204091339),
             1500: (0.9924560677126901, 0.0306775734626419),
             2000: (0.9927236300115263, 0.0300832476714848)},
        -1: {10: (0.917344558057696, 0.2272274280093364),
             20: (0.9596167218117782, 0.1957693958530816),
             200: (0.9961214284318028, 0.1306153145428108),
             500: (0.9984594036825176, 0.114695050994166),
             1000: (0.9992328107842567, 0.1050434828225895),
             1500: (0.9994896813570081, 0.1007227471276981),
             2000: (0.999617679507465, 0.0971748268834355)}}

    keys = [10, 20, 200, 500, 1000, 1500, 2000]
    flavors = [2, 4, 8, 16, -1]

    fig, (axx, axy) = plt.subplots(
        1, 2, tight_layout=True, dpi=100, figsize=(12, 6))

    for n_flavor in flavors:
        points_mu = []
        points_T = []
        for key in keys:
            mu, T = np.array(lf_points[n_flavor][key])
            points_mu.append(1-mu)
            points_T.append(T)

        label_number = parse_flavour_number(n_flavor)
        axx.plot(keys, points_mu, '-x', label=f'$N = {label_number}$')
        axy.plot(keys, points_T, '-x', label=f'$N = {label_number}$')
        axx.set_title(r'X-Component of Lifshitz point $\mathbf{1-\mu}$')
        axy.set_title(r'Y-Component of Lifshitz point $\mathbf{T}$')
        axx.set_xlabel(r'$\Lambda$')
        axy.set_xlabel(r'$\Lambda$')
        axx.set_yscale('log')
        axy.set_yscale('log')
        axx.legend()
        axy.legend()

    # plt.rcParams["font.family"] = "serif"
    # plt.rcParams["font.serif"] = ["Times New Roman"]
    plt.show()


def plots_poster():
    Lambda = 10.0
    fig, ax = plt.subplots(
        1, 1, tight_layout=True, dpi=100, figsize=(8, 6))
    plot_boundaries(Lambda, ax, fig)
    Lambda = 1000.0
    fig, ax = plt.subplots(
        1, 1, tight_layout=True, figsize=(8, 6), dpi=100)
    plot_boundaries(Lambda, ax, fig)
    plt.show()
    plot_liftschitz_points()
