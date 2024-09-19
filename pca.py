from __future__ import annotations

import tkinter
from icecream import ic

from correlFuncs import *
from paramMatching import dispCorrMatr


def transIndep(sample_data):
    s_n = []
    for i in range(1, len(sample_data) + 1):
        str = f"Вибірка {i}"
        if sample_data[str]['var'].get() == 1:
            s_n.append(str)

    buff = np.array([sample_data[s_n[i]]["data"] for i in range(len(s_n))])
    x = buff[0]
    y = buff[1]

    avr_x = average(x)
    avr_y = average(y)

    tg_two_fi = (2 * correlationRelation(x, y) * std_err(x, avr_x) * std_err(y, avr_y)) / (
            std_err(x, avr_x) ** 2 - std_err(y, avr_y) ** 2)

    fi = np.arctan(tg_two_fi) / 2

    x_new = [x[i] * np.cos(fi) + y[i] * np.sin(fi) for i in range(len(x))]
    y_new = [-x[i] * np.sin(fi) + y[i] * np.cos(fi) for i in range(len(x))]

    new_coords = [x_new, y_new]

    for i in range(len(new_coords)):
        sample_data[s_n[i]]["data"] = new_coords[i]

    return fi


def backTransIndep(sample_data, fi):
    ic(fi)
    s_n = []
    for i in range(1, len(sample_data) + 1):
        str = f"Вибірка {i}"
        if sample_data[str]['var'].get() == 1:
            s_n.append(str)

    buff = np.array([sample_data[s_n[i]]["data"] for i in range(len(s_n))])

    x_new = buff[-2]
    y_new = buff[-1]

    x_back = [x_new[i] * np.cos(fi) - y_new[i] * np.sin(fi) for i in range(len(x_new))]
    y_back = [x_new[i] * np.sin(fi) + y_new[i] * np.cos(fi) for i in range(len(y_new))]

    new_coords = [x_back, y_back]

    for i in range(len(new_coords)):
        sample_data[s_n[i]]["data"] = new_coords[i]


def eigenValsVectors(sample_data):
    s_n = []
    for i in range(1, len(sample_data) + 1):
        str = f"Вибірка {i}"
        if sample_data[str]['var'].get() == 1:
            s_n.append(str)

    buff = np.array([sample_data[s_n[i]]["data"] for i in range(len(s_n))])

    disp_corr_matrix = dispCorrMatr(buff)

    eigenvalues, eigenvectors = np.linalg.eig(disp_corr_matrix)

    sorted_indices = np.argsort(eigenvalues)[::-1]
    sorted_eigenvalues = eigenvalues[sorted_indices]
    sorted_eigenvectors = eigenvectors[:, sorted_indices]

    return sorted_eigenvalues, sorted_eigenvectors


# def pca_forw_nDim(sample_data, sorted_eigenvectors, n, sample_menu, sample_checkbuttons):
#     s_n = []
#     for i in range(1, len(sample_data) + 1):
#         str = f"Вибірка {i}"
#         if sample_data[str]['var'].get() == 1:
#             s_n.append(str)
#
#     buff = np.array([sample_data[s_n[i]]["data"] for i in range(len(s_n))])
#
#     disp_corr_matrix = dispCorrMatr(buff)
#
#     buff_trans = np.transpose(sorted_eigenvectors)[:n] @ buff
#     ic(buff_trans)
#
#     for checkbutton in sample_checkbuttons:
#         checkbutton.destroy()
#     sample_checkbuttons.clear()
#     sample_menu.delete(0, 'end')
#
#     sample_data.clear()
#
#     for i in range(n):
#         sample_num = i + 1
#         sample_name = f"Вибірка {sample_num}"
#         sample_var = tkinter.IntVar()
#         sample_data[sample_name] = {"data": buff_trans[i], "var": sample_var}
#
#         sample_menu.add_checkbutton(label=sample_name, variable=sample_var)

def pca_forw_nDim(sample_data, sorted_eigenvectors, n, sample_menu, sample_checkbuttons):
    s_n = []
    for i in range(1, len(sample_data) + 1):
        str = f"Вибірка {i}"
        if sample_data[str]['var'].get() == 1:
            s_n.append(str)

    buff = np.array([sample_data[s_n[i]]["data"] for i in range(len(s_n))])

    ic(sorted_eigenvectors.shape)
    ic(np.transpose(sorted_eigenvectors)[:n])
    ic(np.transpose(sorted_eigenvectors[:n]).shape)
    ic(sorted_eigenvectors[:n].shape)
    ic(buff[:n].shape)
    # buff_trans = np.transpose(buff[:n]) @ np.transpose(sorted_eigenvectors)[:n]
    # buff_trans = np.transpose(sorted_eigenvectors[:n]) @ buff[:n]
    buff_trans = np.transpose(sorted_eigenvectors) @ buff
    ic(buff_trans)
    # buff_trans = np.transpose(buff_trans)
    for checkbutton in sample_checkbuttons:
        checkbutton.destroy()
    sample_checkbuttons.clear()
    sample_menu.delete(0, 'end')

    sample_data.clear()

    for i in range(len(buff)):
        sample_num = i + 1
        sample_name = f"Вибірка {sample_num}"
        sample_var = tkinter.IntVar()
        sample_data[sample_name] = {"data": buff_trans[i], "var": sample_var}

        sample_menu.add_checkbutton(label=sample_name, variable=sample_var)


def pca_back_nDim(sample_data, eigenvectors, n, sample_menu, sample_checkbuttons):
    s_n = []
    for i in range(1, len(sample_data) + 1):
        str = f"Вибірка {i}"
        if sample_data[str]['var'].get() == 1:
            s_n.append(str)

    buff = np.array([sample_data[s_n[i]]["data"] for i in range(len(s_n))])

    ic(eigenvectors)
    buff_back = np.transpose(eigenvectors[:n]) @ buff[:n]

    for checkbutton in sample_checkbuttons:
        checkbutton.destroy()
    sample_checkbuttons.clear()
    sample_menu.delete(0, 'end')

    sample_data.clear()

    for i in range(len(buff)):
        sample_num = i + 1
        sample_name = f"Вибірка {sample_num}"
        sample_var = tkinter.IntVar()
        sample_data[sample_name] = {"data": buff_back[i], "var": sample_var}

        sample_menu.add_checkbutton(label=sample_name, variable=sample_var)

