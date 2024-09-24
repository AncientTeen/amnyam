from __future__ import annotations

import tkinter
from icecream import ic
from sklearn.linear_model import LinearRegression


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

    ic(eigenvalues)
    ic(eigenvectors)

    sorted_indices = np.argsort(eigenvalues)[::-1]
    sorted_eigenvalues = eigenvalues[sorted_indices]
    sorted_eigenvectors = eigenvectors[:, sorted_indices]

    return sorted_eigenvalues, sorted_eigenvectors


def pca_forw_nDim(sample_data, sorted_eigenvectors, sample_menu, sample_checkbuttons):
    s_n = []
    for i in range(1, len(sample_data) + 1):
        str = f"Вибірка {i}"
        if sample_data[str]['var'].get() == 1:
            s_n.append(str)

    buff = np.array([sample_data[s_n[i]]["data"] for i in range(len(s_n))])

    ic(sorted_eigenvectors)

    buff_trans = np.transpose(sorted_eigenvectors) @ buff
    # buff_trans = np.transpose(np.transpose(buff) @ sorted_eigenvectors)
    ic(buff_trans)
    ic(buff_trans.shape)
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
    # buff_back =  np.transpose(np.transpose(buff[:n]) @ eigenvectors[:n])
    buff_back = np.transpose(np.transpose(eigenvectors)[:n]) @ buff[:n]
    # buff_back = np.transpose(np.transpose(buff[:n]) @ eigenvectors[:n])

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


def plane_equation_from_points(buff):
    x1, y1, z1 = buff[0][0], buff[1][0], buff[2][0]
    x2, y2, z2 = buff[0][1], buff[1][1], buff[2][1]
    x3, y3, z3 = buff[0][2], buff[1][2], buff[2][2]
    a1 = x2 - x1
    b1 = y2 - y1
    c1 = z2 - z1
    a2 = x3 - x1
    b2 = y3 - y1
    c2 = z3 - z1
    a = b1 * c2 - b2 * c1
    b = a2 * c1 - a1 * c2
    c = a1 * b2 - b1 * a2
    d = (- a * x1 - b * y1 - c * z1)
    return a, b, c, d
