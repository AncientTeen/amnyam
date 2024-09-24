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

    a, b, c, d = calculate_plane_equation(buff_back)
    print(f"Plane equation: {a}x + {b}y + {c}z + {d} = 0")

    A, a_null = multRegr(buff_back)
    print(f"Regression coefficients: {A}")
    print(f"Intercept: {a_null}")



def calculate_plane_equation(data):
    # Separate the data into X (first two components) and y (third component)
    X = data[:, :2]
    y = data[:, 2]

    # Fit a linear regression model
    model = LinearRegression()
    model.fit(X, y)

    # Get the coefficients and intercept
    a, b = model.coef_
    c = -1  # We assume the plane equation is in the form ax + by + z + d = 0
    d = model.intercept_

    return a, b, c, d

def multRegr(buff, y_sample=1, bound=[1, 10]):
    Y = buff[y_sample - 1]
    X = np.delete(buff, y_sample - 1, 0)
    y_avr = average(Y)
    x_avr = [average(X[i]) for i in range(len(X))]

    Y_null = [Y[i] - y_avr for i in range(len(Y))]
    X_null = [[(X[i][j] - x_avr[i]) for j in range(len(X[0]))] for i in range(len(X))]
    A = np.linalg.inv(X_null @ np.transpose(X_null)) @ X_null @ np.transpose(Y_null)
    a_null = y_avr - sum([A[i] * x_avr[i] for i in range(len(X))])

    e = [Y[i] - a_null - sum([A[j] * X[j][i] for j in range(len(A))]) for i in range(len(Y))]
    S_zal = (1 / (len(Y) - len(buff))) * sum([e[i] ** 2 for i in range(len(e))])
    C = np.linalg.inv(X @ np.transpose(X))
    Y_sq = std_err(Y, y_avr)
    X_sq = [std_err(X[i], x_avr[i]) for i in range(len(X))]

    # Y_low, Y_hat, Y_up = multRegrConfInt(Y, X, A, a_null, C, S_zal)

    return A, a_null
