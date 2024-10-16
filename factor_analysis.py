from __future__ import annotations

import numpy as np
from correl_regress_Mult import *


def generality(matr: list[list[float]]) -> tuple[list[list[float]], float, list[float]]:
    # def generality(matr: list[list[float]]) -> None:

    """Maximum correlation method"""
    h_1 = []
    for i in range(len(matr)):
        max_coff = 0
        for j in range(len(matr[i])):
            if i != j:
                if max_coff < abs(matr[i][j]):
                    max_coff = abs(matr[i][j])
            else:
                continue
        h_1.append(max_coff)

    """The triad method"""
    h_2 = []
    for i in range(len(matr)):
        max_coff_1, max_coff_2 = 0, 0
        k, q = -1, -1  # Indices for the two max values

        # Loop through the row to find two max elements (excluding diagonal)
        for j in range(len(matr[i])):
            if i != j:  # Skip the diagonal
                value = abs(matr[i][j])

                # Update the largest and second largest values
                if value > max_coff_1:
                    max_coff_2 = max_coff_1
                    q = k  # Store previous max index
                    max_coff_1 = value
                    k = j  # Update index of the new max
                elif value > max_coff_2:
                    max_coff_2 = value
                    q = j  # Update index of second max
        h_2.append(abs((max_coff_2 * max_coff_1) / matr[k, q]))

    """Average method"""
    h_3 = []
    for i in range(len(matr)):
        avr = 0
        for j in range(len(matr[i])):
            if i != j:
                avr += abs(matr[i][j])
            else:
                continue
        h_3.append(avr / (len(matr) - 1))

    """Centroid method"""
    h_4 = []

    matr_cpy = matr.copy()

    max_elem = []
    for i, el in enumerate(matr_cpy):
        el = np.delete(el, i)
        max_elem.append(max(abs(el)))

    for i in range(len(matr_cpy)):
        matr_cpy[i][i] = max_elem[i]

    sum_all = 0
    for i in range(len(matr_cpy)):
        for j in range(len(matr_cpy[i])):
            sum_all += abs(matr_cpy[i][j])

    for i in range(len(matr_cpy)):
        sum_row = 0
        for j in range(len(matr_cpy[i])):
            sum_row += abs(matr_cpy[i][j])
        h_4.append((sum_row ** 2) / sum_all)

    """Averroids method"""
    h_5 = []
    sum_all = 0
    for i in range(len(matr)):
        for j in range(len(matr[i])):
            if i != j:
                sum_all += abs(matr[i][j])
            else:
                continue

    for i in range(len(matr)):
        sum_row = 0
        for j in range(len(matr[i])):
            if i != j:
                sum_row += abs(matr[i][j])
            else:
                continue
        h_5.append(((sum_row ** 2) * len(matr)) / (sum_all * (len(matr) - 1)))

    """PCA method"""
    h_6 = []
    eigenvalues, eigenvectors = np.linalg.eig(matr)

    sorted_indices = np.argsort(eigenvalues)[::-1]
    sorted_eigenvalues = eigenvalues[sorted_indices]
    sorted_eigenvectors = eigenvectors[:, sorted_indices]
    # ic(sorted_eigenvalues)
    # ic(np.mean(sorted_eigenvalues))
    w = 0
    for i in range(len(sorted_eigenvalues)):
        if sorted_eigenvalues[i] >= 0.9:
            w += 1

    for i in range(len(matr)):
        sum_sq = 0
        for j in range(w):
            sum_sq += sorted_eigenvectors[j][i] ** 2
        h_6.append(sum_sq)

    # ic(h_1)
    # ic(h_2)
    # ic(h_3)
    # ic(h_4)
    # ic(h_5)
    # ic(h_6)

    h_buff = np.array([h_1, h_2, h_3, h_4, h_5, h_6])
    # h_temp = []
    # for i in range(len(h_buff)):
    #     if max(h_buff[i]) > 1:
    #         continue
    #     else:
    #         h_temp.append(h_buff[i])
    # ic(h_temp)
    # h_buff = h_temp

    f_min_buff = []

    for i in range(len(h_buff)):
        matr_cpy = matr.copy()
        f_min_temp = 0

        for j in range(len(matr)):
            matr_cpy[j][j] = h_buff[i][j]

        eigenvalues, eigenvectors = np.linalg.eig(matr_cpy)

        sorted_indices = np.argsort(eigenvalues)[::-1]
        sorted_eigenvectors = eigenvectors[:, sorted_indices]

        # A = np.transpose(sorted_eigenvectors)
        A = sorted_eigenvectors
        # A_t = sorted_eigenvectors
        A_t = np.transpose(sorted_eigenvectors)

        R_gen = matr_cpy - A @ A_t
        # ic(R_gen)
        for k in range(len(R_gen)):
            for q in range(len(R_gen[k])):
                if k != q:
                    f_min_temp += R_gen[k][q] ** 2

        f_min_buff.append(f_min_temp)

    f_min = min(f_min_buff)
    index_min = np.argmin(f_min_buff)

    h = h_buff[index_min]
    matr_cpy = matr.copy()
    for j in range(len(matr)):
        matr_cpy[j][j] = h[j]

    return matr_cpy, f_min, h


def factor_anal(corr_matr: list[list[float]]) -> tuple[list[list[float]], int, int]:
    itr = 0
    R_matr_buff = []
    f_min_buff = []
    h_buff = []
    A_buff = []
    eigenvals_buff = []
    w = 0

    run = True
    while run:
        if itr < 1:
            R_matr, f_min, h = generality(corr_matr)
            eigenvalues, eigenvectors = np.linalg.eig(R_matr)
            sorted_indices = np.argsort(eigenvalues)[::-1]
            sorted_eigenvectors = eigenvectors[:, sorted_indices]
            A = sorted_eigenvectors
            A_buff.append(A)
            eigenvals_buff.append(eigenvalues)
        else:
            R_matr, f_min, h = generality(R_matr_buff[-1])
            eigenvalues, eigenvectors = np.linalg.eig(R_matr)
            sorted_indices = np.argsort(eigenvalues)[::-1]
            sorted_eigenvectors = eigenvectors[:, sorted_indices]
            # ic(eigenvalues)
            w = len([num for num in eigenvalues if num >= 0])
            A = sorted_eigenvectors
            A_buff.append(A)
            eigenvals_buff.append(eigenvalues)


        if len(R_matr_buff) < 2:
            R_matr_buff.append(R_matr)
            f_min_buff.append(f_min)
            h_buff.append(h)

        else:
            """check_1"""
            f_check = f_min_buff[1] - f_min_buff[0]

            """check_2"""
            a_check = 0
            for i in range(len(A_buff[0])):
                for j in range(len(A_buff[0][0])):
                    a_check += (A_buff[1][i][j] - A_buff[0][i][j]) ** 2

            # """check_3"""

            if f_check > 0 or a_check < 0.0001:
                run = False
            else:
                R_matr_buff[0] = R_matr_buff[1]
                R_matr_buff[1] = R_matr
                f_min_buff[0] = f_min_buff[1]
                f_min_buff[1] = f_min
                h_buff[0] = h_buff[1]
                h_buff[1] = h
                A_buff[0] = A_buff[1]
                A_buff[1] = A
                eigenvals_buff[0] = eigenvals_buff[1]
                eigenvals_buff[1] = eigenvalues


        itr += 1

    A = A_buff[1]
    eigenvals_arr = eigenvals_buff[1]
    # ic(A)
    return A, w, itr, eigenvals_arr
