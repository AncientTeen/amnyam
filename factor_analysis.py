from __future__ import annotations

from correl_regress_Mult import *


def generality(buff: array[array[float]]) -> None:
    corrMatr = corrMatrix(buff)

    """Maximum correlation method"""
    h_1 = []
    for i in range(len(corrMatr)):
        max_coff = 0
        for j in range(len(corrMatr[i])):
            if i != j:
                if max_coff < abs(corrMatr[i][j]):
                    max_coff = abs(corrMatr[i][j])
            else:
                continue
        h_1.append(max_coff)


    """The triad method"""
    h_2 = []
    for i in range(len(corrMatr)):
        max_coff_1 = 0
        max_coff_2 = 0
        k, q = 0, 0
        for j in range(len(corrMatr[i])):
            if i != j:
                if max_coff_1 < abs(corrMatr[i][j]):
                    max_coff_2 = max_coff_1
                    q = k
                    max_coff_1 = abs(corrMatr[i][j])
                    k = j
            else:
                continue
        h_2.append(abs((max_coff_2 * max_coff_1) / corrMatr[k, q]))

    """Average method"""
    h_3 = []
    for i in range(len(corrMatr)):
        avr = 0
        for j in range(len(corrMatr[i])):
            if i != j:
                avr += abs(corrMatr[i][j])
            else:
                continue
        h_3.append(avr / (len(corrMatr) - 1))
    """Centroid method"""
    h_4 = []

    corrMatr_cpy = corrMatr.copy()

    max_elem = []
    for i, el in enumerate(corrMatr):
        el = np.delete(el, i)
        max_elem.append(max(el))

    for i in range(len(corrMatr_cpy)):
        corrMatr_cpy[i][i] = max_elem[i]

    sum_all = 0
    for i in range(len(corrMatr_cpy)):
        for j in range(len(corrMatr_cpy[i])):
            sum_all += abs(corrMatr_cpy[i][j])

    for i in range(len(corrMatr_cpy)):
        sum_row = 0
        for j in range(len(corrMatr_cpy[i])):
            sum_row += abs(corrMatr_cpy[i][j])
        h_4.append((sum_row ** 2) / sum_all)


    """Averroids method"""
    h_5 = []

    sum_all = 0
    for i in range(len(corrMatr)):
        for j in range(len(corrMatr[i])):
            if i != j:
                sum_all += abs(corrMatr[i][j])
            else:
                continue

    for i in range(len(corrMatr)):
        sum_row = 0
        for j in range(len(corrMatr[i])):
            if i != j:
                sum_row += abs(corrMatr[i][j])
            else:
                continue
        h_5.append(((sum_row ** 2) * len(corrMatr)) / (sum_all * (len(corrMatr) - 1)))

    """PCA method"""
    h_6 = []
    eigenvalues, eigenvectors = np.linalg.eig(corrMatr)

    sorted_indices = np.argsort(eigenvalues)[::-1]
    sorted_eigenvalues = eigenvalues[sorted_indices]
    sorted_eigenvectors = eigenvectors[:, sorted_indices]

    w = 0
    for i in range(len(sorted_eigenvalues)):
        if sorted_eigenvalues[i] > 1:
            w += 1

    for i in range(len(corrMatr)):
        sum_sq = 0
        for j in range(w):
            sum_sq += sorted_eigenvectors[j][i] ** 2
        h_6.append(sum_sq)


    ic(h_1)
    ic(h_2)
    ic(h_3)
    ic(h_4)
    ic(h_5)
    ic(h_6)

    pass



