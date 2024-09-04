import numpy as np
from scipy.stats import chi2

from sndDim import *
from icecream import ic


def dispCorrMatr(buff):
    disp_corr_matrix = []
    buff_mean = [average(buff[i]) for i in range(len(buff))]
    buff_std_err = [std_err(buff[i], buff_mean[i]) for i in range(len(buff))]
    ic(buff_mean)
    ic(buff_std_err)
    for i in range(len(buff)):
        row = []
        for j in range(len(buff)):
            if j == i:
                row.append(round(buff_std_err[j] ** 2, 4))
            else:
                row.append(round(buff_std_err[i] * buff_std_err[j] * linearCorrelation(buff[i], buff[j]), 4))
        disp_corr_matrix.append(row)
    return disp_corr_matrix


def checkMatching_dc(samplesGroup):
    g_n = []
    for i in range(1, len(samplesGroup) + 1):
        str = f"Набір {i}"
        if samplesGroup[str]['var'].get() == 1:
            g_n.append(str)
    ic(samplesGroup)
    groupBuff = np.array([samplesGroup[g_n[i]]["samples_group"] for i in range(len(g_n))])
    ic(groupBuff)
    k = len(groupBuff)

    """matching of variance and covariance matrices"""
    dc_matrices = np.array([dispCorrMatr(groupBuff[i]) for i in range(k)])
    N = len(groupBuff) * len(groupBuff[0]) * len(groupBuff[0][0])
    Nd = len(groupBuff[0][0])
    S = (1 / (N - k)) * sum([np.dot((Nd - 1), dc_matrices[i]) for i in range(k)])

    V = sum([((Nd - 1) / 2) * np.log(np.linalg.det(S) / np.linalg.det(dc_matrices[i])) for i in range(k)])
    ic(V)


    plt.figure(figsize=(12, 8))

    n = len(groupBuff[0])

    chi_quant = chi2.ppf(1 - .05, df=n * (n + 1) * (k - 1) / 2)

    plt.title('Перевірка збігу середніх')

    plt.text(0.5, 0.9, 'Збіг дисперсійно-коваріаційних матриць', fontsize=14, ha='center')


    if V <= chi_quant:
        plt.text(0.5, 0.8, f'V = {V:.4f} <= {chi_quant:.4f} --> дисперсійно-коваріаційні матриці збігаються', fontsize=14,
                 ha='center')
        chi_quant_1 = chi2.ppf(1 - .05, df=n)
        V1 = checkMatching_mean1(samplesGroup)
        plt.text(0.5, 0.7, 'Перевірка рівності двох багатовимірних середніх при рівних дк-матрицях', fontsize=14,
                 ha='center')

        if V1 <= chi_quant_1:
            plt.text(0.5, 0.6, f'V = {V1:.4f} <= {chi_quant_1:.4f} --> середні рівні', fontsize=14, ha='center')
        else:
            plt.text(0.5, 0.6, f'V = {V1:.4f} > {chi_quant_1:.4f} --> середні не рівні', fontsize=14, ha='center')
    else:
        plt.text(0.5, 0.8, f'V = {V:.4f} > {chi_quant:.4f} --> дисперсійно-коваріаційні матриці не збігаються', fontsize=14,
                 ha='center')
        chi_quant_2 = chi2.ppf(1 - .05, df=n * (k - 1))
        V2 = checkMatching_mean2(samplesGroup)
        plt.text(0.5, 0.7, 'Перевірка рівності k n-вимірних середніх при розбіжних дк-матрицях', fontsize=14,
                 ha='center')

        if V2 <= chi_quant_2:
            plt.text(0.5, 0.6, f'V = {V2:.4f} <= {chi_quant_2:.4f} --> середні рівні', fontsize=14, ha='center')
        else:
            plt.text(0.5, 0.6, f'V = {V2:.4f} > {chi_quant_2:.4f} --> середні не рівні', fontsize=14, ha='center')


    plt.show()


def checkMatching_mean1(samplesGroup):
    g_n = []
    for i in range(1, len(samplesGroup) + 1):
        str = f"Набір {i}"
        if samplesGroup[str]['var'].get() == 1:
            g_n.append(str)

    groupBuff = np.array([samplesGroup[g_n[i]]["samples_group"] for i in range(len(g_n))])
    k = len(groupBuff)

    if k != 2:
        return None

    dc_matrices = np.array([dispCorrMatr(groupBuff[i]) for i in range(k)])
    V = -(len(groupBuff[0][0]) * 2 - 2 - (len(groupBuff[0]) / 2)) * np.log(
        np.linalg.det(dc_matrices[1]) / np.linalg.det(dc_matrices[0]))
    return V
    ic(V)


def checkMatching_mean2(samplesGroup):
    g_n = []
    for i in range(1, len(samplesGroup) + 1):
        str = f"Набір {i}"
        if samplesGroup[str]['var'].get() == 1:
            g_n.append(str)

    groupBuff = np.array([samplesGroup[g_n[i]]["samples_group"] for i in range(len(g_n))])
    k = len(groupBuff)
    ic(k)
    Nd = len(groupBuff[0][0])
    # ic(len(groupBuff[0]))
    # ic(len(groupBuff[0][0]))
    dc_matrices = np.array([dispCorrMatr(groupBuff[i]) for i in range(k)])

    groupMean = np.array([[np.mean(groupBuff[i][j]) for j in range(len(groupBuff[0]))] for i in range(k)])

    # groupMean = np.transpose(groupMean)
    Xl = np.array([np.transpose(groupBuff[i]) for i in range(k)])
    ic(Xl)
    xd = np.array([[sum(Xl[j][i]) / Nd for i in range(Nd)] for j in range(k)])
    ic(xd)

    ic(Xl[0])

    diffXlxd = np.transpose(Xl[0]) - xd
    Sd = (1 / (Nd - 1)) * (diffXlxd @ np.transpose(diffXlxd))
    ic(Sd)
    x_hat = np.dot(Nd, Sd) @ (np.dot(np.dot(Nd, np.linalg.inv(Sd)), xd))
    ic(x_hat)


def checkMatching_mean2(samplesGroup):
    g_n = []
    for i in range(1, len(samplesGroup) + 1):
        str = f"Набір {i}"
        if samplesGroup[str]['var'].get() == 1:
            g_n.append(str)

    groupBuff = np.array([samplesGroup[g_n[i]]["samples_group"] for i in range(len(g_n))])

    n = len(groupBuff[0])
    k = len(groupBuff)

    xMeanD = []
    Nd = []
    Sd = []
    xMean = None

    listForSd = []
    listForxMeanD = []

    for index in range(len(groupBuff)):
        Nd.append(len(groupBuff[index][0]))
        ValueN = Nd[-1]

        for j in range(n):
            listForxMeanD.append(sum(groupBuff[index][j]) / ValueN)

        xMeanD.append(np.array(listForxMeanD))
        listForxMeanD.clear()

        matrixForSd = np.zeros((n, n))
        for i in range(ValueN):
            for j in range(n):
                listForSd.append(groupBuff[index][j][i])

            vectorForSd = np.array(listForSd)
            matrixForSd += (1.0 / (ValueN - 1.0)) * np.outer(vectorForSd - xMeanD[-1], vectorForSd - xMeanD[-1])
            listForSd.clear()

        Sd.append(matrixForSd)

    sum1X = np.zeros((n, n))
    sum2X = np.zeros((n, 1))

    for i in range(len(groupBuff)):
        sum1X += Nd[i] * np.linalg.inv(Sd[i])
        sum2X += Nd[i] * np.dot(np.linalg.inv(Sd[i]), xMeanD[i].reshape(-1, 1))

    xMean = np.dot(np.linalg.inv(sum1X), sum2X).flatten()

    V = 0

    for i in range(k):
        V += Nd[i] * np.dot((xMeanD[i] - xMean).reshape(1, -1),
                            np.dot(np.linalg.inv(Sd[i]), (xMeanD[i] - xMean).reshape(-1, 1)))

    return V[0][0]
