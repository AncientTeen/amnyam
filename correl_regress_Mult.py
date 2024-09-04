import numpy as np
from sndDim import *
import pingouin as pg
import pandas as pd


def corrMatrix(buff):
    mtrx = np.zeros((len(buff), len(buff)))
    for i in range(len(buff)):
        for j in range(len(buff)):
            if j == i:
                mtrx[i][j] = 1
            else:
                mtrx[i][j] = linearCorrelation(buff[i], buff[j])

    return mtrx


def partCorrCoff(buff):
    transposed_buff = list(map(list, zip(*buff)))

    df = pd.DataFrame(transposed_buff)
    partial_corr_coefficients = {}

    for i in range(len(df.columns)):
        for j in range(i + 1, len(df.columns)):
            sample1 = df.columns[i]
            sample2 = df.columns[j]

            covariates = [s for s in df.columns if s != sample1 and s != sample2]

            partial_corr_result = pg.partial_corr(data=df, x=sample1, y=sample2, covar=covariates)

            key = f"{i + 1} та {j + 1}"
            partial_corr_coefficients[key] = partial_corr_result['r'].values[0]

    return partial_corr_coefficients


def confIntr_partCorrCoff(buff):
    transposed_buff = list(map(list, zip(*buff)))

    df = pd.DataFrame(transposed_buff)
    conf_partial_corr_coefficients = {}

    for i in range(len(df.columns)):
        for j in range(i + 1, len(df.columns)):
            sample1 = df.columns[i]
            sample2 = df.columns[j]

            covariates = [s for s in df.columns if s != sample1 and s != sample2]

            # Calculate partial correlation coefficient and its confidence interval
            conf_partial_corr_result = pg.partial_corr(data=df, x=sample1, y=sample2, covar=covariates)

            # Extract confidence interval values
            ci_values = conf_partial_corr_result['CI95%'].values

            # Store confidence interval values in the dictionary
            key = f"{i + 1} та {j + 1}"
            conf_partial_corr_coefficients[key] = ci_values

    return conf_partial_corr_coefficients


def multCorrCoff(matrix):
    mltCorrCoff = []
    for i in range(len(matrix)):
        mod_matrix = np.delete(matrix, i, 0)
        mod_matrix = np.delete(mod_matrix, i, 1)
        mltCorrCoff.append(np.sqrt(1 - np.linalg.det(matrix) / np.linalg.det(mod_matrix)))

    return mltCorrCoff


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



    Y_low, Y_hat, Y_up = multRegrConfInt(Y, X, A, a_null, C, S_zal)

    return A, a_null, C, S_zal, Y_sq, X_sq, Y_low, Y_hat, Y_up, e, Y


def multRegrConfInt(Y, X, A, a_null, C, S_zal):
    if len(Y) > 10 and len(Y) < 120:
        t_quant = 2.105
    elif len(Y) >= 120:
        t_quant = 1.96

    Y_hat = [(a_null + sum([A[i] * X[i][j] for i in range(len(A))])) for j in range(len(Y))]
    X = np.transpose(X)

    Y_low = [Y_hat[i] - t_quant * np.sqrt(S_zal) * np.sqrt(1 + X[i] @ C @ np.transpose(X[i])) for i in range(len(Y))]
    Y_up = [Y_hat[i] + t_quant * np.sqrt(S_zal) * np.sqrt(1 + X[i] @ C @ np.transpose(X[i])) for i in range(len(Y))]



    return Y_low, Y_hat, Y_up
