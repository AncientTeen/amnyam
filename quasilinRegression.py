import numpy as np

from funcs import *
import math

def quasilinRegr(x, y, fig1, ax1):
    n = len(x)

    fi_average = fi_avr(x, y)
    psi_average = psi_avr(y)
    fi_average_sq = fi_avr_sq(x, y)
    fi_psi_averageProd = fi_psi_avrProd(x, y)

    B_1 = (fi_psi_averageProd - fi_average * psi_average) / (fi_average_sq - fi_average ** 2)
    A_1 = np.exp(psi_average - B_1 * fi_average)

    z = np.log(y)
    z_avr = np.mean(z)
    t = x
    t_avr = np.mean(t)
    t_avr_sq = np.mean(t ** 2)
    t_z_prod = np.mean(t * z)
    sdt_err_t = np.std(t, ddof=1)
    sdt_err_z = np.std(z, ddof=1)
    correl = np.corrcoef(t, z)[0, 1]
    B_corr = correl * (sdt_err_z / sdt_err_t)
    A_corr = np.exp(z_avr - B_corr * t_avr)

    t_matr = np.array([[1, t_avr],
                       [t_avr, t_avr_sq]])

    inverse_t = np.linalg.inv(t_matr)

    t_z_matr = np.array([z_avr, t_z_prod])

    A_B_matr = np.dot(inverse_t, t_z_matr)
    A = np.exp(A_B_matr[0])
    B = A_B_matr[1]



    print(f"a (Matrix Approach): {A}, b (Matrix Approach): {B}")
    print(f"a (Weight Coefficient Approach): {A_1}, b (Weight Coefficient Approach): {B_1}")
    print(f"a (Correlation Approach): {A_corr}, b (Correlation Approach): {B_corr}")

    y_regression = [(A * np.exp(B * x[i])) for i in range(n)]
    ax1.plot(x, y_regression, c='red')

    fig1.canvas.draw()
    fig1.canvas.flush_events()


def confInterQuasilinRegr(x, y, fig1, ax1):
    n = len(x)
    t_quant = 0
    if n <= 30:
        t_quant = 2.05
    elif n > 30 and n <= 60:
        t_quant = 2.025
    elif n > 60 and n <= 120:
        t_quant = 1.99
    elif n > 120:
        t_quant = 1.96

    z = np.log(y)
    z_avr = np.mean(z)
    t = x
    t_avr = np.mean(t)
    t_avr_sq = np.mean(t ** 2)
    t_z_prod = np.mean(t * z)
    sdt_err_t = np.std(t, ddof=1)
    sdt_err_z = np.std(z, ddof=1)
    correl = np.corrcoef(t, z)[0, 1]


    B = correl * (sdt_err_z / sdt_err_t)
    A = np.exp(z_avr - B * t_avr)
    y_regression = [(A * np.exp(B * x[i])) for i in range(n)]

    S_zal = (1 / (n - 2)) * sum([(z[i] - y_regression[i]) ** 2 for i in range(n)])
    # S_zal = (1 / (n - 2)) * sum([(z[i] - (A + B * t[i])) ** 2 for i in range(n)])

    S_a = np.sqrt(S_zal) * np.sqrt(1 / n + t_avr ** 2 / (sdt_err_t ** 2 * (n - 1)))
    S_b = np.sqrt(S_zal) / (sdt_err_t * np.sqrt(n - 1))
    S_z = [np.sqrt(S_zal * (1 / n) + S_b ** 2 * (t[i] - t_avr) ** 2) for i in range(n)]
    y_min = [y_regression[i] - t_quant * S_z[i] for i in range(n)]
    y_max = [y_regression[i] + t_quant * S_z[i] for i in range(n)]
    # y_min = [np.exp(y_regression[i] - t_quant * S_z[i]) for i in range(n)]
    # y_max = [np.exp(y_regression[i] + t_quant * S_z[i]) for i in range(n)]
    ax1.plot(x, y_min, c='blue')
    ax1.plot(x, y_max, c='blue')

    fig1.canvas.draw()
    fig1.canvas.flush_events()



def toleranceLimQuasilinRegr(x, y, fig1, ax1):
    n = len(x)
    t_quant = 0
    if n <= 30:
        t_quant = 2.05
    elif n > 30 and n <= 60:
        t_quant = 2.025
    elif n > 60 and n <= 120:
        t_quant = 1.99
    elif n > 120:
        t_quant = 1.96

    z = np.log(y)
    z_avr = np.mean(z)
    t = x
    t_avr = np.mean(t)
    t_avr_sq = np.mean(t ** 2)
    t_z_prod = np.mean(t * z)
    sdt_err_t = np.std(t, ddof=1)
    sdt_err_z = np.std(z, ddof=1)
    correl = np.corrcoef(t, z)[0, 1]


    B = correl * (sdt_err_z / sdt_err_t)
    A = np.exp(z_avr - B * t_avr)
    y_regression = [(A * np.exp(B * x[i])) for i in range(n)]

    S_zal = np.sqrt((1 / (n - 2)) * sum([(z[i] - y_regression[i]) ** 2 for i in range(n)]))
    # S_zal = (1 / (n - 2)) * sum([(z[i] - (A + B * t[i])) ** 2 for i in range(n)])

    y_min = [y_regression[i] - t_quant * S_zal for i in range(n)]
    y_max = [y_regression[i] + t_quant * S_zal for i in range(n)]
    # y_min = [np.exp(y_regression[i] - t_quant * S_zal) for i in range(n)]
    # y_max = [np.exp(y_regression[i] + t_quant * S_zal) for i in range(n)]

    ax1.plot(x, y_min, c='green')
    ax1.plot(x, y_max, c='green')

    fig1.canvas.draw()
    fig1.canvas.flush_events()

def predictNewObsQuasilinRegr(x, y, fig1, ax1):
    n = len(x)
    t_quant = 0
    if n <= 30:
        t_quant = 2.05
    elif n > 30 and n <= 60:
        t_quant = 2.025
    elif n > 60 and n <= 120:n = len(x)
    t_quant = 0
    if n <= 30:
        t_quant = 2.05
    elif n > 30 and n <= 60:
        t_quant = 2.025
    elif n > 60 and n <= 120:
        t_quant = 1.99
    elif n > 120:
        t_quant = 1.96

    z = np.log(y)
    z_avr = np.mean(z)
    t = x
    t_avr = np.mean(t)
    t_avr_sq = np.mean(t ** 2)
    t_z_prod = np.mean(t * z)
    sdt_err_t = np.std(t, ddof=1)
    sdt_err_z = np.std(z, ddof=1)
    correl = np.corrcoef(t, z)[0, 1]


    B = correl * (sdt_err_z / sdt_err_t)
    A = np.exp(z_avr - B * t_avr)
    y_regression = [(A * np.exp(B * x[i])) for i in range(n)]

    S_zal = np.sqrt((1 / (n - 2)) * sum([(z[i] - y_regression[i]) ** 2 for i in range(n)]))
    # S_zal = (1 / (n - 2)) * sum([(z[i] - (A + B * t[i])) ** 2 for i in range(n)])

    S_b = np.sqrt(S_zal) / (sdt_err_t * np.sqrt(n - 1))
    S_z = [np.sqrt(S_zal ** 2 * (1 + 1 / n) + S_b ** 2 * (t[i] - t_avr) ** 2) for i in range(n)]
    y_min = [y_regression[i] - t_quant * S_z[i] for i in range(n)]
    y_max = [y_regression[i] + t_quant * S_z[i] for i in range(n)]

    ax1.plot(x, y_min, c='black')
    ax1.plot(x, y_max, c='black')

    fig1.canvas.draw()
    fig1.canvas.flush_events()


def inervalABQuaslinRegr(x, y):
    n = len(x)
    t_quant = 0
    if n <= 30:
        t_quant = 2.05
    elif n > 30 and n <= 60:
        t_quant = 2.025
    elif n > 60 and n <= 120:
        t_quant = 1.99
    elif n > 120:
        t_quant = 1.96

    z = np.log(y)
    z_avr = np.mean(z)
    t = x
    t_avr = np.mean(t)
    t_avr_sq = np.mean(t ** 2)
    t_z_prod = np.mean(t * z)
    sdt_err_t = np.std(t, ddof=1)
    sdt_err_z = np.std(z, ddof=1)
    correl = np.corrcoef(t, z)[0, 1]
    b = correl * (sdt_err_z / sdt_err_t)
    a = np.exp(z_avr - b * t_avr)
    y_regression = [(a * np.exp(b * x[i])) for i in range(n)]

    # S_zal = (1 / (n - 2)) * sum([(z[i] - y_regression[i]) ** 2 for i in range(n)])
    S_zal = (1 / (n - 2)) * sum([(z[i] - (a + b * t[i])) ** 2 for i in range(n)])
    S_a = np.sqrt(S_zal) * np.sqrt(1 / n + t_avr ** 2 / (sdt_err_t ** 2 * (n - 1)))
    S_b = np.sqrt(S_zal) / (sdt_err_t * np.sqrt(n - 1))

    # S_zal = np.sum(np.sqrt(z - (a_reg + b * t), 2) / (n - 2))
    # S_a = S_zal * np.sqrt((1 / n) + (t ** 2 / sdt_err_t / sdt_err_t / (n - 1)))
    # S_b = S_zal / np.sqrt(sdt_err_t) / np.sqrt(n - 1)

    a_inf = round(a - t_quant * S_a, 4)
    a_sup = round(a + t_quant * S_a, 4)
    b_inf = round(b - t_quant * S_b, 4)
    b_sup = round(b + t_quant * S_b, 4)

    return a_inf, a, a_sup, b_inf, b, b_sup