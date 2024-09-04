import numpy as np

from correlFuncs import *


def linRegressionMNK(x, y, fig1, ax1):
    x_avr = average(x)
    y_avr = average(y)
    std_err_x = std_err(x, x_avr)
    std_err_y = std_err(y, y_avr)
    correl = linearCorrelation(x, y)

    b = correl * (std_err_y / std_err_x)
    a = y_avr - b * x_avr

    y_regression = [a + b * x[i] for i in range(len(x))]

    ax1.plot(x, y_regression, c='red')

    fig1.canvas.draw()
    fig1.canvas.flush_events()


def regrCoffMNK(x, y):
    x_avr = average(x)
    y_avr = average(y)
    std_err_x = std_err(x, x_avr)
    std_err_y = std_err(y, y_avr)
    correl = linearCorrelation(x, y)

    b = round(correl * (std_err_y / std_err_x), 4)
    a = round(y_avr - b * x_avr, 4)

    return a, b


def linRegressionTeyla(x, y, fig1, ax1):
    n = len(x)
    b_arr = []
    for i in range(n - 1):
        for j in range(i + 1, n):
            b_arr.append((y[j] - y[i]) / (x[j] - x[i]))
    b = statistics.median(b_arr)
    a_arr = [(y[i] - b * x[i]) for i in range(n)]
    a = statistics.median(a_arr)

    y_regression = [a + b * x[i] for i in range(len(x))]

    ax1.plot(x, y_regression, c='green')

    fig1.canvas.draw()
    fig1.canvas.flush_events()


def regrCoffTeyla(x, y):
    n = len(x)
    b_arr = []
    for i in range(n - 1):
        for j in range(i + 1, n):
            b_arr.append((y[j] - y[i]) / (x[j] - x[i]))
    b = round(statistics.median(b_arr), 4)
    a_arr = [(y[i] - b * x[i]) for i in range(n)]
    a = round(statistics.median(a_arr), 4)
    return a, b


def determinatCoff(x, y):
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
    x_avr = average(x)
    y_avr = average(y)
    std_err_x = std_err(x, x_avr)
    std_err_y = std_err(y, y_avr)
    correl = linearCorrelation(x, y)

    b = correl * (std_err_y / std_err_x)
    a = y_avr - b * x_avr


    S_zal = np.sqrt(sum([(y[i] - a - b * x[i]) ** 2 for i in range(n)]) / (n - 2))


    r = linearCorrelation(x, y)
    R21 = round(r ** 2 * 100, 4)
    R2 = round((1 - S_zal ** 2 / std_err_y ** 2) * 100, 4)

    return R2


def inervalABLinRegr(x, y):
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
    x_avr = average(x)
    y_avr = average(y)
    std_err_x = std_err(x, x_avr)
    std_err_y = std_err(y, y_avr)
    correl = linearCorrelation(x, y)

    b = correl * (std_err_y / std_err_x)
    a = y_avr - b * x_avr

    S_zal = 0
    for i in range(n):
        S_zal += (y[i] - a - b * x[i]) ** 2

    S_zal = np.sqrt(S_zal / (n - 2))

    S_a = S_zal * np.sqrt(1 / n + x_avr ** 2 / (std_err_x ** 2 * (n - 1)))
    S_b = S_zal / (std_err_x * np.sqrt(n - 1))

    a_inf = round(a - t_quant * S_a, 4)
    a_sup = round(a + t_quant * S_a, 4)
    b_inf = round(b - t_quant * S_b, 4)
    b_sup = round(b + t_quant * S_b, 4)

    return a_inf, a_sup, b_inf, b_sup


def toleranceLim(x, y, fig1, ax1):
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
    x_avr = average(x)
    y_avr = average(y)
    std_err_x = std_err(x, x_avr)
    std_err_y = std_err(y, y_avr)
    correl = linearCorrelation(x, y)

    b = correl * (std_err_y / std_err_x)
    a = y_avr - b * x_avr

    y_regression = [a + b * x[i] for i in range(len(x))]



    S_zal = std_err_y * np.sqrt((1 - correl ** 2) * ((n - 1) / (n - 2)))

    y_min = [y_regression[i] - t_quant * S_zal for i in range(n)]
    print('len(y_min)', len(y_min))
    y_max = [y_regression[i] + t_quant * S_zal for i in range(n)]
    print('len(y_max)', len(y_max))


    ax1.plot(x, y_min, c='red')
    ax1.plot(x, y_max, c='red')

    fig1.canvas.draw()
    fig1.canvas.flush_events()


def predictNewObs(x, y, fig1, ax1):
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

    x_avr = average(x)
    y_avr = average(y)
    std_err_x = std_err(x, x_avr)
    std_err_y = std_err(y, y_avr)
    correl = linearCorrelation(x, y)

    b = correl * (std_err_y / std_err_x)
    a = y_avr - b * x_avr

    S_zal = std_err_y * np.sqrt((1 - correl ** 2) * ((n - 1) / (n - 2)))
    S_b = S_zal / (std_err_x * np.sqrt(n - 1))

    S_y0 = [np.sqrt(S_zal ** 2 * (1 + 1 / n) + S_b ** 2 * (x[i] - x_avr) ** 2) for i in range(n)]
    mx = max(x)
    mn = min(x)
    h = (max(x) - min(x)) / n
    x_for_s = [(min(x) + (i) * h) for i in range(n)]
    # S_y0 = [np.sqrt(S_zal ** 2 * (1 + 1 / n) + S_b ** 2 * (x_for_s[i]  - x_avr) ** 2) for i in range(n)]


    y_regression = [a + b * x[i] for i in range(len(x))]

    y_inf = [y_regression[i] - t_quant * S_y0[i] for i in range(n)]

    y_sup = [y_regression[i] + t_quant * S_y0[i] for i in range(n)]




    ax1.plot(np.sort(x), np.sort(y_inf), c='green')
    ax1.plot(np.sort(x), np.sort(y_sup), c='green')

    fig1.canvas.draw()
    fig1.canvas.flush_events()




def confInterLinRegr(x, y, fig1, ax1):
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

    x_avr = average(x)
    y_avr = average(y)
    std_err_x = std_err(x, x_avr)
    std_err_y = std_err(y, y_avr)
    correl = linearCorrelation(x, y)

    b = correl * (std_err_y / std_err_x)
    a = y_avr - b * x_avr

    S_zal = std_err_y * np.sqrt((1 - correl ** 2) * ((n - 1) / (n - 2)))
    S_b = S_zal / (std_err_x * np.sqrt(n - 1))

    h = (max(x) - min(x)) / n
    x_for_s = [(min(x) + (i) * h) for i in range(n)]


    S_y0 = [np.sqrt(S_zal ** 2 * (1 / n) + S_b ** 2 * (x[i] - x_avr) ** 2) for i in range(n)]
    # S_y0 = [np.sqrt(S_zal ** 2 * (1 / n) + S_b ** 2 * (x_for_s[i] - x_avr) ** 2) for i in range(n)]
    y_regression = [a + b * x[i] for i in range(len(x))]
    y_inf = [(y_regression[i] - t_quant * S_y0[i]) for i in range(n)]
    y_sup = [(y_regression[i] + t_quant * S_y0[i]) for i in range(n)]
    ax1.plot(np.sort(x), np.sort(y_inf), c='green')
    ax1.plot(np.sort(x), np.sort(y_sup), c='green')
    # ax1.scatter(x, y_inf, c='green', s=1)
    # ax1.scatter(x, y_sup, c='green', s=1)

    fig1.canvas.draw()
    fig1.canvas.flush_events()

def adequateReproductionRegr(x, y):
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

    x_avr = average(x)
    y_avr = average(y)
    std_err_x = std_err(x, x_avr)
    std_err_y = std_err(y, y_avr)
    correl = linearCorrelation(x, y)

    b = correl * (std_err_y / std_err_x)
    a = y_avr - b * x_avr

    S_zal = 0
    for i in range(n):
        S_zal += (y[i] - a - b * x[i]) ** 2

    S_zal = np.sqrt(S_zal / (n - 2))

    f = round((S_zal ** 2) / (std_err_y ** 2), 4)
    # f = correl ** 2

    return f