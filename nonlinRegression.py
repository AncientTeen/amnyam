from correlFuncs import *


def parabolicRegression(x, y, fig1, ax1):
    x_avr = average(x)
    x_avr_sq = average_sq(x)
    x_avr_cub = average_cub(x)
    x_avr_frth = average_frth(x)
    y_avr = average(y)
    avr_abrakadabra = average_abrakadabra(x, y)
    std_err_x = std_err(x, x_avr)
    std_err_y = std_err(y, y_avr)
    correl = linearCorrelation(x, y)

    b = ((x_avr_frth - x_avr_sq ** 2) * correl * std_err_x * std_err_y - (
            x_avr_cub - x_avr_sq * x_avr) * avr_abrakadabra) / (
                std_err_x ** 2 * (x_avr_frth - x_avr_sq ** 2) - (x_avr_cub - x_avr_sq * x_avr) ** 2)

    c = (std_err_x ** 2 * avr_abrakadabra - (x_avr_cub - x_avr_sq * x_avr) * correl * std_err_x * std_err_y) / (
            std_err_x ** 2 * (x_avr_frth - x_avr_sq ** 2) - (x_avr_cub - x_avr_sq * x_avr) ** 2)

    a = y_avr - b * x_avr - c * x_avr_sq

    y_regression = [(a + b * x[i] + c * x[i] ** 2) for i in range(len(x))]

    ax1.plot(x, y_regression, c='red')

    fig1.canvas.draw()
    fig1.canvas.flush_events()


def inervalNonLinRegr(x, y):
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

    n = len(x)
    x_avr = average(x)
    x_avr_sq = average_sq(x)
    x_avr_cub = average_cub(x)
    x_avr_frth = average_frth(x)
    y_avr = average(y)
    avr_abrakadabra = average_abrakadabra(x, y)
    avr_abrakadabra1 = average_abrakadabra1(x, y)
    std_err_x = std_err(x, x_avr)
    std_err_y = std_err(y, y_avr)
    correl = linearCorrelation(x, y)

    fi1 = [x[i] - x_avr for i in range(n)]
    fi2 = [(x[i] ** 2 - (x_avr_cub - x_avr_sq * x_avr) / std_err_x - x_avr_sq) for i in range(n)]

    a1 = y_avr

    b1 = round(avr_abrakadabra1 / (std_err_x ** 2), 4)

    fi2_avr_sq = average_sq(fi2)
    fi2_y_avr_prod = averageProduct(fi2, y)

    c1 = round(fi2_y_avr_prod / fi2_avr_sq, 4)

    S_zal = 0
    for i in range(n):
        S_zal += (y[i] - a1 - b1 * fi1[i] - c1 * fi2[i]) ** 2

    S_zal = S_zal / (n - 3)
    S_a = S_zal / np.sqrt(n)
    S_b = S_zal / (std_err_x * np.sqrt(n))
    S_c = S_zal / np.sqrt(n * fi2_avr_sq)

    a_inf = round(a1 - t_quant * S_a, 4)
    a_sup = round(a1 + t_quant * S_a, 4)
    b_sup = round(b1 + t_quant * S_b, 4)
    b_inf = round(b1 - t_quant * S_b, 4)
    c_sup = round(c1 + t_quant * S_c, 4)
    c_inf = round(c1 - t_quant * S_c, 4)

    # return a, a1, b, b1, c, c1
    return a_inf, a1, a_sup, b_inf, b1, b_sup, c_inf, c1, c_sup


def toleranceLimNonLin(x, y, fig1, ax1):
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
    x_avr_sq = average_sq(x)
    x_avr_cub = average_cub(x)
    x_avr_frth = average_frth(x)
    y_avr = average(y)
    avr_abrakadabra = average_abrakadabra(x, y)
    avr_abrakadabra1 = average_abrakadabra1(x, y)
    std_err_x = std_err(x, x_avr)
    std_err_y = std_err(y, y_avr)
    correl = linearCorrelation(x, y)

    fi1 = [x[i] - x_avr for i in range(n)]
    fi2 = [(x[i] ** 2 - (x_avr_cub - x_avr_sq * x_avr) / std_err_x - x_avr_sq) for i in range(n)]

    a1 = y_avr

    b1 = avr_abrakadabra1 / (std_err_x ** 2)

    fi2_avr_sq = average_sq(fi2)
    fi2_y_avr_prod = averageProduct(fi2, y)

    c1 = fi2_y_avr_prod / fi2_avr_sq

    S_zal = 0
    for i in range(n):
        S_zal += (y[i] - a1 - b1 * fi1[i] - c1 * fi2[i]) ** 2

    S_zal = S_zal / (n - 3)

    y_min = [(a1 + b1 * fi1[i] + c1 * fi2[i] - t_quant * S_zal) for i in range(len(x))]
    y_max = [(a1 + b1 * fi1[i] + c1 * fi2[i] + t_quant * S_zal) for i in range(len(x))]

    ax1.plot(x, y_min, c='red')
    ax1.plot(x, y_max, c='red')

    fig1.canvas.draw()
    fig1.canvas.flush_events()


def confInterNonLinRegr(x, y, fig1, ax1):
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
    x_avr_sq = average_sq(x)
    x_avr_cub = average_cub(x)
    x_avr_frth = average_frth(x)
    y_avr = average(y)
    avr_abrakadabra = average_abrakadabra(x, y)
    avr_abrakadabra1 = average_abrakadabra1(x, y)
    std_err_x = std_err(x, x_avr)
    std_err_y = std_err(y, y_avr)
    correl = linearCorrelation(x, y)

    b = round(((x_avr_frth - x_avr_sq ** 2) * correl * std_err_x * std_err_y - (
            x_avr_cub - x_avr_sq * x_avr) * avr_abrakadabra) / (
                      std_err_x ** 2 * (x_avr_frth - x_avr_sq ** 2) - (x_avr_cub - x_avr_sq * x_avr) ** 2), 4)

    c = round((std_err_x ** 2 * avr_abrakadabra - (x_avr_cub - x_avr_sq * x_avr) * correl * std_err_x * std_err_y) / (
            std_err_x ** 2 * (x_avr_frth - x_avr_sq ** 2) - (x_avr_cub - x_avr_sq * x_avr) ** 2), 4)

    a = round(y_avr - b * x_avr - c * x_avr_sq, 4)

    fi1 = [x[i] - x_avr for i in range(n)]
    fi2 = [(x[i] ** 2 - (x_avr_cub - x_avr_sq * x_avr) / std_err_x - x_avr_sq) for i in range(n)]

    a1 = y_avr

    b1 = avr_abrakadabra1 / (std_err_x ** 2)

    fi2_avr_sq = average_sq(fi2)
    fi2_y_avr_prod = averageProduct(fi2, y)

    c1 = fi2_y_avr_prod / fi2_avr_sq

    S_zal = 0
    for i in range(n):
        S_zal += (y[i] - a1 - b1 * fi1[i] - c1 * fi2[i]) ** 2

    S_zal = S_zal / (n - 3)
    S_a = S_zal / np.sqrt(n)
    S_b = S_zal / (std_err_x * np.sqrt(n))
    S_c = S_zal / np.sqrt(n * fi2_avr_sq)

    S_y = [(np.sqrt((1 / n) * S_zal ** 2 + S_b ** 2 * fi1[i] ** 2 + S_c ** 2 * fi2[i] ** 2)) for i in range(n)]

    y_regression = [(a + b * x[i] + c * x[i] ** 2) for i in range(n)]

    y_low = [(y_regression[i] - t_quant * S_y[i]) for i in range(len(x))]
    y_high = [(y_regression[i] + t_quant * S_y[i]) for i in range(len(x))]

    ax1.plot(x, y_low, c='green')
    ax1.plot(x, y_high, c='green')

    fig1.canvas.draw()
    fig1.canvas.flush_events()


def predictNewObsNonLin(x, y, fig1, ax1):
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
    x_avr_sq = average_sq(x)
    x_avr_cub = average_cub(x)
    x_avr_frth = average_frth(x)
    y_avr = average(y)
    avr_abrakadabra = average_abrakadabra(x, y)
    avr_abrakadabra1 = average_abrakadabra1(x, y)
    std_err_x = std_err(x, x_avr)
    std_err_y = std_err(y, y_avr)
    correl = linearCorrelation(x, y)

    b = round(((x_avr_frth - x_avr_sq ** 2) * correl * std_err_x * std_err_y - (
            x_avr_cub - x_avr_sq * x_avr) * avr_abrakadabra) / (
                      std_err_x ** 2 * (x_avr_frth - x_avr_sq ** 2) - (x_avr_cub - x_avr_sq * x_avr) ** 2), 4)

    c = round((std_err_x ** 2 * avr_abrakadabra - (x_avr_cub - x_avr_sq * x_avr) * correl * std_err_x * std_err_y) / (
            std_err_x ** 2 * (x_avr_frth - x_avr_sq ** 2) - (x_avr_cub - x_avr_sq * x_avr) ** 2), 4)

    a = round(y_avr - b * x_avr - c * x_avr_sq, 4)

    fi1 = [x[i] - x_avr for i in range(n)]
    fi2 = [(x[i] ** 2 - (x_avr_cub - x_avr_sq * x_avr) / std_err_x - x_avr_sq) for i in range(n)]

    a1 = y_avr

    b1 = avr_abrakadabra1 / (std_err_x ** 2)

    fi2_avr_sq = average_sq(fi2)
    fi2_y_avr_prod = averageProduct(fi2, y)

    c1 = fi2_y_avr_prod / fi2_avr_sq

    S_zal = 0
    for i in range(n):
        S_zal += (y[i] - a1 - b1 * fi1[i] - c1 * fi2[i]) ** 2

    S_zal = S_zal / (n - 3)

    S_zal1 = std_err_y * np.sqrt((1 - correl ** 2) * ((n - 1) / (n - 2)))

    S_a = S_zal / np.sqrt(n)
    S_b = S_zal / (std_err_x * np.sqrt(n))
    S_c = S_zal / np.sqrt(n * fi2_avr_sq)

    # S_y = [(np.sqrt((1 + 1 / n) * S_zal1 ** 2 + S_b ** 2 * fi1[i] ** 2 + S_c ** 2 * fi2[i] ** 2)) for i in range(n)]
    S_y = [((S_zal / np.sqrt(n)) * np.sqrt(n + 1 + fi1[i] ** 2 / std_err_x ** 2 + fi2[i] ** 2 / fi2_avr_sq)) for i in
           range(n)]

    y_regression = [(a + b * x[i] + c * x[i] ** 2) for i in range(n)]

    y_low = [(y_regression[i] - t_quant * S_y[i]) for i in range(len(x))]
    y_high = [(y_regression[i] + t_quant * S_y[i]) for i in range(len(x))]

    ax1.plot(x, y_low, c='black')
    ax1.plot(x, y_high, c='black')

    fig1.canvas.draw()
    fig1.canvas.flush_events()


def determinatCoffNonLinRegr(x, y):
    n = len(x)
    x_avr = average(x)
    x_avr_sq = average_sq(x)
    x_avr_cub = average_cub(x)
    x_avr_frth = average_frth(x)
    y_avr = average(y)
    avr_abrakadabra = average_abrakadabra(x, y)
    avr_abrakadabra1 = average_abrakadabra1(x, y)
    std_err_x = std_err(x, x_avr)
    std_err_y = std_err(y, y_avr)
    correl = linearCorrelation(x, y)

    fi1 = [x[i] - x_avr for i in range(n)]
    fi2 = [(x[i] ** 2 - (x_avr_cub - x_avr_sq * x_avr) / std_err_x - x_avr_sq) for i in range(n)]

    a1 = y_avr

    b1 = avr_abrakadabra1 / (std_err_x ** 2)

    fi2_avr_sq = average_sq(fi2)
    fi2_y_avr_prod = averageProduct(fi2, y)

    c1 = fi2_y_avr_prod / fi2_avr_sq

    S_zal = 0
    for i in range(n):
        S_zal += (y[i] - a1 - b1 * fi1[i] - c1 * fi2[i]) ** 2

    S_zal = S_zal / (n - 3)

    R2 = abs(round((1 - S_zal ** 2 / std_err_y ** 2) * 100, 4))

    return R2


def tValsForParams(x, y):
    n = len(x)
    x_avr = average(x)
    x_avr_sq = average_sq(x)
    x_avr_cub = average_cub(x)
    y_avr = average(y)
    avr_abrakadabra1 = average_abrakadabra1(x, y)
    std_err_x = std_err(x, x_avr)
    std_err_y = std_err(y, y_avr)

    fi1 = [x[i] - x_avr for i in range(n)]
    fi2 = [(x[i] ** 2 - (x_avr_cub - x_avr_sq * x_avr) / std_err_x - x_avr_sq) for i in range(n)]

    a1 = y_avr

    b1 = avr_abrakadabra1 / (std_err_x ** 2)

    fi2_avr_sq = average_sq(fi2)
    fi2_y_avr_prod = averageProduct(fi2, y)

    c1 = fi2_y_avr_prod / fi2_avr_sq

    S_zal = 0
    for i in range(n):
        S_zal += (y[i] - a1 - b1 * fi1[i] - c1 * fi2[i]) ** 2

    S_zal = S_zal / (n - 3)

    t_a = abs(round((a1 / S_zal) * np.sqrt(n), 4))
    t_b = abs(round((b1 * std_err_x / S_zal) * np.sqrt(n), 4))
    t_c = abs(round((c1 / S_zal) * np.sqrt(n * fi2_avr_sq), 4))

    return t_a, t_b, t_c
