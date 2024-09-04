import numpy as np
from matplotlib.backends._backend_tk import NavigationToolbar2Tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from paramFuncs import *
from icecream import ic


def funcReversMatr(arr, n):
    A = [[0 for i in range(n)] for j in range(n)]

    for i in range(n):
        for j in range(n):
            A[i][j] = arr[i][j]

    arr_extended = [[0 for i in range(n + n)] for j in range(n)]
    for i in range(n):
        for j in range(n):
            arr_extended[i][j] = arr[i][j]
    for i in range(n):
        arr_extended[i][i + n] = 1

    d = arr_extended[0][0]
    for i in range(n + n):
        if d == 0.0:
            sys.exit('Divide by zero detected!1')
        arr_extended[0][i] /= d

    for i in range(n):

        if arr_extended[i][i] == 0.0:
            return 0

        for j in range(i + 1, n):
            ratio = arr_extended[j][i] / arr_extended[i][i]

            for k in range(n + n):
                arr_extended[j][k] = arr_extended[j][k] - ratio * arr_extended[i][k]

            d = arr_extended[i][i]
            for q in range(n + n):
                if d == 0.0:
                    sys.exit('Divide by zero detected!1')
                arr_extended[i][q] /= d

    d = arr_extended[n - 1][n - 1]
    for i in range(n + n):
        if d == 0.0:
            sys.exit('Divide by zero detected!1')
        arr_extended[n - 1][i] /= d

    for i in range(n - 1, -1, -1):

        if arr_extended[i][i] == 0.0:
            return 0

        for j in range(i - 1, -1, -1):
            ratio = arr_extended[j][i] / arr_extended[i][i]

            for k in range(n + n):
                arr_extended[j][k] = arr_extended[j][k] - ratio * arr_extended[i][k]

            d = arr_extended[i][i]
            for q in range(n + n):
                if d == 0.0:
                    sys.exit('Divide by zero detected!1')
                arr_extended[i][q] /= d

    A_rvrs = [[0 for i in range(n)] for j in range(n)]

    for i in range(n):
        for j in range(n):
            A_rvrs[i][j] = round(arr_extended[i][j + n], 4)

    dot = np.dot(A, A_rvrs)

    return A_rvrs


def erf(val):
    # a1 = 0.0705230784
    # a2 = 0.0422820123
    # a3 = 0.0092705272
    # a4 = 0.0001520143
    # a5 = 0.0002765672
    # a6 = 0.0000430638

    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911

    sign = np.sign(val)
    x = abs(val)

    t = 1.0 / (1.0 + p * x)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * np.exp(-x * x)

    # t = 1 / (1 + a1 * val + a2 * (val ** 2) + a3 * (val ** 3) + a4 * (val ** 4) + a5 * (val ** 5) + a6 * (val ** 6))
    # res = 1.0 - t * np.exp(-x * x)

    # return sign * np.sqrt(res)
    return sign * y


def f1(k):
    res = k ** 2 - (0.5 * (1 - ((-1) ** k)))
    return res


def f2(k):
    res = 5 * (k ** 2) + 22 - (7.5 * (1 - ((-1) ** k)))
    return res


def K_zet(zet, n):
    res = 0
    for i in range(1, 500):
        res += ((-1) ** i) * np.exp(-2 * (i ** 2) * (zet ** 2)) * (
                1 - ((2 * i ** 2 * zet) / (3 * np.sqrt(n))) - (1 / (18 * n)) * (
            (f1(i) - 4 * (f1(i) + 3) * i ** 2 * zet ** 2 + 8 * i ** 4 * zet ** 4)) + (
                        (i ** 2 * zet) / (27 * np.sqrt(n ** 3))) * (
                        (f2(i) ** 2 / 5) - ((4 * (f2(i) + 45) * i ** 2 * zet ** 2) / 15) + 8 * i ** 4 * zet ** 4))

    res = 1 + 2 * res
    return res


def exp_distr(l, x):
    return 1 - np.exp(-l * x)


def exp_up(l, x):
    n = len(x)
    return 1 - np.exp(-l * x) + np.sqrt((l * np.exp(-l * x)) ** 2 * ((x ** 2 * np.exp(-2 * l * x) * l ** 2) / n)) * 1.96


def exp_low(l, x):
    n = len(x)
    return 1 - np.exp(-l * x) - np.sqrt((l * np.exp(-l * x)) ** 2 * ((x ** 2 * np.exp(-2 * l * x) * l ** 2) / n)) * 1.96


def norm_distr(m, sq, x):
    return 0.5 * (1 + erf(((x - m) / (np.sqrt(2) * sq))))


def weib_distr(alf, beta, x):
    return 1 - np.exp(-(x ** beta) / alf)


def uni_distr(a, b, x):
    return (x - a) / (b - a)


def gammaFunc(x):
    if x < 0.5:
        return (np.pi) / (np.sin(np.pi * x) * gammaFunc(1 - x))

    res = 1
    while x > 1.5:
        x -= 1
        res *= x
    return res * np.sqrt(2 * np.pi) * np.exp(-x + 0.5 * np.log(x) + (1 / (12 * x + (1 / (10 * x)))))


def stat_mom(arr, k):
    res = 0
    for i in range(len(arr)):
        res += arr[i] ** k
    res = res / len(arr)
    return res


def centre_mom(arr, k):
    st_m = stat_mom(arr, 1)

    res = 0
    for i in range(len(arr)):
        res += (arr[i] - st_m) ** k

    res = res / len(arr)
    return res


def sort_key(item):
    return item['data']


def empr_func(data, t):
    sorted_data = np.sort(data)

    count_x = np.sum(sorted_data <= t)

    res = count_x / len(sorted_data)

    return res


def removeAnomalous(x, y):
    n = len(x)

    arr_x = []
    arr_y = []
    avr_x = average(x)
    avr_y = average(y)
    stdErr_x = std_err(x, avr_x)
    stdErr_y = std_err(y, avr_y)
    exCf_x = excessCoef(x, avr_x)
    exCf_y = excessCoef(y, avr_y)
    cntrExCf_x = contrExcessCoef(exCf_x)
    cntrExCf_y = contrExcessCoef(exCf_y)

    t_x = 1.2 + 3.6 * (1 - cntrExCf_x) * np.log10(n / 10)
    t_y = 1.2 + 3.6 * (1 - cntrExCf_y) * np.log10(n / 10)

    a_x = avr_x - t_x * stdErr_x
    b_x = avr_x + t_x * stdErr_x
    a_y = avr_y - t_y * stdErr_y
    b_y = avr_y + t_y * stdErr_y

    for i in range(n):
        if x[i] < a_x or x[i] > b_x:
            continue
        arr_x.append(x[i])
        arr_y.append(y[i])

    for i in range(n):
        if y[i] < a_y or y[i] > b_y:
            continue
        arr_y.append(y[i])
        arr_y.append(x[i])

    return arr_x, arr_y


def logs(x, y):
    # x = np.sort(x)
    # y = np.sort(y)
    log_x = []
    log_y = []

    for i in range(len(x)):
        if x[i] < 0:
            x_val = round((np.log10(x[i] + abs(x[0]) + 0.01)), 4)
            log_x.append(x_val)
        else:
            x_val = round(np.log10(x[i]), 5)
            log_x.append(x_val)

    for i in range(len(y)):
        if y[i] < 0:
            y_val = round((np.log10(y[i] + abs(y[0]) + 0.01)), 4)
            log_y.append(y_val)
        else:
            y_val = round(np.log10(y[i]), 5)
            log_y.append(y_val)

    return log_x, log_y


def fi(x):
    return x


def psi(x):
    return np.log(x)


def omega(x):
    return x ** 2


def fi_avr(x, y):
    n = len(x)
    numerator = sum([(fi(x[i]) * omega(y[i])) for i in range(n)])
    denominator = sum([omega(y[i]) for i in range(n)])

    return numerator / denominator


def psi_avr(y):
    n = len(y)
    numerator = sum([(psi(y[i]) * omega(y[i])) for i in range(n)])
    denominator = sum([omega(y[i]) for i in range(n)])

    return numerator / denominator


def fi_avr_sq(x, y):
    n = len(x)
    numerator = sum([(fi(x[i]) ** 2 * omega(y[i])) for i in range(n)])
    denominator = sum([omega(y[i]) for i in range(n)])

    return numerator / denominator


def fi_psi_avrProd(x, y):
    n = len(x)
    numerator = sum([(fi(x[i]) * psi(y[i]) * omega(y[i])) for i in range(n)])
    denominator = sum([omega(y[i]) for i in range(n)])

    return numerator / denominator


def logarimization(sample_data):
    s_n = []
    for i in range(1, len(sample_data) + 1):
        str = f"Вибірка {i}"
        if sample_data[str]['var'].get() == 1:
            s_n.append(str)

    buff = np.array([sample_data[s_n[i]]["data"] for i in range(len(s_n))])
    log_buff = []
    for i in range(len(buff)):
        log_arr = []
        for j in range(len(buff[0])):
            if buff[i][0] < 0:
                x = round((np.log10(buff[i][j] + abs(buff[i][j]) + 0.01)), 4)
                log_arr.append(x)
            else:
                x = round(np.log10(buff[i][j]), 5)
                log_arr.append(x)
        log_buff.append(log_arr)

    for i in range(len(s_n)):
        sample_data[s_n[i]]["data"] = log_buff[i]


def standartization(sample_data):
    s_n = []
    for i in range(1, len(sample_data) + 1):
        str = f"Вибірка {i}"
        if sample_data[str]['var'].get() == 1:
            s_n.append(str)

    buff = np.array([sample_data[s_n[i]]["data"] for i in range(len(s_n))])
    stn_buff = []
    for i in range(len(buff)):
        stn_arr = []
        avr = average(buff[i])
        avrSq = std_err(buff[i], avr)
        for j in range(len(buff[0])):
            x = round(((buff[i][j] - avr) / avrSq), 4)
            stn_arr.append(x)

        stn_buff.append(stn_arr)

    for i in range(len(s_n)):
        sample_data[s_n[i]]["data"] = stn_buff[i]




def removeAnomals(sample_data):
    s_n = []
    for i in range(1, len(sample_data) + 1):
        str = f"Вибірка {i}"
        if sample_data[str]['var'].get() == 1:
            s_n.append(str)

    buff = np.array([sample_data[s_n[i]]["data"] for i in range(len(s_n))])

    data = np.transpose(buff)
    threshold = 1.5
    q1 = np.percentile(data, 25, axis=0)
    q3 = np.percentile(data, 75, axis=0)
    iqr = q3 - q1
    lower_bound = q1 - threshold * iqr
    upper_bound = q3 + threshold * iqr
    mask = np.all((data >= lower_bound) & (data <= upper_bound), axis=1)

    cleaned_samples_t = data[mask]
    cleaned_samples = np.transpose(cleaned_samples_t)
    ic(len(cleaned_samples_t))
    for i in range(len(s_n)):
        sample_data[s_n[i]]["data"] = cleaned_samples[i]
