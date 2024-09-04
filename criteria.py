from tkinter import *
from scipy.stats import chi2
from funcs import *
from tkinter import *

from scipy.stats import chi2

from funcs import *


def avr_match_dep(sample_data, T3):
    s_n_buff = []
    s_check = 0
    for i in range(1, len(sample_data) + 1):
        s_n = f"Вибірка {i}"
        state = sample_data[s_n]["var"].get()
        if state == 1:
            s_n_buff.append(s_n)
            s_check += 1
    if s_check == 2:
        arr_buff = [sample_data[s_n_buff[0]]["data"], sample_data[s_n_buff[1]]["data"]]
        if len(arr_buff[0]) != len(arr_buff[1]):
            T3.delete('1.0', END)
            T3.insert(END, 'Вибірки мають різні розміри\n')

        # arr_buff[0] = shellSort(arr_buff[0], len(arr_buff[0]))
        # arr_buff[1] = shellSort(arr_buff[1], len(arr_buff[1]))

        degrees_of_freedom = len(arr_buff[0]) - 2
        quantile = 0
        if degrees_of_freedom > 18 and degrees_of_freedom <= 25:
            quantile = 2.08
        elif degrees_of_freedom > 25 and degrees_of_freedom <= 60:
            quantile = 2.03
        elif degrees_of_freedom > 60:
            quantile = 1.96

        z = []
        for i in range(len(arr_buff[0])):
            z.append(arr_buff[0][i] - arr_buff[1][i])

        z = np.sort(z)
        avr_z = average(z)
        avrSq_z = std_err(z, avr_z)
        t = (avr_z * (len(arr_buff[0]) ** (1 / 2))) / avrSq_z

        T3.delete('1.0', END)
        T3.insert(END, 'Перевірка збігу середніх\n')
        if abs(t) < quantile:
            T3.insert(END, f"t = {round(abs(t), 4)} < {quantile}, отже, середні збігаються\n")
        else:
            T3.insert(END, f"t = {round(abs(t), 4)} > {quantile}, отже, середні не збігаються\n")

        # T3.insert(END, f"lib t-test dep = {stats.ttest_rel(arr_buff[0], arr_buff[1])}\n")



    else:
        T3.delete('1.0', END)
        T3.insert(END, 'Більше 2-ох вибірок\n')

def sign_test(sample_data, T3):
    s_n_buff = []
    s_check = 0
    for i in range(1, len(sample_data) + 1):
        s_n = f"Вибірка {i}"
        state = sample_data[s_n]["var"].get()
        if state == 1:
            s_n_buff.append(s_n)
            s_check += 1
    if s_check == 2:
        arr_buff = [sample_data[s_n_buff[0]]["data"], sample_data[s_n_buff[1]]["data"]]

        if len(arr_buff[0]) != len(arr_buff[1]):
            T3.delete('1.0', END)
            T3.insert(END, 'Вибірки мають різні розміри\n')

        z = []
        for i in range(len(arr_buff[0])):
            z.append(arr_buff[0][i] - arr_buff[1][i])

        u = []
        for i in range(len(z)):
            if z[i] > 0:
                u.append(1)

        n = len(u)
        S = sum(u)
        alpha0 = 0
        S_star = 0

        if n <= 15:
            for l in range(n - S + 1):
                numerator = np.math.factorial(l)
                denominator = np.math.factorial(n) * np.math.factorial(n - l)
                alpha0 += numerator / denominator
            alpha0 *= 2 ** (-n)
        else:
            S_star = (2 * S - 1 - n) / np.sqrt(n)

        T3.delete('1.0', END)
        T3.insert(END, 'Критерій знаків:\n')
        if n <= 15:
            if abs(alpha0) < 1.64:
                T3.insert(END, f"v = {round(abs(alpha0), 4)} < 1.64, отже, вибірки однорідні\n")
            else:
                T3.insert(END, f"v = {round(abs(alpha0), 4)} > 1.64, отже, вибірки неоднорідні\n")
        else:
            if abs(S_star) < 1.64:
                T3.insert(END, f"v = {round(abs(S_star), 4)} < 1.64, отже, вибірки однорідні\n")
            else:
                T3.insert(END, f"v = {round(abs(S_star), 4)} > 1.64, отже, вибірки неоднорідні\n")

    else:
        T3.delete('1.0', END)
        T3.insert(END, 'Більше 2-ох вибірок\n')


def q_test(sample_data, T3):
    s_n_buff = []
    s_check = 0
    for i in range(1, len(sample_data) + 1):
        s_n = f"Вибірка {i}"
        state = sample_data[s_n]["var"].get()
        if state == 1:
            s_n_buff.append(s_n)
            s_check += 1

    arr_buff = []
    for i in range(len(s_n_buff)):
        arr_buff.append(sample_data[s_n_buff[i]]["data"])

    k = len(arr_buff)

    degrees_of_freedom = k - 1
    quantile = 0
    if degrees_of_freedom == 1:
        quantile = 0.0039
    elif degrees_of_freedom == 2:
        quantile = 0.103
    elif degrees_of_freedom == 3:
        quantile = 0.352
    elif degrees_of_freedom == 4:
        quantile = 0.711
    elif degrees_of_freedom == 5:
        quantile = 1.15
    elif degrees_of_freedom == 6:
        quantile = 1.64

    # sort_arr_buff = []
    # for i in range(k):
    #     sort_arr_buff.append(np.sort(arr_buff[i]))

    avr_buff = []
    for i in range(k):
        avr_buff.append(average(arr_buff[i]))

    bin_arr_buff = []
    for i in range(k):
        bin_buff = []
        for j in range(len(arr_buff[i])):
            if arr_buff[i][j] - avr_buff[i] > 0:
                bin_buff.append(1)
            else:
                bin_buff.append(0)

        bin_arr_buff.append(bin_buff)

    t_buff = []
    t_mean = 0
    for i in range(k):
        t = 0
        for j in range(len(bin_arr_buff[i])):
            if bin_arr_buff[i][j] == 1:
                t += 1
        t_mean += t
        t_buff.append(t)
    t_mean = t_mean / k

    u_sum = sum(t_buff)
    u_sq_sum = 0
    for i in range(len(arr_buff[0])):
        u_sq = 0
        for j in range(k):
            if bin_arr_buff[j][i] == 1:
                u_sq += 1
        u_sq_sum += u_sq ** 2

    t_disp = 0
    for i in range(k):
        t_disp += (t_buff[i] - t_mean) ** 2

    Q = (k * (k - 1) * t_disp) / (k * u_sum - u_sq_sum)

    T3.delete('1.0', END)
    T3.insert(END, 'Q критерій\n')
    # T3.insert(END, f"Q value = {round(Q, 4)}\n")
    if Q <= quantile:
        T3.insert(END, f"Q value = {round(Q, 4)} <= {quantile}, отже, вибірки однорідні\n")
    else:
        T3.insert(END, f"Q value = {round(Q, 4)} > {quantile}, отже, вибірки неоднорідні\n")


def abbe_test(sample_data, T3):
    s_n_buff = []
    s_check = 0
    for i in range(1, len(sample_data) + 1):
        s_n = f"Вибірка {i}"
        state = sample_data[s_n]["var"].get()
        if state == 1:
            s_n_buff.append(s_n)
            s_check += 1
    if s_check == 2:
        arr_buff = [sample_data[s_n_buff[0]]["data"], sample_data[s_n_buff[1]]["data"]]

        if len(arr_buff[0]) != len(arr_buff[1]):
            T3.delete('1.0', END)
            T3.insert(END, 'Вибірки мають різні розміри\n')

        # arr_buff[0] = shellSort(arr_buff[0], len(arr_buff[0]))
        # arr_buff[1] = shellSort(arr_buff[1], len(arr_buff[1]))

        n1 = len(arr_buff[0])
        n2 = len(arr_buff[1])

        z = []
        for i in range(len(arr_buff[0])):
            z.append(arr_buff[0][i] - arr_buff[1][i])

        d_sq = 0
        for i in range(len(z) - 1):
            d_sq += (z[i + 1] - z[i]) ** 2
        d_sq = d_sq / (len(z) - 1)

        z = np.sort(z)

        avr_z = average(z)
        disp_z = std_err(z, avr_z) ** 2

        q = d_sq / (2 * disp_z)

        U = (q - 1) * np.sqrt((len(z) ** 2 - 1) / (len(z) - 2))

        df = 1
        q = 0.95
        n = n1 + n2
        chi2_critical = chi2.ppf(q, df)
        abbe_critical = round((df * (n - 1)) / (n - 1) * (1 - chi2_critical / (df * (n - 1))), 4)

        T3.delete('1.0', END)
        T3.insert(END, 'Критерій Аббе:\n')
        if abs(U) < abbe_critical:
            T3.insert(END, f"U value = {round(abs(U), 4)} < {abbe_critical}, отже, вибірки однорідні\n")
        else:
            T3.insert(END, f"U value = {round(abs(U), 4)} > {abbe_critical}, отже, вибірки неоднорідні\n")


    else:
        T3.delete('1.0', END)
        T3.insert(END, 'Більше 2-ох вибірок\n')


def avr_match_indep(sample_data, T3):
    s_n_buff = []
    s_check = 0
    for i in range(1, len(sample_data) + 1):
        s_n = f"Вибірка {i}"
        state = sample_data[s_n]["var"].get()
        if state == 1:
            s_n_buff.append(s_n)
            s_check += 1
    if s_check == 2:
        arr_buff = [sample_data[s_n_buff[0]]["data"], sample_data[s_n_buff[1]]["data"]]
        arr_buff[0] = np.sort(arr_buff[0])
        arr_buff[1] = np.sort(arr_buff[1])

        x_avr = average(arr_buff[0])
        y_avr = average(arr_buff[1])
        n1 = len(arr_buff[0])
        n2 = len(arr_buff[1])

        degrees_of_freedom = n1 - 2
        quantile = 0
        if degrees_of_freedom > 18 and degrees_of_freedom <= 25:
            quantile = 2.08
        elif degrees_of_freedom > 25 and degrees_of_freedom <= 60:
            quantile = 2.03
        elif degrees_of_freedom > 60:
            quantile = 1.96

        xSq = std_err(arr_buff[0], x_avr)
        ySq = std_err(arr_buff[1], y_avr)

        if len(arr_buff[0]) + len(arr_buff[1]) > 25:

            t = (x_avr - y_avr) / ((xSq ** 2 / n1 + ySq ** 2 / n2) ** (1 / 2))


        else:
            xSq_avr = xSq ** 2 / n1
            ySq_avr = ySq ** 2 / n2

            t = ((x_avr - y_avr) / ((((n1 - 1) * xSq_avr + (n2 - 1) * ySq_avr) / (n1 + n2 - 2)) ** (1 / 2))) * (
                    ((n1 * n2) / (n1 + n2)) ** (1 / 2))

        T3.delete('1.0', END)
        T3.insert(END, 'Перевірка збігу середніх\n')
        if abs(t) < quantile:
            T3.insert(END, f"t = {round(abs(t), 4)} < {quantile}, отже, середні збігаються\n")
        else:
            T3.insert(END, f"t = {round(abs(t), 4)} > {quantile}, отже, середні не збігаються\n")
        # T3.insert(END, f"lib t-test indep = {stats.ttest_ind(arr_buff[0], arr_buff[1])}\n")


    else:
        T3.delete('1.0', END)
        T3.insert(END, 'Більше 2-ох вибірок\n')


def f_test(sample_data, T3):
    s_n_buff = []
    s_check = 0
    for i in range(1, len(sample_data) + 1):
        s_n = f"Вибірка {i}"
        state = sample_data[s_n]["var"].get()
        if state == 1:
            s_n_buff.append(s_n)
            s_check += 1
    if s_check == 2:
        arr_buff = [sample_data[s_n_buff[0]]["data"], sample_data[s_n_buff[1]]["data"]]
        if len(arr_buff[0]) != len(arr_buff[1]):
            T3.delete('1.0', END)
            T3.insert(END, 'Вибірки мають різні розміри\n')
        arr_buff[0] = np.sort(arr_buff[0])
        arr_buff[1] = np.sort(arr_buff[1])

        n1 = len(arr_buff[0])
        n2 = len(arr_buff[1])

        degrees_of_freedom_1 = n1 - 1
        degrees_of_freedom_2 = n2 - 1
        quantile = 0
        if degrees_of_freedom_1 > 20 and degrees_of_freedom_2 > 20 and degrees_of_freedom_1 <= 24 and degrees_of_freedom_2 <= 24:
            quantile = 2.05
        elif degrees_of_freedom_1 > 24 and degrees_of_freedom_2 > 24 and degrees_of_freedom_1 <= 60 and degrees_of_freedom_2 <= 60:
            quantile = 1.75
        elif degrees_of_freedom_1 > 60 and degrees_of_freedom_2 > 60:
            quantile = 1.265

        x_avr = average(arr_buff[0])
        y_avr = average(arr_buff[1])

        s_x = std_err(arr_buff[0], x_avr) ** 2
        s_y = std_err(arr_buff[1], y_avr) ** 2

        f = 0
        if s_x >= s_y:
            f = s_x / s_y
        else:
            f = s_y / s_x

        T3.delete('1.0', END)
        T3.insert(END, 'Перевірка збігу дисперсій\n')
        # T3.insert(END, f"f = {round(f, 4)}\n")
        if f < quantile:
            T3.insert(END, f"f = {round(f, 4)} < {quantile}, отже, дисперсії збігаються\n")
        else:
            T3.insert(END, f"f = {round(f, 4)} > {quantile}, отже, дисперсії не збігаються\n")


    else:
        T3.delete('1.0', END)
        T3.insert(END, 'Більше 2-ох вибірок\n')


def Bartlett_test(sample_data, T3):
    s_n_buff = []
    s_check = 0
    for i in range(1, len(sample_data) + 1):
        s_n = f"Вибірка {i}"
        state = sample_data[s_n]["var"].get()
        if state == 1:
            s_n_buff.append(s_n)
            s_check += 1

    arr_buff = []
    for i in range(len(s_n_buff)):
        arr_buff.append(sample_data[s_n_buff[i]]["data"])

    k = len(arr_buff)
    for i in range(k):
        arr_buff[i] = np.sort(arr_buff[i])

    degrees_of_freedom = k - 1
    quantile = 0
    if degrees_of_freedom == 1:
        quantile = 0.0039
    elif degrees_of_freedom == 2:
        quantile = 0.103
    elif degrees_of_freedom == 3:
        quantile = 0.352
    elif degrees_of_freedom == 4:
        quantile = 0.711
    elif degrees_of_freedom == 5:
        quantile = 1.15
    elif degrees_of_freedom == 6:
        quantile = 1.64

    avr_buff = []
    for i in range(k):
        avr_buff.append(average(arr_buff[i]))

    disp_buff = []
    for i in range(k):
        disp_buff.append(std_err(arr_buff[i], avr_buff[i]) ** 2)

    numerator = 0
    denominator = 0
    for i in range(k):
        numerator += (len(arr_buff[i]) - 1) * disp_buff[i]
        denominator += len(arr_buff[i]) - 1

    S_2 = numerator / denominator

    B = 0
    for i in range(k):
        B += (len(arr_buff[i]) - 1) * np.log(disp_buff[i] / S_2)

    B *= -1

    C = 0
    f = 0
    d_s = 0
    for i in range(k):
        f += 1 / (len(arr_buff[i]) - 1)
        d_s += len(arr_buff[i]) - 1

    C = 1 + (1 / (3 * (k - 1))) * (f - 1 / d_s)

    xi_val = B / C

    T3.delete('1.0', END)
    T3.insert(END, 'Критерій Бартлетта\n')
    # T3.insert(END, f"xi value = {round(xi_val, 4)}\n")
    if xi_val < quantile:
        T3.insert(END,
                  f"xi_val = {round(xi_val, 4)} < {quantile}, отже, дисперсії збігаються\n")
    else:
        T3.insert(END,
                  f"xi_val = {round(xi_val, 4)} > {quantile}, отже, дисперсії  не збігаються\n")

    # T3.insert(END, f"lib t-test dep = {stats.ttest_rel(arr_buff[0], arr_buff[1])}\n")


def anova_func(sample_data, T3):
    s_n_buff = []
    s_check = 0
    for i in range(1, len(sample_data) + 1):
        s_n = f"Вибірка {i}"
        state = sample_data[s_n]["var"].get()
        if state == 1:
            s_n_buff.append(s_n)
            s_check += 1

    arr_buff = []
    for i in range(len(s_n_buff)):
        arr_buff.append(sample_data[s_n_buff[i]]["data"])

    k = len(arr_buff)
    for i in range(k):
        arr_buff[i] = np.sort(arr_buff[i])
    avr_buff = []
    for i in range(k):
        avr_buff.append(average(arr_buff[i]))

    disp_buff = []
    for i in range(k):
        disp_buff.append(std_err(arr_buff[i], avr_buff[i]) ** 2)

    gen_mean = 0
    gen_vol = 0
    for i in range(k):
        gen_mean += len(arr_buff[i]) * avr_buff[i]
        gen_vol += len(arr_buff[i])
    gen_mean = gen_mean / gen_vol

    """intergroup variation"""
    intrgrp_var = 0
    for i in range(k):
        intrgrp_var += len(arr_buff[i]) * ((avr_buff[i] - gen_mean) ** 2)
    intrgrp_var = intrgrp_var / (k - 1)

    """Variation within each sample"""
    var_wthin_each_sample = 0
    for i in range(k):
        var_wthin_each_sample += len(arr_buff[i] - 1) * disp_buff[i]
    var_wthin_each_sample = var_wthin_each_sample / (gen_vol - k)

    f_statistic = intrgrp_var / var_wthin_each_sample

    degrees_of_freedom_1 = k - 1
    degrees_of_freedom_2 = gen_vol - k
    quantile = 0
    if degrees_of_freedom_1 == 1 and degrees_of_freedom_2 > 19 and degrees_of_freedom_2 <= 24:
        quantile = 4.32
    elif degrees_of_freedom_1 == 1 and degrees_of_freedom_2 > 24 and degrees_of_freedom_2 <= 60:
        quantile = 4.13
    elif degrees_of_freedom_1 == 1 and degrees_of_freedom_2 > 60:
        quantile = 3.92
    elif degrees_of_freedom_1 == 2 and degrees_of_freedom_2 > 19 and degrees_of_freedom_2 <= 24:
        quantile = 3.45
    elif degrees_of_freedom_1 == 2 and degrees_of_freedom_2 > 24 and degrees_of_freedom_2 <= 60:
        quantile = 3.27
    elif degrees_of_freedom_1 == 2 and degrees_of_freedom_2 > 60:
        quantile = 3.08

    T3.delete('1.0', END)
    T3.insert(END, 'Однофакторний дисперсійний аналіз\n')
    if f_statistic < quantile:
        T3.insert(END, f"F value = {round(f_statistic, 4)} < {quantile}, отже, середні збігаються\n")
    else:
        T3.insert(END, f"F value = {round(f_statistic, 4)} > {quantile}, отже, середні не збігаються\n")
    # T3.insert(END, f"lib t-test dep = {stats.ttest_rel(arr_buff[0], arr_buff[1])}\n")


"""Kolmogorov_Smirnov_test that really work"""


def Kolmogorov_Smirnov_test(sample_data, T3):
    # Сортуємо вибірки
    s_n_buff = []
    s_check = 0
    for i in range(1, len(sample_data) + 1):
        s_n = f"Вибірка {i}"
        state = sample_data[s_n]["var"].get()
        if state == 1:
            s_n_buff.append(s_n)
            s_check += 1
    if s_check == 2:
        arr_buff = [sample_data[s_n_buff[0]]["data"], sample_data[s_n_buff[1]]["data"]]

        sample1 = np.sort(arr_buff[0])
        sample2 = np.sort(arr_buff[1])

        # Обчислюємо емпіричні функції розподілу для обох вибірок
        n1 = sample1.size
        n2 = sample2.size

        # Об'єднуємо вибірки без дублікатів для порівняння
        all_values = np.unique(np.concatenate([sample1, sample2]))
        # Обчислюємо емпіричні функції розподілу на об'єднаних даних
        cdf1_all = np.searchsorted(sample1, all_values, side='right') / n1
        cdf2_all = np.searchsorted(sample2, all_values, side='right') / n2

        # Знаходимо максимальну абсолютну різницю між емпіричними функціями розподілу
        max_difference = np.max(np.abs(cdf1_all - cdf2_all))

        # Обчислюємо статистику Смірнова
        smirnov_statistic = np.sqrt((n1 * n2) / (n1 + n2)) * max_difference

        T3.delete('1.0', END)
        T3.insert(END, 'Критерій Смирнова-Колмогорова:\n')
        if smirnov_statistic < 0.05:
            T3.insert(END, f"kstets = {round(smirnov_statistic, 4)} < 0.05, отже, вибірки однорідні\n")
        else:
            T3.insert(END, f"kstets = {round(smirnov_statistic, 4)} > 0.05, отже, вибірки неоднорідні\n")

        # T3.insert(END, f"lib kstest = {stats.kstest(arr_buff[0], arr_buff[1])}\n")




    else:
        T3.delete('1.0', END)
        T3.insert(END, 'Більше 2-ох вибірок\n')


def wilcoxon_test(sample_data, T3):
    s_n_buff = []
    s_check = 0
    for i in range(1, len(sample_data) + 1):
        s_n = f"Вибірка {i}"
        state = sample_data[s_n]["var"].get()
        if state == 1:
            s_n_buff.append(s_n)
            s_check += 1
    if s_check == 2:
        rank_dict = {}
        samp_diff = []
        arr_buff = [sample_data[s_n_buff[0]]["data"], sample_data[s_n_buff[1]]["data"]]
        arr_buff[0] = np.sort(arr_buff[0])
        arr_buff[1] = np.sort(arr_buff[1])
        gen_samp = []
        for i in range(len(arr_buff)):
            for j in range(len(arr_buff[i])):
                gen_samp.append(arr_buff[i][j])
                samp_diff.append({'data': arr_buff[i][j], 'samp': i})

        print("samp diff")
        print(samp_diff)
        sorted_samp_diff = sorted(samp_diff, key=sort_key)
        print("sorted_samp_diff")
        print(sorted_samp_diff)

        samps = []
        for i in range(len(sorted_samp_diff)):
            samps.append(sorted_samp_diff[i]['samp'])

        gen_samp = np.sort(gen_samp)

        for i in range(len(gen_samp)):
            if i > 0 and gen_samp[i - 1] == gen_samp[i]:
                s = 2 * i + 1
                count = 2
                for j in range(i + 1, len(gen_samp)):
                    if gen_samp[j - 1] != gen_samp[j]:
                        break
                    else:
                        count += j + 1
                        s += 1

                for k in range(i - 1, i - 1 + count):
                    rank_dict[k] = {'data': gen_samp[k], 'rank': s / count, 'samp': samps[k]}
            else:
                rank_dict[i] = {'data': gen_samp[i], 'rank': i + 1, 'samp': samps[i]}

        print(rank_dict)

        w_0 = 0
        w_1 = 0

        for i in range(len(gen_samp)):
            if rank_dict[i]['samp'] == 0:
                w_0 += rank_dict[i]['rank']
        for i in range(len(gen_samp)):
            if rank_dict[i]['samp'] == 1:
                w_1 += rank_dict[i]['rank']

        w = min(w_0, w_1)
        mean_w = len(arr_buff[0]) * (len(gen_samp) + 1) / 2

        disp_w = len(arr_buff[0]) * len(arr_buff[1]) * (len(gen_samp) + 1) / 12

        wil_val = (w - mean_w) / np.sqrt(disp_w)

        T3.delete('1.0', END)
        T3.insert(END, 'Критерій Вілкоксона:\n')
        # T3.insert(END, f'w = {round(w, 4)}\n')
        #
        # T3.insert(END, f'mean_w = {round(mean_w, 4)}\n')
        # T3.insert(END, f'disp_w = {round(disp_w, 4)}\n')
        # T3.insert(END, f'w-value = {round(abs(wil_val), 4)}\n')
        if abs(wil_val) < 1.96:
            T3.insert(END, f"w-value = {round(abs(wil_val), 4)} < 1.96, отже, вибірки однорідні\n")
        else:
            T3.insert(END, f"w-value = {round(abs(wil_val), 4)} > 1.96, отже, вибірки неоднорідні\n")

        # T3.insert(END, f'lib w-value = {wilcoxon(arr_buff[0], arr_buff[1])}\n')



    else:
        T3.delete('1.0', END)
        T3.insert(END, 'Більше 2-ох вибірок\n')


def u_test(sample_data, T3):
    s_n_buff = []
    s_check = 0
    for i in range(1, len(sample_data) + 1):
        s_n = f"Вибірка {i}"
        state = sample_data[s_n]["var"].get()
        if state == 1:
            s_n_buff.append(s_n)
            s_check += 1
    if s_check == 2:
        rank_dict = {}
        samp_diff = []
        arr_buff = [sample_data[s_n_buff[0]]["data"], sample_data[s_n_buff[1]]["data"]]
        arr_buff[0] = np.sort(arr_buff[0])
        arr_buff[1] = np.sort(arr_buff[1])
        gen_samp = []
        for i in range(len(arr_buff)):
            for j in range(len(arr_buff[i])):
                gen_samp.append(arr_buff[i][j])
                samp_diff.append({'data': arr_buff[i][j], 'samp': i})

        sorted_samp_diff = sorted(samp_diff, key=sort_key)

        samps = []
        for i in range(len(sorted_samp_diff)):
            samps.append(sorted_samp_diff[i]['samp'])

        gen_samp = np.sort(gen_samp)
        for i in range(len(gen_samp)):
            if i > 0 and gen_samp[i - 1] == gen_samp[i]:
                s = 2 * i + 1
                count = 2
                for j in range(i + 1, len(gen_samp)):
                    if gen_samp[j - 1] != gen_samp[j]:
                        break
                    else:
                        count += j + 1
                        s += 1

                for k in range(i - 1, i - 1 + count):
                    rank_dict[k] = {'data': gen_samp[k], 'rank': s / count, 'samp': samps[k]}
            else:
                rank_dict[i] = {'data': gen_samp[i], 'rank': i + 1, 'samp': samps[i]}

        print(rank_dict)

        w_0 = 0
        w_1 = 0

        for i in range(len(gen_samp)):
            if rank_dict[i]['samp'] == 0:
                w_0 += rank_dict[i]['rank']
        for i in range(len(gen_samp)):
            if rank_dict[i]['samp'] == 1:
                w_1 += rank_dict[i]['rank']

        w = min(w_0, w_1)

        u = len(arr_buff[0]) * len(arr_buff[1]) + len(arr_buff[0]) * (len(arr_buff[0]) - 1) / 2 - abs(w)
        mean_u = len(arr_buff[0]) * len(arr_buff[1]) / 2
        disp_u = len(arr_buff[0]) * len(arr_buff[1]) * (len(gen_samp) + 1) / 12
        u_val = (u - mean_u) / np.sqrt(disp_u)

        T3.delete('1.0', END)
        T3.insert(END, 'U-критерій Манна–Уїтні:\n')
        # T3.insert(END, f'u-value = {round(abs(u_val), 4)}\n')
        if abs(u_val) < 1.96:
            T3.insert(END, f"u-value = {round(abs(u_val), 4)} < 1.96, отже, вибірки однорідні\n")
        else:
            T3.insert(END, f"u-value = {round(abs(u_val), 4)} > 1.96, отже, вибірки неоднорідні\n")
        # T3.insert(END, f'lib u_val = {mannwhitneyu(arr_buff[0], arr_buff[1])}\n')

    else:
        T3.delete('1.0', END)
        T3.insert(END, 'Більше 2-ох вибірок\n')


def diff_mean_rank(sample_data, T3):
    s_n_buff = []
    s_check = 0
    for i in range(1, len(sample_data) + 1):
        s_n = f"Вибірка {i}"
        state = sample_data[s_n]["var"].get()
        if state == 1:
            s_n_buff.append(s_n)
            s_check += 1
    if s_check == 2:
        rank_dict = {}
        samp_diff = []
        arr_buff = [sample_data[s_n_buff[0]]["data"], sample_data[s_n_buff[1]]["data"]]
        arr_buff[0] = np.sort(arr_buff[0])
        arr_buff[1] = np.sort(arr_buff[1])
        gen_samp = []
        for i in range(len(arr_buff[0])):
            gen_samp.append(arr_buff[0][i])
            samp_diff.append({'data': arr_buff[0][i], 'samp': 0})
        for j in range(len(arr_buff[1])):
            gen_samp.append(arr_buff[1][j])

            samp_diff.append({'data': arr_buff[1][j], 'samp': 1})

        sorted_samp_diff = sorted(samp_diff, key=sort_key)

        samps = []
        for i in range(len(sorted_samp_diff)):
            samps.append(sorted_samp_diff[i]['samp'])

        gen_samp = np.sort(gen_samp)
        for i in range(len(gen_samp)):
            if i > 0 and gen_samp[i - 1] == gen_samp[i]:
                s = 2 * i + 1
                count = 2
                for j in range(i + 1, len(gen_samp)):
                    if gen_samp[j - 1] != gen_samp[j]:
                        break
                    else:
                        count += j + 1
                        s += 1

                for k in range(i - 1, i - 1 + count):
                    rank_dict[k] = {'data': gen_samp[k], 'rank': s / count, 'samp': samps[k]}
            else:
                rank_dict[i] = {'data': gen_samp[i], 'rank': i + 1, 'samp': samps[i]}

        print(rank_dict)

        r1_mean = 0
        r2_mean = 0

        for i in range(len(gen_samp)):
            if rank_dict[i]['samp'] == 0:
                r1_mean += rank_dict[i]['rank']
        for i in range(len(gen_samp)):
            if rank_dict[i]['samp'] == 1:
                r2_mean += rank_dict[i]['rank']

        r1_mean = r1_mean / len(arr_buff[0])
        r2_mean = r2_mean / len(arr_buff[0])

        v = (r1_mean - r2_mean) / (
                len(gen_samp) * np.sqrt((len(gen_samp) + 1) / (12 * len(arr_buff[0]) * len(arr_buff[1]))))

        T3.delete('1.0', END)
        T3.insert(END, 'Різниця середніх рангів:\n')
        # T3.insert(END, f'v = {round(abs(v), 4)}\n')
        if abs(v) < 1.64:
            T3.insert(END, f"v = {round(abs(v), 4)} < 1.64, отже, вибірки однорідні\n")
        else:
            T3.insert(END, f"v = {round(abs(v), 4)} > 1.64, отже, вибірки неоднорідні\n")

    else:
        T3.delete('1.0', END)
        T3.insert(END, 'Більше 2-ох вибірок\n')


def h_test(sample_data, T3):
    s_n_buff = []
    s_check = 0
    for i in range(1, len(sample_data) + 1):
        s_n = f"Вибірка {i}"
        state = sample_data[s_n]["var"].get()
        if state == 1:
            s_n_buff.append(s_n)
            s_check += 1

    rank_dict = {}
    samp_diff = []
    arr_buff = []
    for i in range(len(s_n_buff)):
        arr_buff.append(sample_data[s_n_buff[i]]["data"])
    k_s = s_check

    degrees_of_freedom = k_s - 1
    quantile = 0
    if degrees_of_freedom == 1:
        quantile = 0.0039
    elif degrees_of_freedom == 2:
        quantile = 0.103
    elif degrees_of_freedom == 3:
        quantile = 0.352
    elif degrees_of_freedom == 4:
        quantile = 0.711
    elif degrees_of_freedom == 5:
        quantile = 1.15
    elif degrees_of_freedom == 6:
        quantile = 1.64

    k = len(arr_buff)
    for i in range(k):
        arr_buff[i] = np.sort(arr_buff[i])

    gen_samp = []
    for i in range(k_s):
        for j in range(len(arr_buff[i])):
            gen_samp.append(arr_buff[i][j])
            samp_diff.append({'data': arr_buff[i][j], 'samp': i})

    n = len(gen_samp)

    sorted_samp_diff = sorted(samp_diff, key=sort_key)

    samps = []
    for i in range(len(sorted_samp_diff)):
        samps.append(sorted_samp_diff[i]['samp'])

    gen_samp = np.sort(gen_samp)
    for i in range(n):
        if i > 0 and gen_samp[i - 1] == gen_samp[i]:
            s = 2 * i + 1
            count = 2
            for j in range(i + 1, n):
                if gen_samp[j - 1] != gen_samp[j]:
                    break
                else:
                    count += j + 1
                    s += 1

            for k in range(i - 1, i - 1 + count):
                rank_dict[k] = {'data': gen_samp[k], 'rank': s / count, 'samp': samps[k]}
        else:
            rank_dict[i] = {'data': gen_samp[i], 'rank': i + 1, 'samp': samps[i]}

    w_buff = [0 for i in range(k_s)]

    for i in range(k_s):
        for j in range(n):
            if rank_dict[j]['samp'] == i:
                w_buff[i] += rank_dict[j]['rank']

    w_mean = [(w_buff[i] / len(arr_buff[i])) for i in range(k_s)]

    h = 0
    for i in range(k_s):
        h += ((w_mean[i] - (n + 1) / 2) ** 2) / (
                (n + 1) * (n - len(arr_buff[i])) / (12 * len(arr_buff[i]))) * (1 - len(arr_buff[i]) / n)

    T3.delete('1.0', END)
    T3.insert(END, 'H-критерій:\n')
    # T3.insert(END, f'H value = {round(h, 4)}\n')
    if h < quantile:
        T3.insert(END, f"H value = {round(h, 4)} < {quantile}, отже, вибірки однорідні\n")
    else:
        T3.insert(END, f"H value = {round(h, 4)} > {quantile}, отже, вибірки неоднорідні\n")
