from paramFuncs import *


def linearCorrelation(x, y):
    n = len(x)
    avr_x = average(x)
    avr_y = average(y)
    std_err_x = std_err(x, avr_x)
    std_err_y = std_err(y, avr_y)
    avrProd = averageProduct(x, y)

    correl = round((n / (n - 1)) * ((avrProd - avr_x * avr_y) / (std_err_x * std_err_y)), 4)
    # correl = (n / (n - 1)) * ((avrProd - avr_x * avr_y) / (std_err_x * std_err_y))
    return correl


def confInterCorrel(x, y):
    n = len(x)
    t = 0
    if n - 1 < 120:
        if n - 1 == 69:
            t = 1.995
        elif n - 1 == 24:
            t = 2.06
    else:
        t = 1.96

    correl = linearCorrelation(x, y)

    inf = round(correl + ((correl * (1 - correl ** 2)) / (2 * n)) - ((t * (1 - correl ** 2)) / (np.sqrt(n - 1))), 4)
    sup = round(correl + ((correl * (1 - correl ** 2)) / (2 * n)) + ((t * (1 - correl ** 2)) / (np.sqrt(n - 1))), 4)
    return inf, sup


def densitySnDimFunc(x, y):
    n = len(x)
    if n < 100:
        b = round(n ** (1 / 2))
        if b % 2 == 0:
            b -= 1
    else:
        b = round(n ** (1 / 3))
        if b % 2 == 0:
            b -= 1

    m_val_x = np.arange(min(x), max(x), (max(x) - min(x)) / b)
    m_val_y = np.arange(min(y), max(y), (max(y) - min(y)) / b)

    m_val_x = np.append(m_val_x, max(x))
    m_val_y = np.append(m_val_y, max(y))

    avr_val_samp_x = [[0 for i in range(b)] for j in range(b)]
    avr_val_samp_y = [[0 for i in range(b)] for j in range(b)]
    for i in range(1, b + 1):
        for j in range(1, b + 1):
            min_x = 0
            max_x = 0
            min_y = 0
            max_y = 0
            for k in range(n):
                if m_val_x[i - 1] <= x[k] < m_val_x[i] and m_val_y[j - 1] <= y[k] < m_val_y[j]:
                    if min_x > x[k]:
                        min_x = x[k]
                    if min_y > y[k]:
                        min_y = y[k]
                    if max_x < x[k]:
                        max_x = x[k]
                    if max_y < y[k]:
                        max_x = y[k]
            avr_val_samp_x[i - 1][j - 1] = (max_x + min_x) / 2
            avr_val_samp_y[i - 1][j - 1] = (max_y + min_y) / 2

    avr_x = average(x)
    avr_y = average(y)
    std_err_x = std_err(x, avr_x)
    std_err_y = std_err(y, avr_y)
    correl = linearCorrelation(x, y)

    func_val = [[0 for i in range(b)] for j in range(b)]
    for i in range(b):
        for j in range(b):
            func_val[i][j] = 1 / (2 * np.pi * std_err_x * std_err_y * np.sqrt(1 - correl ** 2)) * np.exp(
                (-1 / (2 * (1 - correl ** 2))) * (((avr_val_samp_x[i][j] - avr_x) / std_err_x) ** 2 - 2 * correl * (
                        (avr_val_samp_x[i][j] - avr_x) / std_err_x) * ((avr_val_samp_y[i][j] - avr_y) / std_err_y) + (
                                                          (avr_val_samp_y[i][j] - avr_y) / std_err_y) ** 2))

    return func_val


def adequateReproduction(x, y):
    try:
        n = len(x)
        if n < 100:
            b = round(n ** (1 / 2))
            if b % 2 == 0:
                b -= 1
        else:
            b = round(n ** (1 / 3))
            if b % 2 == 0:
                b -= 1

        p_val_arr = [[0 for i in range(b)] for j in range(b)]

        m_val_x = np.arange(min(x), max(x), (max(x) - min(x)) / b)
        m_val_y = np.arange(min(y), max(y), (max(y) - min(y)) / b)

        m_val_x = np.append(m_val_x, max(x))
        m_val_y = np.append(m_val_y, max(y))
        print(f"m_val_x: {m_val_x}")
        print(f"m_val_y: {m_val_y}")
        for i in range(1, b + 1):
            for j in range(1, b + 1):
                count = 0
                for k in range(n):
                    if m_val_x[i - 1] <= x[k] < m_val_x[i] and m_val_y[j - 1] <= y[k] < m_val_y[j]:
                        count += 1
                p_val_arr[i - 1][j - 1] = count / n

        dens_func_val = densitySnDimFunc(x, y)

        xi_val = 0
        for i in range(b):
            for j in range(b):
                xi_val += (p_val_arr[i][j] - dens_func_val[i][j]) ** 2 / (dens_func_val[i][j] * m_val_x[j] * m_val_y[j])
        xi_val = round(xi_val, 4)
        return xi_val

    except Exception as e:
        print("The error is: ", e)


def t_corr_test(x, y):
    n = len(x)
    correl = linearCorrelation(x, y)

    t_corr = round((correl * np.sqrt(n - 2)) / np.sqrt(1 - correl ** 2), 4)
    return t_corr


def correlationRelation(x, y):
    try:
        n = len(x)
        if n < 100:
            b = round(n ** (1 / 2))
            if b % 2 == 0:
                b -= 1
        else:
            b = round(n ** (1 / 3))
            if b % 2 == 0:
                b -= 1

        m_val_y = np.arange(min(y), max(y), (max(y) - min(y)) / b)
        m_val_y = np.append(m_val_y, max(y))
        print(f'm_val_y {m_val_y}')
        print(f'm_ax_y {max(y)}')
        avr_y = average(y)
        y_class_avr = []
        m_arr = []

        for j in range(1, b + 1):
            count = 0
            y_avr = 0
            for k in range(n):
                if m_val_y[j - 1] <= y[k] <= m_val_y[j]:
                    y_avr += y[k]
                    count += 1
            y_avr = y_avr / count
            m_arr.append(count)
            y_class_avr.append(y_avr)

        p_coff_numerator = 0
        for i in range(b):
            p_coff_numerator += m_arr[i] * (y_class_avr[i] - avr_y) ** 2

        p_coff_denominator = 0
        for j in range(1, b + 1):
            for k in range(n):
                if m_val_y[j - 1] <= y[k] <= m_val_y[j]:
                    p_coff_denominator += (y[k] - avr_y) ** 2

        p_coff = round(np.sqrt(p_coff_numerator / p_coff_denominator), 4)

        return p_coff
    except Exception as e:
        print("The error is: ", e)




def tValCorrRel(x, y):
    try:
        n = len(x)
        if n < 100:
            b = round(n ** (1 / 2))
            if b % 2 == 0:
                b -= 1
        else:
            b = round(n ** (1 / 3))
            if b % 2 == 0:
                b -= 1
        p_coff = correlationRelation(x, y)

        f_val = round((p_coff ** 2 / (1 - p_coff ** 2)) * ((n - b) / (b - 1)), 4)
        return f_val
    except Exception as e:
        print("The error is: ", e)



def fValCorrRel(x, y):
    try:
        n = len(x)
        p_coff = correlationRelation(x, y)

        t_val = round(p_coff * np.sqrt(n - 2) / np.sqrt(1 - p_coff ** 2), 4)
        return t_val
    except Exception as e:
        print("The error is: ", e)

def spearmanCoff(x, y):
    n = len(x)

    rank_x = [sorted(x).index(elem) + 1 for elem in x]
    rank_y = [sorted(y).index(elem) + 1 for elem in y]

    d = [rx - ry for rx, ry in zip(rank_x, rank_y)]

    tau_val = round(1 - (6 * sum([d_i ** 2 for d_i in d])) / (n * (n ** 2 - 1)), 4)

    return tau_val

def kendallCoff(x, y):
    n = len(x)
    concordant, discordant = 0, 0

    for i in range(n - 1):
        for j in range(i + 1, n):
            if (x[i] < x[j] and y[i] < y[j]) or (x[i] > x[j] and y[i] > y[j]):
                concordant += 1
            elif (x[i] < x[j] and y[i] > y[j]) or (x[i] > x[j] and y[i] < y[j]):
                discordant += 1

    tau = round((concordant - discordant) / (concordant + discordant), 4)

    return tau


def combTable2x2(x, y):
    n = len(x)
    avr_x = average(x)
    avr_y = average(y)
    bin_arr_x = []
    bin_arr_y = []

    for i in range(n):
        if x[i] - avr_x > 0:
            bin_arr_x.append(1)
        else:
            bin_arr_x.append(0)
        if y[i] - avr_y > 0:
            bin_arr_y.append(1)
        else:
            bin_arr_y.append(0)

    n00 = 0
    n01 = 0
    n10 = 0
    n11 = 0

    for i in range(n):
        if bin_arr_x[i] == 0 and bin_arr_y[i] == 0:
            n00 += 1
        elif bin_arr_x[i] == 1 and bin_arr_y[i] == 0:
            n01 += 1
        elif bin_arr_x[i] == 0 and bin_arr_y[i] == 1:
            n10 += 1
        elif bin_arr_x[i] == 1 and bin_arr_y[i] == 1:
            n11 += 1

    m0 = n00 + n10
    m1 = n01 + n11
    n0 = n00 + n01
    n1 = n11 + n10

    return n00, n01, n10, n11, m0, m1, n0, n1

def fecherInd(x, y):
    n = len(x)
    avr_x = average(x)
    avr_y = average(y)
    bin_arr_x = []
    bin_arr_y = []

    for i in range(n):
        if x[i] - avr_x > 0:
            bin_arr_x.append(1)
        else:
            bin_arr_x.append(0)
        if y[i] - avr_y > 0:
            bin_arr_y.append(1)
        else:
            bin_arr_y.append(0)
    # print('combinations table')
    # print(bin_arr_x)
    # print(len(bin_arr_x))
    # print(bin_arr_y)
    # print(len(bin_arr_y))

    n00 = 0
    n01 = 0
    n10 = 0
    n11 = 0

    for i in range(n):
        if bin_arr_x[i] == 0 and bin_arr_y[i] == 0:
            n00 += 1
        elif bin_arr_x[i] == 1 and bin_arr_y[i] == 0:
            n01 += 1
        elif bin_arr_x[i] == 0 and bin_arr_y[i] == 1:
            n10 += 1
        elif bin_arr_x[i] == 1 and bin_arr_y[i] == 1:
            n11 += 1
    # print(f'n00 = {n00}')
    # print(f'n01 = {n01}')
    # print# print(f'n11 = {n11}')

    I = (n00 + n11 - n10 - n01) / (n00 + n11 + n01 + n10)

    return I


def fiCoff(x, y):
    n = len(x)
    avr_x = average(x)
    avr_y = average(y)
    bin_arr_x = []
    bin_arr_y = []

    for i in range(n):
        if x[i] - avr_x > 0:
            bin_arr_x.append(1)
        else:
            bin_arr_x.append(0)
        if y[i] - avr_y > 0:
            bin_arr_y.append(1)
        else:
            bin_arr_y.append(0)

    n00 = 0
    n01 = 0
    n10 = 0
    n11 = 0

    for i in range(n):
        if bin_arr_x[i] == 0 and bin_arr_y[i] == 0:
            n00 += 1
        elif bin_arr_x[i] == 1 and bin_arr_y[i] == 0:
            n01 += 1
        elif bin_arr_x[i] == 0 and bin_arr_y[i] == 1:
            n10 += 1
        elif bin_arr_x[i] == 1 and bin_arr_y[i] == 1:
            n11 += 1

    m0 = n00 + n10
    m1 = n01 + n11
    n0 = n00 + n01
    n1 = n11 + n10
    fi = round((n00 * n11 - n01 * n10) / np.sqrt(n0 * n1 * m0 * m1), 4)

    return fi


def yuleCoff(x, y):
    try:
        n = len(x)
        avr_x = average(x)
        avr_y = average(y)
        bin_arr_x = []
        bin_arr_y = []

        for i in range(n):
            if x[i] - avr_x > 0:
                bin_arr_x.append(1)
            else:
                bin_arr_x.append(0)
            if y[i] - avr_y > 0:
                bin_arr_y.append(1)
            else:
                bin_arr_y.append(0)

        n00 = 0
        n01 = 0
        n10 = 0
        n11 = 0

        for i in range(n):
            if bin_arr_x[i] == 0 and bin_arr_y[i] == 0:
                n00 += 1
            elif bin_arr_x[i] == 1 and bin_arr_y[i] == 0:
                n01 += 1
            elif bin_arr_x[i] == 0 and bin_arr_y[i] == 1:
                n10 += 1
            elif bin_arr_x[i] == 1 and bin_arr_y[i] == 1:
                n11 += 1

        Y = round((np.sqrt(n00 * n11) - np.sqrt(n10 * n01)) / (np.sqrt(n00 * n11) + np.sqrt(n01 * n10)), 4)
        Q = round((2 * Y) / ((1 + Y) ** 2), 4)

        """check for significant and for confidence interval"""
        sQ = (1 / 2) * (1 - Q ** 2) * np.sqrt(1 / n00 + 1 / n01 + 1 / n10 + 1 / n11)
        sY = (1 / 4) * (1 - Y ** 2) * np.sqrt(1 / n00 + 1 / n01 + 1 / n10 + 1 / n11)
        uQ = abs(round(Q / sQ, 4))
        uY = abs(round(Y / sY, 4))

        return Q, Y, uQ, uY
    except Exception as e:
        print("The error is: ", e)



def pirsonCorr(x, y, n, m):
    N_len = len(x)
    if N_len < 100:
        b = round(N_len ** (1 / 2))
        if b % 2 == 0:
            b -= 1
    else:
        b = round(N_len ** (1 / 3))
        if b % 2 == 0:
            b -= 1

    m_val_x = np.arange(min(x), max(x), (max(x) - min(x)) / m)
    m_val_y = np.arange(min(y), max(y), (max(y) - min(y)) / n)

    m_val_x = np.append(m_val_x, max(x))
    m_val_y = np.append(m_val_y, max(y))

    combTable = [[0 for i in range(m)] for j in range(n)]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            count = 0
            for k in range(N_len):
                if m_val_x[j - 1] <= x[k] < m_val_x[j] and m_val_y[i - 1] <= y[k] < m_val_y[i]:
                    count += 1

            combTable[i - 1][j - 1] = count

    m_arr = []
    for i in range(m):
        columnSum = 0
        for j in range(n):
            columnSum += combTable[j][i]
        m_arr.append(columnSum)

    n_arr = []
    for i in range(n):
        rowSum = 0
        for j in range(m):
            rowSum += combTable[i][j]
        n_arr.append(rowSum)

    N = [[0 for i in range(m)] for j in range(n)]
    for i in range(n):
        for j in range(m):
            N[i][j] = (n_arr[i] * m_arr[j]) / N_len

    chiValCoff = 0
    for i in range(n):
        for j in range(m):
            chiValCoff += (combTable[i][j] - N[i][j]) ** 2 / N[i][j]

    chi = round(np.sqrt(chiValCoff ** 2 / (N_len + chiValCoff ** 2)), 4)

    return chi


def kendallCorr(x, y, n, m):
    if n != m:
        return 1

    N_len = len(x)
    if N_len < 100:
        b = round(N_len ** (1 / 2))
        if b % 2 == 0:
            b -= 1
    else:
        b = round(N_len ** (1 / 3))
        if b % 2 == 0:
            b -= 1

    m_val_x = np.arange(min(x), max(x), (max(x) - min(x)) / m)
    m_val_y = np.arange(min(y), max(y), (max(y) - min(y)) / n)

    m_val_x = np.append(m_val_x, max(x))
    m_val_y = np.append(m_val_y, max(y))

    combTable = [[0 for i in range(m)] for j in range(n)]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            count = 0
            for k in range(N_len):
                if m_val_x[j - 1] <= x[k] < m_val_x[j] and m_val_y[i - 1] <= y[k] < m_val_y[i]:
                    count += 1

            combTable[i - 1][j - 1] = count

    m_arr = []
    for i in range(m):
        columnSum = 0
        for j in range(n):
            columnSum += combTable[j][i]
        m_arr.append(columnSum)

    n_arr = []
    for i in range(n):
        rowSum = 0
        for j in range(m):
            rowSum += combTable[i][j]
        n_arr.append(rowSum)

    P = 0
    for i in range(n):
        n_sum_external = 0
        for j in range(m):
            n_sum_external += combTable[i][j]
            n_sum_internal = 0
            for k in range(i + 1, n):
                for l in range(j + 1, m):
                    n_sum_internal += combTable[k][l]

            P += n_sum_external * n_sum_internal

    Q = 0
    for i in range(n):
        n_sum_external = 0
        for j in range(m):
            n_sum_external += combTable[i][j]
            n_sum_internal = 0
            for k in range(i + 1, n):
                for l in range(1, j - 1):
                    n_sum_internal += combTable[k][l]

            Q += n_sum_external * n_sum_internal

    T1 = 0
    for i in range(n):
        T1 += n_arr[i] * (n_arr[i] - 1)
    T1 = T1 * (1 / 2)

    T2 = 0
    for j in range(m):
        T1 += m_arr[j] * (m_arr[j] - 1)
    T2 = T2 * (1 / 2)

    tau = (P - Q) / np.sqrt(((1 / 2) * N_len * (N_len - 1) - T1) * ((1 / 2) * N_len * (N_len - 1) - T2))

    return tau


def stewardCorr(x, y, n, m):
    if n == m:
        return 1

    N_len = len(x)
    if N_len < 100:
        b = round(N_len ** (1 / 2))
        if b % 2 == 0:
            b -= 1
    else:
        b = round(N_len ** (1 / 3))
        if b % 2 == 0:
            b -= 1

    m_val_x = np.arange(min(x), max(x), (max(x) - min(x)) / m)
    m_val_y = np.arange(min(y), max(y), (max(y) - min(y)) / n)

    m_val_x = np.append(m_val_x, max(x))
    m_val_y = np.append(m_val_y, max(y))

    combTable = [[0 for i in range(m)] for j in range(n)]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            count = 0
            for k in range(N_len):
                if m_val_x[j - 1] <= x[k] < m_val_x[j] and m_val_y[i - 1] <= y[k] < m_val_y[i]:
                    count += 1

            combTable[i - 1][j - 1] = count

    m_arr = []
    for i in range(m):
        columnSum = 0
        for j in range(n):
            columnSum += combTable[j][i]
        m_arr.append(columnSum)

    n_arr = []
    for i in range(n):
        rowSum = 0
        for j in range(m):
            rowSum += combTable[i][j]
        n_arr.append(rowSum)

    P = 0
    for i in range(n):
        n_sum_external = 0
        for j in range(m):
            n_sum_external += combTable[i][j]
            n_sum_internal = 0
            for k in range(i + 1, n):
                for l in range(j + 1, m):
                    n_sum_internal += combTable[k][l]

            P += n_sum_external * n_sum_internal

    Q = 0
    for i in range(n):
        n_sum_external = 0
        for j in range(m):
            n_sum_external += combTable[i][j]
            n_sum_internal = 0
            for k in range(i + 1, n):
                for l in range(1, j - 1):
                    n_sum_internal += combTable[k][l]

            Q += n_sum_external * n_sum_internal

    tau = (2 * (P - Q) * min(m, n)) / (N_len ** 2 * (min(m, n) - 1))


    return tau

