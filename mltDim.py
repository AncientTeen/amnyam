import numpy as np
import pandas as pd
import scipy.stats as ss
import statsmodels.api as sm
from scipy.stats import chi2

from correl_regress_Mult import *
from paramMatching import dispCorrMatr
from factor_analysis import factor_anal
from pca import plane_equation_from_points
from claster_class_analysis import distance_matrix, single_linkage_clustering

sorted_eigenvectors = None

eps, Y_viz = 0, 0


def visualization(sample_data, s_n, root, e, Y):
    viz_tabs = Notebook(root)
    plot1 = Frame(viz_tabs)
    viz_tabs.add(plot1, text='Паралельні координати')
    plot2 = Frame(viz_tabs)
    viz_tabs.add(plot2, text='Бульбашкова діаграма')
    plot3 = Frame(viz_tabs)
    viz_tabs.add(plot3, text='3D графік')
    plot4 = Frame(viz_tabs)
    viz_tabs.add(plot4, text='Діагностична діаграма')
    viz_tabs.place(x=250, y=40)

    """parallel coordinates"""
    fig1, ax1 = plt.subplots(figsize=(8, 4), dpi=100)

    plt.grid(color='grey', linestyle='-', linewidth=0.5)

    buff = np.array([sample_data[s_n[i]]["data"] for i in range(len(s_n))])
    # ic(buff)

    min_vals = np.min(buff, axis=1)
    max_vals = np.max(buff, axis=1)
    data_normalized = (buff - min_vals[:, np.newaxis]) / (max_vals[:, np.newaxis] - min_vals[:, np.newaxis])

    num_samples = data_normalized.shape[0]
    num_features = data_normalized.shape[1]

    x_label = [i + 1 for i in range(num_samples)]

    for i in range(num_features):
        plt.plot(x_label, data_normalized[:, i], color='green')

    plt.xticks(x_label)
    plt.title('Паралельні Координати')

    multi_series = FigureCanvasTkAgg(fig1, master=plot1)
    multi_series.get_tk_widget().pack()
    toolbar1 = NavigationToolbar2Tk(multi_series, plot1, pack_toolbar=False)
    toolbar1.update()
    toolbar1.pack()

    """ bubble plot
        only if we have 3 or 4 samples"""
    if len(buff) == 3 or len(buff) == 4:
        fig2, ax2 = plt.subplots(figsize=(8, 4), dpi=100)

        if len(buff) == 3:
            x = buff[0]
            y = buff[1]

            z = buff[2]
            z_norm = []
            for i in range(len(z)):
                if z[i] < 0:
                    z_norm.append(np.sqrt((z[i] + abs(z[i]) + 0.01) / np.pi))
                else:
                    z_norm.append(np.sqrt(z[i] / np.pi))

            plt.scatter(x, y, s=z_norm)

            bubble_chart = FigureCanvasTkAgg(fig2, master=plot2)
            bubble_chart.get_tk_widget().pack()
            toolbar2 = NavigationToolbar2Tk(bubble_chart, plot2, pack_toolbar=False)
            toolbar2.update()
            toolbar2.pack()

    """3D scatter plot
        only if we have 3 or 4 samples"""
    if len(buff) == 3 or len(buff) == 4:
        fig3 = plt.figure(figsize=(8, 4), dpi=100)
        ax3 = fig3.add_subplot(projection='3d')

        if len(buff) == 3:

            x = buff[0]
            y = buff[1]
            z = buff[2]
            ax3.scatter(x, y, z, alpha=1)

            tripleD_plot = FigureCanvasTkAgg(fig3, master=plot3)
            tripleD_plot.get_tk_widget().pack()
            toolbar3 = NavigationToolbar2Tk(tripleD_plot, plot3, pack_toolbar=False)
            toolbar3.update()
            toolbar3.pack()
        elif len(buff) == 4:
            x = buff[0]
            y = buff[1]
            z = buff[2]
            d = buff[3]

            ax3.scatter(x, y, z, s=d * 0.25, alpha=1)

            tripleD_plot = FigureCanvasTkAgg(fig3, master=plot3)
            tripleD_plot.get_tk_widget().pack()
            toolbar3 = NavigationToolbar2Tk(tripleD_plot, plot3, pack_toolbar=False)
            toolbar3.update()
            toolbar3.pack()

    """diagnostic diagram"""
    fig4, ax4 = plt.subplots(figsize=(8, 4), dpi=100)
    ax4.scatter(Y, e, c='red')
    diagnostic_plot = FigureCanvasTkAgg(fig4, master=plot4)
    diagnostic_plot.get_tk_widget().pack()
    toolbar4 = NavigationToolbar2Tk(diagnostic_plot, plot4, pack_toolbar=False)
    toolbar4.update()
    toolbar4.pack()




def outputDataMlt(sample_data, s_n, root, y_sample=1, regBound=[1, 10]):
    global eps, Y_viz

    ic(s_n)
    tabControl = Notebook(root)

    tab1 = Frame(tabControl)
    tab2 = Frame(tabControl)
    tab3 = Frame(tabControl)
    tab4 = Frame(tabControl)
    tab5 = Frame(tabControl)
    tab6 = Frame(tabControl)
    tab7 = Frame(tabControl)

    tabControl.add(tab1, text='Об\'єкти')
    tabControl.add(tab2, text='Протокол')
    tabControl.add(tab3, text='Кореляція')
    tabControl.add(tab4, text='Регресія')
    tabControl.add(tab5, text='МГК')
    tabControl.add(tab6, text='Різноманіття')
    tabControl.add(tab7, text='Факторний аналіз')

    tabControl.place(x=10, y=550)

    T1 = Text(master=tab1, height=10, width=140)
    T2 = Text(master=tab2, height=10, width=140)
    T3 = Text(master=tab3, height=10, width=140)
    T4 = Text(master=tab4, height=10, width=140)
    T5 = Text(master=tab5, height=10, width=140)
    T6 = Text(master=tab6, height=10, width=140)
    T7 = Text(master=tab7, height=10, width=140)

    T1.pack()
    T2.pack()
    T3.pack()
    T4.pack()
    T5.pack()
    T6.pack()
    T7.pack()


    buff = np.round(np.array([sample_data[s_n[i]]["data"] for i in range(len(s_n))]), 4)
    # result = single_linkage_clustering(buff, metric=1, num_clusters=2)

    # ic(result)
    try:
        for i in range(len(buff)):
            T1.insert(END, f"x{i + 1}: {buff[i]}\n")

        buff_mean = [average(buff[i]) for i in range(len(buff))]
        for i in range(len(buff_mean)):
            T2.insert(END, f"Середнє x{i + 1}: {round(buff_mean[i], 4)}\n")

        T2.insert(END, f"\n")
        buff_std_err = [std_err(buff[i], buff_mean[i]) for i in range(len(buff))]
        for i in range(len(buff)):
            T2.insert(END, f"Сер.квадратич. x{i + 1}: {round(buff_std_err[i], 4)}\n")

        disp_corr_matrix = dispCorrMatr(buff)

        T2.insert(END, f"\n")
        T2.insert(END, f"Дисперсійно-Коваріаційна матриця:\n")
        T2.insert(END, f"\t")
        for i in range(len(buff)):
            T2.insert(END, f"x{i + 1}\t\t")
        T2.insert(END, f"\n")
        for i in range(len(buff)):
            T2.insert(END, f"x{i + 1}\t")
            for j in range(len(buff)):
                T2.insert(END, f"{disp_corr_matrix[i][j]}\t\t")
            T2.insert(END, f"\n")

        T2.insert(END, f"\n")
        T2.insert(END, f"Кореляційна матриця:\n")
        T2.insert(END, f"\t")
        corrMatr = corrMatrix(buff)

        for i in range(len(buff)):
            T2.insert(END, f"x{i + 1}\t\t")
        T2.insert(END, f"\n")
        for i in range(len(buff)):
            T2.insert(END, f"x{i + 1}\t")
            for j in range(len(buff)):
                T2.insert(END, f"{corrMatr[i][j]}\t\t")
            T2.insert(END, f"\n")

        partial_corr_coefficients = partCorrCoff(buff)
        # Print the partial correlation coefficients for each pair
        t_quant = 4.698
        if len(buff[0]) > 10 and len(buff[0]) < 120:
            t_quant = 2.105
        elif len(buff[0]) >= 120:
            t_quant = 1.96
    except Exception as e:
        print("The error is: ", e)

    """Кореляція"""
    try:
        T3.insert(END, f"Часткові коефіцієнти кореляції\n")
        t_par_corr = {}
        for pair, coeff in partial_corr_coefficients.items():
            T3.insert(END, f"Частковий коефіцієнт кореляції між {pair}: {coeff:.4f}\n")
            # ic(pair)
            t_par_corr[pair] = coeff * np.sqrt(len(buff[0]) - len(buff) - 2) / np.sqrt(1 - coeff ** 2)
        T3.insert(END, f"\n")

        T3.insert(END, f"Перевірка значущості коефіцієнтів\n")
        for pair, t_coeff in t_par_corr.items():
            if abs(t_coeff) > t_quant:
                T3.insert(END, f"Коефіцієнт між {pair} незначущий: |{t_coeff:.4f}| > {t_quant}\n")
            else:
                T3.insert(END, f"Коефіцієнт між {pair} значущий: |{t_coeff:.4f}| < {t_quant}\n")
        T3.insert(END, f"\n")

        сonf_partial_corr_coefficients = confIntr_partCorrCoff(buff)
        T3.insert(END, f"Довірчі інтервали для коефіцієнтів\n")
        for pair, conf_intr in сonf_partial_corr_coefficients.items():
            # ic(conf_intr)
            T3.insert(END, f"Довірчий інтервал між {pair}: [{conf_intr[0][0]}, {conf_intr[0][1]}]\n")
        T3.insert(END, f"\n")

        T3.insert(END, f"Множинні коефіцієнти кореляції\n")
        multiCorrelationCoff = multCorrCoff(corrMatr)
        f = []
        for i in range(len(multiCorrelationCoff)):
            T3.insert(END, f"r({i + 1}) = {multiCorrelationCoff[i]:.4f}\n")
            f.append(((len(buff[0]) - len(buff) - 1) / len(buff)) * (
                    multiCorrelationCoff[i] ** 2 / (1 - multiCorrelationCoff[i] ** 2)))

        T3.insert(END, f"\n")

        dfn = len(buff)
        dfd = len(buff[0]) - len(buff) - 1
        f_quant = ss.f.ppf(0.95, dfn=dfn, dfd=dfd)

        T3.insert(END, f"Перевірка значущості коефіцієнтів\n")
        for i in range(len(f)):
            if f[i] >= f_quant:
                T3.insert(END, f"r({i + 1}): f = {f[i]:.4f} >= {f_quant:.4f} --> Коефіцієнт значущий\n")
            else:
                T3.insert(END, f"r({i + 1}): f = {f[i]:.4f} < {f_quant:.4f} --> Коефіцієнт незначущий\n")
        T3.insert(END, f"\n")
    except Exception as e:
        print("The error is: ", e)

    """Регресія"""
    try:
        T4.insert(END, f"Лінійна регресія\n")
        A, a_null, C, S_zal, Y_sq, X_sq, Y_low, Y_hat, Y_up, eps, Y_viz = multRegr(buff, y_sample)

        alpha = 0.05
        df = len(buff[0]) - len(buff)
        quantile_lower = chi2.ppf((1 - (1 - alpha)) / 2, df)
        quantile_upper = chi2.ppf((1 + (1 - alpha)) / 2, df)

        lower_bound = S_zal * (len(buff[0]) - len(buff)) / quantile_upper
        upper_bound = S_zal * (len(buff[0]) - len(buff)) / quantile_lower

        T4.insert(END, f"x{y_sample} = {a_null:.4f}")
        xi = 0
        for i in range(len(A)):
            xi += 1
            if i + 1 == y_sample:
                xi += 1
            T4.insert(END, f" + ({A[i]:.4f}) * x{xi}")
        T4.insert(END, f"\n")
        T4.insert(END, f"Оцінка залишкової дисперсії: {S_zal:.4f}\n")
        T4.insert(END, f"\n")

        T4.insert(END, f"\nІнтервальна оцінка параметрів\n")
        t_buff = [t_quant * S_zal * np.sqrt(C[i][i]) for i in range(len(C))]
        xi = 0
        for i in range(len(A)):
            xi += 1
            if i + 1 == y_sample:
                xi += 1
            T4.insert(END, f"a{xi}: [{(A[i] - t_buff[i]):.4f}, {(A[i] + t_buff[i]):4f}]\n")
        T4.insert(END, f"\n")

        T4.insert(END, f"Значущість параметрів регресії\n")
        t_a = [A[i] / np.sqrt(S_zal * C[i][i]) for i in range(len(A))]
        xi = 0
        for i in range(len(A)):
            xi += 1
            if i + 1 == y_sample:
                xi += 1

            if abs(t_a[i]) <= t_quant:
                T4.insert(END, f"t{xi}: |{t_a[i]:.4f}| <= {t_quant} --> параметр а{xi} є незначущим\n")
            else:
                T4.insert(END, f"t{xi}: |{t_a[i]:.4f}| > {t_quant} --> параметр а{xi} є значущим\n")

        T4.insert(END, f"\nСтандартизована оцінка параметрів регресії\n")
        std_A = [(A[i] * X_sq[i]) / Y_sq for i in range(len(A))]
        xi = 0
        for i in range(len(A)):
            xi += 1
            if i + 1 == y_sample:
                xi += 1
            T4.insert(END, f"a{xi} = {std_A[i]:.4f}\n")
        T4.insert(END, f"\n")

        T4.insert(END, f"Коефіцієнт детермінації\n")
        R_square = multiCorrelationCoff[y_sample - 1]
        T4.insert(END, f"R^2 = {R_square * 100:.4f}%\n")
        T4.insert(END, f"\n")

        T4.insert(END, f"Перевірка значущості відтворення\n")
        f_repr = (R_square / (1 - R_square)) * (dfd / dfn)
        if f_repr >= f_quant:
            T4.insert(END, f"f = {f_repr:.4f} >= {f_quant:.4f} --> Регресійна модель значуща\n")
        else:
            T4.insert(END, f"f = {f_repr:.4f} < {f_quant:.4f} --> Регресійна модель не значуща\n")
        T4.insert(END, f"\n")

        T4.insert(END, f"Толерантні межі для залишкової дисперсії\n")
        T4.insert(END, f"{lower_bound:.4f} <= {S_zal:.4f} <= {upper_bound:.4f}\n")
        T4.insert(END, f"\n")

        low = regBound[0] - 1
        up = regBound[1]
        # ic(low, up)
        T4.insert(END, f"Довірчий інтервал для значення регресії\n")
        for i in range(low, up):
            T4.insert(END, f"{Y_low[i]:.4f} <= {Y_hat[i]:.4f} <= {Y_up[i]:.4f}\n")
        T4.insert(END, f"\n")

    except Exception as e:
        print("The error is: ", e)

    """Метод головних компонент"""
    try:
        sorted_eigenvalues, sorted_eigenvectors = eigenValsVectors(sample_data)
        T5.insert(END, f"Результати МГК:\n")

        for i in range(len(buff)):
            T5.insert(END, f"\tx{i + 1}`\t")
        T5.insert(END, f"\n")
        for i in range(len(buff)):
            T5.insert(END, f"x{i + 1}\t")
            for j in range(len(buff)):
                T5.insert(END, f"{sorted_eigenvectors[i][j]:.4f}\t\t")
            T5.insert(END, f"\n")

        T5.insert(END, f"\nВласні числа: \t")
        for i in range(len(buff)):
            T5.insert(END, f"{sorted_eigenvalues[i]:.4f}\t\t")
        T5.insert(END, f"\n")

        T5.insert(END, f"\n% на напрям:  \t")
        percentDir = [sorted_eigenvalues[i] / len(sorted_eigenvalues) for i in range(len(sorted_eigenvalues))]
        for i in range(len(buff)):
            T5.insert(END, f"{percentDir[i] * 100:.4f}%\t\t")
        T5.insert(END, f"\n")

        T5.insert(END, f"\nНакопичений %:\t")
        perc = 0
        for i in range(len(buff)):
            perc += percentDir[i] * 100
            T5.insert(END, f"{perc:.4f}%\t\t")
        T5.insert(END, f"\n\n")
    except Exception as e:
        print("The error is: ", e)

    """Різноманіття"""
    try:
        if len(buff) == 3:
            a, b, c, d = plane_equation_from_points(buff)
            A, a_null, C, S_zal, Y_sq, X_sq, Y_low, Y_hat, Y_up, eps, Y_viz = multRegr(buff, y_sample)

            T6.insert(END, f"Рівняння площини:\n")
            T6.insert(END, f"{a:.4f} * x1 + ({b:.4f}) * x2 + ({c:.4f}) * x3 + ({-d:.4f}) = 0:\n")
            T6.insert(END, f"x1 = {d / a:.4f} + ({-b / a:.4f}) * x2 + ({-c / a:.4f}) * x3\n\n")

            T6.insert(END, f"Лінійна регресія:\n")
            T6.insert(END, f"x{y_sample} = {a_null:.4f}")
            xi = 0
            for i in range(len(A)):
                xi += 1
                if i + 1 == y_sample:
                    xi += 1
                T6.insert(END, f" + ({A[i]:.4f}) * x{xi}")

        else:
            pass
    except Exception as e:
        print("The error is: ", e)

    # correlation_matrix = np.array([
    #     [1, 0.559, 0.165, 0.089, -0.035, 0.066, 0.293, -0.149, -0.071, 0.118, 0.396],
    #     [0.559, 1, 0.076, -0.083, 0.016, 0.201, 0.613, -0.223, 0.016, -0.057, 0.783],
    #     [0.165, 0.076, 1, -0.211, -0.278, -0.286, -0.086, -0.192, -0.349, -0.308, -0.205],
    #     [0.089, -0.083, -0.211, 1, 0.227, -0.053, -0.376, 0.503, 0.128, 0.291, 0.035],
    #     [-0.035, 0.016, -0.278, 0.227, 1, 0.796, 0.503, 0.184, 0.832, 0.738, 0.505],
    #     [0.066, 0.201, -0.286, -0.053, 0.796, 1, 0.508, -0.137, 0.872, 0.675, 0.686],
    #     [0.293, 0.613, -0.086, -0.376, 0.503, 0.508, 1, 0.096, 0.495, 0.400, 0.783],
    #     [-0.149, -0.223, -0.192, 0.503, 0.184, -0.137, 0.096, 1, 0.261, 0.252, -0.005],
    #     [-0.071, 0.016, -0.349, 0.128, 0.832, 0.872, 0.495, 0.261, 1, 0.761, 0.614],
    #     [0.118, -0.057, -0.308, 0.291, 0.738, 0.675, 0.400, 0.252, 0.761, 1, 0.405],
    #     [0.396, 0.783, -0.205, 0.035, 0.505, 0.686, 0.783, -0.005, 0.614, 0.405, 1]
    # ])

    # ic(factor_anal(correlation_matrix))


    """Факторний аналіз"""
    try:
        A, w, itr = factor_anal(corrMatr)
        T7.insert(END, f"Матриця факторного відображення:\n")

        for i in range(w):
            T7.insert(END, f"\tF{i + 1}\t")
        T7.insert(END, f"\n")
        for i in range(len(A)):
            T7.insert(END, f"x{i + 1}\t")
            for j in range(w):
                T7.insert(END, f"{A[i][j]:.4f}\t\t")
            T7.insert(END, f"\n")

    except Exception as e:
        print("The error is: ", e)

    return eps, Y_viz


def pca_two_frwd(sample_data):
    global fi
    fi = transIndep(sample_data)


def pca_two_bck(sample_data):
    backTransIndep(sample_data, fi)


def pcaFnDim(sample_data, sample_menu, sample_checkbuttons):
    global sorted_eigenvectors

    sorted_eigenvalues, sorted_eigenvectors = eigenValsVectors(sample_data)
    pca_forw_nDim(sample_data, sorted_eigenvectors, sample_menu, sample_checkbuttons)


def pcaBnDim(sample_data, sample_menu, sample_checkbuttons):
    user_input1 = simpledialog.askstring("N", "Введіть кількість компонент:")
    if user_input1 != "":
        n = int(user_input1)
        # sorted_eigenvalues, sorted_eigenvectors = eigenValsVectors(sample_data)

        pca_back_nDim(sample_data, sorted_eigenvectors, n, sample_menu, sample_checkbuttons)
