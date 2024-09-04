from tkinter import *
from tkinter import simpledialog
from tkinter.ttk import Notebook

import matplotlib.pyplot as plt

from linRegression import *
from nonlinRegression import *
from quasilinRegression import *
from pca import *

x = None
y = None
global tab7
global T3
global root2
""""""


def display_function_list(event, root, fig1, ax1, x, y):
    if event.button == 3:  # Check if it's a right-click event
        function_list = ["Лінійна регресія", "Метод Тейла", "Толерантні межі", "Прогноз нового спостереження",
                         "Довірчі інтервали", "Параболічна регресія", "Толеларнтні межі для параболічної регресії",
                         "Довірчі інтервали для параболічної регресії",
                         "Прогноз нового спостереження для параболічної регресії", "Квазілінійна регресія",
                         "Толеларнтні межі для Квазілінійної регресії",
                         "Довірчі інтервали для Квазілінійної регресії",
                         "Прогноз нового спостереження для квазілінійної регресії",
                         "Очистити"]  # Add your function names here

        # Create a new Toplevel window
        function_list_window = Toplevel(root)

        # Create a Listbox and populate it with function options
        listbox = Listbox(function_list_window, selectmode=SINGLE, width=40, height=15)
        for function_name in function_list:
            listbox.insert(END, function_name)

        # Function to handle item selection from the Listbox
        def select_function(event):
            selected_index = listbox.curselection()
            if selected_index:
                selected_function = function_list[selected_index[0]]

                # Depending on the selected function, you can call the corresponding functionality
                if selected_function == "Лінійна регресія":
                    # Call your function 1
                    linRegressionMNK(x, y, fig1, ax1)
                    pass
                elif selected_function == "Метод Тейла":
                    linRegressionTeyla(x, y, fig1, ax1)
                    pass
                elif selected_function == "Толерантні межі":
                    toleranceLim(x, y, fig1, ax1)
                    pass
                elif selected_function == "Прогноз нового спостереження":
                    predictNewObs(x, y, fig1, ax1)
                    pass
                elif selected_function == "Довірчі інтервали":
                    confInterLinRegr(x, y, fig1, ax1)
                    pass
                elif selected_function == "Параболічна регресія":
                    parabolicRegression(x, y, fig1, ax1)
                    pass
                elif selected_function == "Толеларнтні межі для параболічної регресії":
                    toleranceLimNonLin(x, y, fig1, ax1)
                    pass
                elif selected_function == "Довірчі інтервали для параболічної регресії":
                    confInterNonLinRegr(x, y, fig1, ax1)
                    pass
                elif selected_function == "Прогноз нового спостереження для параболічної регресії":
                    predictNewObsNonLin(x, y, fig1, ax1)
                    pass
                elif selected_function == "Квазілінійна регресія":
                    quasilinRegr(x, y, fig1, ax1)
                    pass
                elif selected_function == "Толеларнтні межі для Квазілінійної регресії":
                    toleranceLimQuasilinRegr(x, y, fig1, ax1)
                    pass

                elif selected_function == "Довірчі інтервали для Квазілінійної регресії":
                    confInterQuasilinRegr(x, y, fig1, ax1)
                    pass
                elif selected_function == "Прогноз нового спостереження для квазілінійної регресії":
                    predictNewObsQuasilinRegr(x, y, fig1, ax1)
                    pass
                elif selected_function == "Очистити":
                    variation_series(x, y, root)
                    pass
                # Add more conditions as needed

                # Close the function_list_window after a selection is made
                function_list_window.destroy()

        # Bind the item selection event to the select_function function
        listbox.bind('<ButtonRelease-1>', select_function)

        listbox.pack()


# Modify your existing variation_series function
def variation_series(x, y, root):
    fig1, ax1 = plt.subplots(figsize=(5, 4), dpi=100)

    plt.grid(color='grey', linestyle='-', linewidth=0.5)
    plt.xlabel('X')
    plt.ylabel('Y')

    ax1.scatter(x=x, y=y, s=5)

    # Bind the left mouse click event to display points
    def on_left_click(event):
        print(f"Left click at ({event.xdata}, {event.ydata})")

    # Bind the right mouse click event to display_function_list
    def on_right_click(event):
        display_function_list(event, root, fig1, ax1, x, y)

    fig1.canvas.mpl_connect('button_press_event', on_left_click)
    fig1.canvas.mpl_connect('button_press_event', on_right_click)

    variation_series_out = FigureCanvasTkAgg(fig1, master=root)
    variation_series_out.get_tk_widget().place(x=50, y=40)
    toolbar = NavigationToolbar2Tk(variation_series_out, root, pack_toolbar=False)
    toolbar.update()
    toolbar.place(x=150, y=450)




def SndDimDens(x, y, root, cls=None):

    if cls:
        b = cls
    else:
        if len(x) < 100:
            b = round((len(x) ** (1 / 2)))
            if b % 2 == 0:
                b -= 1
        else:
            b = round((len(x) ** (1 / 3)))
            if b % 2 == 0:
                b -= 1

    fig2, ax2 = plt.subplots(figsize=(5, 4), dpi=100)

    n = len(x)

    plt.hist2d(x, y, bins=(b, b), cmap=plt.cm.jet, weights=np.ones_like(x) / n)
    plt.colorbar()

    variation_series = FigureCanvasTkAgg(fig2, master=root)
    variation_series.get_tk_widget().place(x=800, y=40)
    toolbar = NavigationToolbar2Tk(variation_series, root, pack_toolbar=False)
    toolbar.update()
    toolbar.place(x=900, y=450)



def classes(x, y, root):
    global cls
    user_input = simpledialog.askstring("Classes", "Введіть кількість класів")

    if user_input != "":
        cls = int(user_input)

    print("x: ", x)
    print("y: ", y)
    SndDimDens(x, y, root, cls)


def remAnom(x, y, root, cls=None):
    x, y = removeAnomalous(x, y)

    variation_series(x, y, root)
    SndDimDens(x, y, root, cls)


def logData(x, y, root, cls=None):
    print('len(x) ', len(x))
    print('len(y) ', len(y))
    x, y = logs(x, y)
    print('len(x) ', len(x))
    print('len(y) ', len(y))
    variation_series(x, y, root)
    SndDimDens(x, y, root, cls)
    outputDataScdDim(x, y, root)


def outputDataScdDim(x, y, root):

    tabControl = Notebook(root)

    tab1 = Frame(tabControl)
    tab2 = Frame(tabControl)
    tab3 = Frame(tabControl)
    tab4 = Frame(tabControl)
    tab5 = Frame(tabControl)
    tab6 = Frame(tabControl)
    global tab7
    tab7 = Frame(tabControl)

    tabControl.add(tab1, text='Об\'єкти')
    tabControl.add(tab2, text='Протокол')
    tabControl.add(tab3, text='Кореляційний аналіз')
    tabControl.add(tab4, text='Лінійний регреійний аналіз')
    tabControl.add(tab5, text='Нелінійний регреійний аналіз')
    tabControl.add(tab6, text='Квазілінійний регреійний аналіз')
    tabControl.add(tab7, text='Таблиця NxM')

    tabControl.place(x=10, y=500)

    T1 = Text(master=tab1, height=10, width=100)
    T2 = Text(master=tab2, height=10, width=100)
    T3 = Text(master=tab3, height=10, width=100)
    T4 = Text(master=tab4, height=10, width=100)
    T5 = Text(master=tab5, height=10, width=100)
    T6 = Text(master=tab6, height=10, width=100)
    T7 = Text(master=tab7, height=10, width=100)

    T1.pack()
    T2.pack()
    T3.pack()
    T4.pack()
    T5.pack()
    T6.pack()
    T7.pack()

    x_sort = np.sort(x)
    y_sort = np.sort(y)

    T1.insert(END, f"x: {x_sort}")
    T1.insert(END, f"\ny: {y_sort}")

    n_x = len(x)
    n_y = len(y)
    n = len(x)

    if n < 100:
        b = round(n ** (1 / 2))
        if b % 2 == 0:
            b -= 1
    else:
        b = round(n ** (1 / 3))
        if b % 2 == 0:
            b -= 1
    try:
        T2.insert(END, 'Характеристика\t\t\t' + 'INF\t\t' + 'Значення\t\t' + 'SUP\t\t' + 'SKV\n')

        avr_x = average(x)
        avr_y = average(y)
        std_err_x = std_err(x, avr_x)
        std_err_y = std_err(y, avr_y)
        avrIntr_x = confInterAvr(x)
        avrIntr_y = confInterAvr(y)
        T2.insert(END,
                  'Середнє X: \t\t\t' + str(avrIntr_x[0]) + '\t\t' + str(avr_x) + '\t\t' + str(avrIntr_x[1]) + '\t\t' +
                  str(round((std_err_x / n_x) ** (1 / 2), 4)) + '\n')

        std_err_x_inter = confInterStdErr(x)
        T2.insert(END,
                  'Сер.квадратич. X: \t\t\t' + str(std_err_x_inter[0]) + '\t\t' + str(std_err_x) + '\t\t' + str(
                      std_err_x_inter[1]) + '\t\t' +
                  str(round(std_err_x * (2 / (n_x - 1)) ** (1 / 4), 4)) + '\n')

        T2.insert(END,
                  'Середнє Y: \t\t\t' + str(avrIntr_y[0]) + '\t\t' + str(avr_y) + '\t\t' + str(avrIntr_y[1]) + '\t\t' +
                  str(round((std_err_y / n_y) ** (1 / 2), 4)) + '\n')

        std_err_y_inter = confInterStdErr(y)
        T2.insert(END,
                  'Сер.квадратич. Y: \t\t\t' + str(std_err_y_inter[0]) + '\t\t' + str(std_err_y) + '\t\t' + str(
                      std_err_y_inter[1]) + '\t\t' +
                  str(round(std_err_y * (2 / (n_y - 1)) ** (1 / 4), 4)) + '\n')

        """Correlation analysis"""
        correl = linearCorrelation(x, y)
        correl_inter = confInterCorrel(x, y)
        T3.insert(END, '\t\t\t' + 'INF\t\t' + 'Значення\t\t' + 'SUP\t\t' + 'SKV\n')

        T3.insert(END,
                  'Лін. Кореляція: \t\t\t' + str(correl_inter[0]) + '\t\t' + str(correl) + '\t\t' + str(
                      correl_inter[1]) + '\t\t' +
                  str(round((1 - correl ** 2) / np.sqrt(n_x - 1), 4)) + '\n')

        xi_val = adequateReproduction(x, y)
        xi_val_quant = 0
        v_xi_val = b * b - 2
        if v_xi_val > 19 and v_xi_val < 24:
            xi_val_quant = 30.2
        elif v_xi_val > 25 and v_xi_val < 30:
            xi_val_quant = 37.35
        elif v_xi_val > 40 and v_xi_val < 50:
            xi_val_quant = 57.5
        elif v_xi_val > 60 and v_xi_val <= 90:
            xi_val_quant = 90.95

        if xi_val <= xi_val_quant:
            T3.insert(END,
                      '\n\u03C7^2: ' + str(xi_val) + ' <= ' + str(xi_val_quant) + ' --> Відтворення адекватне' + '\n')
        else:
            T3.insert(END,
                      '\n\u03C7^2: ' + str(xi_val) + ' > ' + str(xi_val_quant) + ' --> Відтворення неадекватне' + '\n')

        t_corr = t_corr_test(x, y)
        T3.insert(END, '\nПеревірка на значимість коефіцієнта кореляції: ' + str(t_corr))

        p_coff = correlationRelation(x, y)
        T3.insert(END, '\n\nКоефіцієнт кореляційного відношення: ' +
                  '\u03C1 = ' + str(p_coff))

        t_val = tValCorrRel(x, y)
        t_val_quant = 0
        v_t_val = n - 2
        if v_t_val > 10 and v_t_val < 60:
            t_val_quant = 2.115
        elif v_t_val > 60:
            t_val_quant = 1.96
        T3.insert(END, '\n\nt-статистика для перевірки зв’язку поміж випадковими величинами:\n')

        if t_val == None:
            pass
        elif t_val <= t_val_quant:
            T3.insert(END, 't: ' + str(t_val) + ' <= ' + str(
                t_val_quant) + ' --> Кореляційний зв’язок поміж випадковими величинами не наявний' + '\n')
        else:
            T3.insert(END, 't: ' + str(t_val) + ' > ' + str(
                t_val_quant) + ' --> Кореляційний зв’язок поміж випадковими величинами наявний' + '\n')

        f_val = fValCorrRel(x, y)
        v1_f_val = b - 1
        v2_f_val = n - b
        f_val_quant = 0

        if v1_f_val > 1 and v1_f_val < 5 and v2_f_val > 10 and v2_f_val < 100:
            f_val_quant = 3.665
        elif v1_f_val > 5 and v1_f_val <= 10 and v2_f_val > 60:
            f_val_quant = 2.08

        T3.insert(END, '\nf-статистика для перевірки значущості кореляційного відношення:\n')
        if f_val == None:
            pass
        elif f_val <= f_val_quant:
            T3.insert(END, 'f: ' + str(f_val) + ' <= ' + str(
                f_val_quant) + ' --> Кореляційне відношення не значуще' + '\n')
        else:
            T3.insert(END, 'f: ' + str(f_val) + ' > ' + str(
                f_val_quant) + ' --> Кореляційне відношення значуще' + '\n')

        spear_val = spearmanCoff(x, y)
        t_val_spearman = abs(round((spear_val * np.sqrt(n - 2)) / (np.sqrt(1 - spear_val ** 2)), 4))
        std_err_spearman = np.sqrt((1 - spear_val ** 2) / (n - 2))
        t_spearman = 0
        if (n - 2) > 10 and (n - 2) < 60:
            t_spearman = 2.115
        elif (n - 2) > 60:
            t_spearman = 1.96

        T3.insert(END, '\nКоефіцієнт Спірмена:\n' +
                  '\u03C4: ' + str(spear_val))
        if t_val_spearman <= t_spearman:
            T3.insert(END, '\nt: ' + str(t_val_spearman) + ' <= ' + str(
                t_spearman) + ' --> Оцінка коефіцієнту не значуща' + '\n')
        else:
            T3.insert(END, '\nt: ' + str(t_val_spearman) + ' > ' + str(
                t_spearman) + ' --> Оцінка коефіцієнту  значуща' + '\n')
        T3.insert(END,
                  'Інтервальне оцінювання коефіцієнту: \n' + f'[{round(spear_val - t_spearman * std_err_spearman, 4)}, {round(spear_val + t_spearman * std_err_spearman, 4)}]')

        kend_val = kendallCoff(x, y)
        u_val_kendall = round(((3 * kend_val) / np.sqrt(2 * (2 * n + 5))) * np.sqrt(n * (n - 1)), 4)
        std_err_u_val_kendall = np.sqrt((4 * n + 10) / (9 * (n ** 2 - n)))

        u_kendall = 1.96

        T3.insert(END, '\n\nКоефіцієнт Кендалла:\n' +
                  '\u03C4: ' + str(kend_val))
        if u_val_kendall <= u_kendall:
            T3.insert(END, '\nt: ' + str(u_val_kendall) + ' <= ' + str(
                u_kendall) + ' --> Оцінка коефіцієнту не значуща' + '\n')
        else:
            T3.insert(END, '\nt: ' + str(u_val_kendall) + ' > ' + str(
                u_kendall) + ' --> Оцінка коефіцієнту  значуща' + '\n')
        T3.insert(END,
                  'Інтервальне оцінювання коефіцієнту: \n' + f'[{round(kend_val - t_spearman * std_err_u_val_kendall, 4)}, {round(kend_val + t_spearman * std_err_u_val_kendall, 4)}]')

        n00, n01, n10, n11, m0, m1, n0, n1 = combTable2x2(x, y)
        T3.insert(END, '\n\nТаблиця сполучень 2х2:\n')
        T3.insert(END, 'Y \ X \t 0 \t 1\n')
        T3.insert(END, f'0 \t {n00} \t {n01} \t {n0}\n')
        T3.insert(END, f'1 \t {n10} \t {n11} \t {n1}\n')
        T3.insert(END, f' \t {m0} \t {m1} \t {n_x}\n')

        fecher_val = fecherInd(x, y)
        T3.insert(END, '\nІндекс Фехнера:\n' +
                  'I: ' + str(fecher_val))

        fi_val = fiCoff(x, y)
        T3.insert(END, '\n\nКоефіцієнт сполучень Фі:\n' +
                  'Fi: ' + str(fi_val))

        fi_val_test = round(n * fi_val ** 2, 4)
        T3.insert(END, '\nПеревірка значущості коефіцієнта Фі:\n')
        if fi_val_test >= 3.84:
            T3.insert(END, f'\u03C7^2: {fi_val_test} >= 3.84 --> оцінка коефіцієнта є значущою\n')
        else:
            T3.insert(END, f'\u03C7^2: {fi_val_test} < 3.84 --> оцінка коефіцієнта є не значущою\n')

        try:
            yule_Q, yule_Y, uQ, uY = yuleCoff(x, y)
            T3.insert(END, '\nКоефіцієнт зв\'язку Юла:\n' +
                      'Q, Y: ' + str(yule_Q) + '\t' + str(yule_Y) + '\n')
            if uQ <= 1.96 and uY <= 1.96:
                T3.insert(END, f'uQ: {uQ} <= 1.96 та uY: {uY} <= 1.96 --> коефіцієнти значущі\n')
            elif uQ > 1.96 and uY <= 1.96:
                T3.insert(END, f'uQ: {uQ} > 1.96 та uY: {uY} <= 1.96 --> коефіцієнти не значущі\n')
            elif uQ <= 1.96 and uY > 1.96:
                T3.insert(END, f'uQ: {uQ} <= 1.96 та uY: {uY} > 1.96 --> коефіцієнти не значущі\n')
            elif uQ > 1.96 and uY > 1.96:
                T3.insert(END, f'uQ: {uQ} > 1.96 та uY: {uY} > 1.96 --> коефіцієнти не значущі\n')
        except Exception as e:
            print("The error is: ", e)

        T3.insert(END, '\n\n\n\n\n\n')

        """Linear regression analysis"""
        T4.insert(END, 'Оцінка параметрів:\n')
        a_mnk, b_mnk = regrCoffMNK(x, y)
        T4.insert(END, f'a: {a_mnk}\n')
        T4.insert(END, f'b: {b_mnk}\n')
        T4.insert(END, 'Рівняння:\n')
        T4.insert(END, f'y = {b_mnk} * x + {a_mnk}\n\n')
        T4.insert(END, 'Метод Тейла:\n')
        a_teyl, b_teyl = regrCoffTeyla(x, y)
        T4.insert(END, f'a: {a_teyl}\n')
        T4.insert(END, f'b: {b_teyl}\n')
        T4.insert(END, 'Рівняння:\n')
        T4.insert(END, f'y = {b_teyl} * x + {a_teyl}\n\n')
        T4.insert(END, 'Дослідження точності оцінок параметрів:\n')
        a_inf, a_sup, b_inf, b_sup = inervalABLinRegr(x, y)
        T4.insert(END, f'a: [{a_inf}; {a_sup}]\n')
        T4.insert(END, f'b: [{b_inf}; {b_sup}]\n\n')
        T4.insert(END, 'Коефцієнт детермінації:\n')
        R2 = determinatCoff(x, y)
        T4.insert(END, f'R^2: {R2} %\n')
        T4.insert(END, 'Адекватність відтворення моделі регресії:\n')
        f_quant = 0
        if n - 1 <= 30:
            f_quant = 2.12
        elif n - 1 > 30 and n - 1 <= 60:
            f_quant = 1.68
        elif n - 1 > 60 and n - 1 <= 120:
            f_quant = 1.44
        elif n - 1 > 120:
            f_quant = 1
        f = adequateReproductionRegr(x, y)
        if f <= f_quant:
            T4.insert(END, f'f = {f} <= {f_quant} --> відтворення адекватне\n')
        else:
            T4.insert(END, f'f = {f} > {f_quant} --> відтворення неадекватне\n')

        """Nonlinear regression analysis"""
        T5.insert(END, 'Оцінка параметрів:\n')
        a_inf, a, a_sup, b_inf, b, b_sup, c_inf, c, c_sup = inervalNonLinRegr(x, y)
        T5.insert(END, f'a: {a}\n')
        T5.insert(END, f'b: {b}\n')
        T5.insert(END, f'c: {c}\n')
        T5.insert(END, 'Рівняння:\n')
        T5.insert(END, f'y = {a} + {b} * x + {c} * x**2\n\n')
        T5.insert(END, 'Довірче оцінювання параметрів:\n')
        T5.insert(END, f'a: [{a_inf}; {a_sup}]\n')
        T5.insert(END, f'b: [{b_inf}; {b_sup}]\n')
        T5.insert(END, f'c: [{c_inf}; {c_sup}]\n\n')

        t_quant = 0
        if n - 3 <= 30:
            t_quant = 2.05
        elif n - 3 > 30 and n - 3 <= 60:
            t_quant = 2.025
        elif n - 3 > 60 and n - 3 <= 120:
            t_quant = 1.99
        elif n - 3 > 120:
            t_quant = 1.96
        T5.insert(END, 'Оцінка точності та значущості параметрів:\n')
        t_a, t_b, t_c = tValsForParams(x, y)

        if t_a <= t_quant:
            T5.insert(END, f't_a = {t_a} <= {t_quant}  --> оцінка значуща\n')
        else:
            T5.insert(END, f't_a = {t_a} > {t_quant}  --> оцінка не значуща\n')

        if t_b <= t_quant:
            T5.insert(END, f't_b = {t_b} <= {t_quant}  --> оцінка значуща\n')
        else:
            T5.insert(END, f't_b = {t_a} > {t_quant}  --> оцінка не значуща\n')

        if t_c <= t_quant:
            T5.insert(END, f't_c = {t_c} <= {t_quant}  --> оцінка значуща\n')
        else:
            T5.insert(END, f't_c = {t_c} > {t_quant}  --> оцінка не значуща\n')

        T5.insert(END, '\nКоефцієнт детермінації:\n')
        R2_regr = determinatCoffNonLinRegr(x, y)
        T5.insert(END, f'R^2: {R2_regr} %\n')
        T5.insert(END, '\nАдекватність відтворення моделі регресії:\n')
        f_nonlin = adequateReproductionRegr(x, y)
        if f_nonlin <= f_quant:
            T5.insert(END, f'f = {f_nonlin} <= {f_quant} --> відтворення адекватне\n')
        else:
            T5.insert(END, f'f = {f_nonlin} > {f_quant} --> відтворення неадекватне\n')

        """Quasilinear regression analysis"""
        T6.insert(END, 'Оцінка параметрів:\n')
        a_inf, a, a_sup, b_inf, b, b_sup = inervalABQuaslinRegr(x, y)
        T6.insert(END, f'a: {round(a, 4)}\n')
        T6.insert(END, f'b: {round(b, 4)}\n')
        T6.insert(END, 'Рівняння:\n')
        T6.insert(END, f'y = {round(a, 4)} * exp({round(b, 4)} * x)\n\n')
        T6.insert(END, 'Довірче оцінювання параметрів:\n')
        T6.insert(END, f'a: [{round(a_inf, 4)}; {round(a_sup, 4)}]\n')
        T6.insert(END, f'b: [{round(b_inf, 4)}; {round(b_sup, 4)}]\n\n')

    except Exception as e:
        # By this way we can know about the type of error occurring
        print("The error is: ", e)


def n_m_values(x, y):
    user_input_n = simpledialog.askstring("n", "Введіть значення n")

    global n_val, m_val
    if user_input_n != "":
        n_val = int(user_input_n)

    user_input_m = simpledialog.askstring("m", "Введіть значення m")
    if user_input_m != "":
        m_val = int(user_input_m)

    T7 = Text(master=tab7, height=10, width=100)
    T7.grid(row=11, column=0)

    chi_val = pirsonCorr(x, y, n_val, m_val)
    tau_kendall = kendallCorr(x, y, n_val, m_val)

    tau_steward = stewardCorr(x, y, n_val, m_val)

    T7.insert(END, f"Коефіцієнт сполучень Пірсона: {chi_val}")
    if tau_kendall == 1:
        T7.insert(END, f"\nМіра зв’язку Кендалла: n != m")
    else:
        T7.insert(END, f"\nМіра зв’язку Кендалла: {tau_kendall}")
    if tau_steward == 1:
        T7.insert(END, f"\nМіра зв’язку Стюарда: n == m")
    else:
        T7.insert(END, f"\nМіра зв’язку Стюарда: {tau_steward}")



