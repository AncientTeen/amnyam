from tkinter import filedialog as fd
import seaborn as sns
from criteria import *
from mltDim import *
from paramMatching import *
from pca import *



"""hello amnyam"""
"""hello amnyam111"""

cls = None
# arr = None

sample_data = {}
samplesGroup = {}
sample_checkbuttons = []
y_sample = 1
regBound = [1, 10]
fi = None
eigenvectors = None


def prevent_menu_close(event):
    pass


def activate_sample(sample_name):
    sample_data[sample_name]['var'].set(1)


def openFile():
    root.filename = fd.askopenfilename(initialdir="/", title="Select file", filetypes=[('All Files', '*.*'),
                                                                                       ('Python Files', '*.py'),
                                                                                       ('Text Document', '*.txt'),
                                                                                       ('CSV files', "*.csv")])
    global array
    t = []
    if root.filename.split('.')[1] == 'txt':
        array = np.loadtxt(root.filename, dtype='float')

        # ic(array)
        # ic(array.shape)
        # ic(len(array.shape))
        # ic(len(array))

        if len(array.shape) == 1:
            t.append(array)
        elif len(array.shape) > 1:
            for i in range(len(array[0])):
                t_buff = []
                for j in range(len(array)):
                    t_buff.append(array[j][i])
                t.append(t_buff)

    for i in range(len(t)):
        sample_num = len(sample_data) + 1
        sample_name = f"Вибірка {sample_num}"
        sample_var = tkinter.IntVar()
        arr = []
        for j in range(len(t[0])):
            arr.append(t[i][j])
        sample_data[sample_name] = {"data": arr, "var": sample_var}

        checkbutton = Checkbutton(text=sample_name, variable=sample_var)
        sample_checkbuttons.append(checkbutton)

        sample_menu.add_checkbutton(label=sample_name, variable=sample_var)


def createHist(sample, master):
    fig, ax = plt.subplots(figsize=(5, 4), dpi=100)

    if cls:
        b = cls
    else:
        if len(sample) < 100:
            b = round((len(sample) ** (1 / 2)))
            if b % 2 == 0:
                b -= 1
        else:
            b = round((len(sample) ** (1 / 3)))
            if b % 2 == 0:
                b -= 1

    plt.grid(color='grey', linestyle='--', linewidth=0.5)
    plt.xlabel('Варіанти')
    plt.ylabel('Частоти')

    plt.title('Відносні частоти')

    plt.hist(sample, bins=b, edgecolor="black", color='blue', weights=np.ones_like(sample) / len(sample))

    hist = FigureCanvasTkAgg(fig, master=master)
    hist.get_tk_widget().place(x=50, y=40)
    toolbar = NavigationToolbar2Tk(hist, master, pack_toolbar=False)
    toolbar.update()
    toolbar.place(x=150, y=450)


def create_distribution_function(sample, master):
    arr = np.sort(sample)

    fig, ax = plt.subplots(figsize=(5, 4), dpi=100)

    plt.grid(color='grey', linestyle='--', linewidth=0.5)

    n = len(arr)

    t = 0
    if n - 1 < 120:
        if n - 1 == 69:
            t = 1.995
        elif n - 1 == 24:
            t = 2.06
    else:
        t = 1.96

    if cls:
        b = cls
    else:
        if n < 100:
            b = round((n ** (1 / 2)))
        else:
            b = round((n ** (1 / 3)))

    s_y = np.arange(1, n + 1) / n
    ax.scatter(x=arr, y=s_y, s=5)
    sns.histplot(arr, element="step", fill=False,
                 cumulative=True, stat="density", common_norm=False, bins=b, color='red')

    plt.title('Функція розподілу')
    plt.xlabel('')
    plt.ylabel('')

    distr_func = FigureCanvasTkAgg(fig, master=master)
    distr_func.get_tk_widget().place(x=800, y=40)
    toolbar = NavigationToolbar2Tk(distr_func, master, pack_toolbar=False)
    toolbar.update()
    toolbar.place(x=900, y=450)


def outputData(sample, master):
    tabControl = Notebook(master)

    tab1 = Frame(tabControl)
    tab2 = Frame(tabControl)
    tab3 = Frame(tabControl)

    tabControl.add(tab1, text='Об\'єкти')
    tabControl.add(tab2, text='Протокол')
    tabControl.add(tab3, text='Оцінка критеріїв')

    tabControl.place(x=10, y=500)

    T1 = Text(master=tab1, height=10, width=150)

    T2 = Text(master=tab2, height=10, width=150)

    T3 = Text(master=tab3, height=10, width=150)
    T1.pack()
    T2.pack()
    T3.pack()

    T1.delete('1.0', END)
    T2.delete('1.0', END)

    arr = np.sort(sample)

    T1.insert(END, f'data: {arr}')

    n = len(arr)
    T2.insert(END, 'Характеристика\t\t\t' + 'INF\t\t' + 'Значення\t\t' + 'SUP\t\t' + 'SKV\n')

    avr = average(arr)
    sq = std_err(arr, avr)
    avrIntr = confInterAvr(arr)
    T2.insert(END,
              'Середнє значення: \t\t\t' + str(avrIntr[0]) + '\t\t' + str(avr) + '\t\t' + str(avrIntr[1]) + '\t\t' +
              str(round(sq / (len(arr) ** (1 / 2)), 4)) + '\n')

    md = medium(arr)
    T2.insert(END, 'Медіана: \t\t\t' + str(md) + '\t\t' + str(md) + '\t\t' + str(md) + '\t\t' + '0.0000\n')

    avrSq = std_err(arr, avr)
    sqIntr = confInterStdErr(arr)
    T2.insert(END,
              'Сер. квадратичне: \t\t\t' + str(sqIntr[0]) + '\t\t' + str(avrSq) + '\t\t' + str(sqIntr[1]) + '\t\t' +
              str(round(sq * (2 / (n - 1)) ** (1 / 4), 4)) + '\n')

    assmCf = assymCoef(arr, avr)
    assmIntrCof = confInterAssym(arr)
    T2.insert(END, 'Коефіцієнт асиметрії: \t\t\t' + str(assmIntrCof[0]) + '\t\t' + str(assmCf) + '\t\t' + str(
        assmIntrCof[1]) + '\t\t' +
              str(round((6 * (n - 2) / ((n + 1) * (n + 3))) ** (1 / 2), 4)) + '\n')

    exCf = excessCoef(arr, avr)
    exIntrCof = confInterExcess(arr)
    T2.insert(END, 'Коефіцієнт ексцесу: \t\t\t' + str(exIntrCof[0]) + '\t\t' + str(exCf) + '\t\t' + str(
        exIntrCof[1]) + '\t\t' +
              str(round((24 * n * (n - 1) ** 2 / ((n - 3) * (n - 2) * (n + 3) * (n + 5))) ** (1 / 2), 4)) + '\n')

    shftSq = 0

    for i in range(n):
        shftSq += arr[i] ** 2 - avr ** 2
    shftSq = round((shftSq / n) ** (1 / 2), 4)

    shftExCf = 0
    for i in range(n):
        shftExCf += (arr[i] - avr) ** 4
    shftExCf = shftExCf / (n * (shftSq ** 4))

    cntrExCf = contrExcessCoef(exCf)
    cExIntrCof = confInterContrEx(arr)

    T2.insert(END, 'Коефіцієнт контрексцесу: \t\t\t' + str(cExIntrCof[0]) + '\t\t' + str(cntrExCf) + '\t\t' + str(
        cExIntrCof[1]) + '\t\t' +
              str(round(((abs(shftExCf) / (29 * n)) ** (1 / 2)) * ((abs(shftExCf ** 2 - 1)) ** (3 / 4)), 4)) + '\n')

    prsCf = pirsonCoef(avrSq, avr)
    variatIntrCof = confInterVariation(arr)
    if prsCf is None or prsCf < 10 or prsCf > 10:
        T2.insert(END, 'Коефіцієнт Варіації: \t\t\t\t\t' + str(prsCf) + '\n')
    else:
        T2.insert(END, 'Коефіцієнт Варіації: \t\t\t' + str(variatIntrCof[0]) + '\t\t' + str(prsCf) + '\t\t' + str(
            variatIntrCof[1]) + '\t\t' +
                  str(round(prsCf * (((1 + 2 * prsCf) / (2 * n)) ** (1 / 2)), 4)) + '\n')

    trcnAvr = truncatedAverage(arr)
    T2.insert(END, 'Усічене середнє: \t\t\t\t\t' + str(trcnAvr) + '\n')

    mdWlsh = mediumWalsh(arr)
    T2.insert(END, 'Медіана Уолша: \t\t\t\t\t' + str(mdWlsh) + '\n')

    mdAbsMss = mediumAbsMiss(arr, md)
    T2.insert(END, 'Медіана абс. відхилень: \t\t\t\t\t' + str(mdAbsMss) + '\n')

    nonParamCfVar = nonParamCoefVar(mdAbsMss, md)
    T2.insert(END, 'Непарам. коеф. варіації: \t\t\t\t\t' + str(nonParamCfVar) + '\n')

    T2.insert(END, '-----------------------------\n')
    T2.insert(END, 'Квантилі : \n')
    T2.insert(END, '0.05: \t' + str(round(np.quantile(arr, 0.05), 3)) + '\n')
    T2.insert(END, '0.1: \t' + str(round(np.quantile(arr, 0.1), 3)) + '\n')
    T2.insert(END, '0.25: \t' + str(round(np.quantile(arr, 0.25), 3)) + '\n')
    T2.insert(END, '0.5: \t' + str(round(np.quantile(arr, 0.5), 3)) + '\n')
    T2.insert(END, '0.75: \t' + str(round(np.quantile(arr, 0.75), 3)) + '\n')
    T2.insert(END, '0.9: \t' + str(round(np.quantile(arr, 0.9), 3)) + '\n')
    T2.insert(END, '0.95: \t' + str(round(np.quantile(arr, 0.95), 3)) + '\n')


def regrSample():
    while True:
        user_input1 = simpledialog.askstring("r", "Введіть ознаку для регресії")
        if user_input1 != "":
            global y_sample
            y_sample = int(user_input1)
        if y_sample > len(sample_data):
            continue
        else:
            break


def regrSlice():
    while True:
        user_input1 = simpledialog.askstring("low", "Введіть нижню межу")
        user_input2 = simpledialog.askstring("up", "Введіть верхню")
        if user_input1 != "":
            low = int(user_input1)
            regBound[0] = low
        if user_input2 != "":
            up = int(user_input2)
            regBound[1] = up

        if low > up or low < 1 or up > len(sample_data["Вибірка 1"]["data"]):

            continue
        else:
            break


def scatterplotMatrix(sample_data):
    s_n = []
    for i in range(1, len(sample_data) + 1):
        str = f"Вибірка {i}"
        if sample_data[str]['var'].get() == 1:
            s_n.append(str)

    buff = np.array([sample_data[s_n[i]]["data"] for i in range(len(s_n))])

    if len(buff[0]) < 100:
        b = round((len(buff[0]) ** (1 / 2)))
        if b % 2 == 0:
            b -= 1
    else:
        b = round((len(buff[0]) ** (1 / 3)))
        if b % 2 == 0:
            b -= 1

    fig5, ax5 = plt.subplots(len(buff), len(buff), figsize=(16, 9))
    for i in range(len(buff)):
        for j in range(len(buff)):
            if i == j:
                ax5[i][j].hist(buff[i], bins=b, edgecolor="black", color='blue',
                               weights=np.ones_like(buff[i]) / len(buff[i]))
            elif i > j:
                ax5[i][j].scatter(x=buff[i], y=buff[j], c='red')
            elif j > i:
                ax5[i][j].text(0.5, 0.5, f'Corr: {linearCorrelation(buff[i], buff[j]):.4f}', ha='center', va='center',
                               fontsize=14)
    plt.tight_layout()
    fig5.show()


def showSample():
    dim1_check = var_1d.get()
    dim2_check = var_2d.get()
    dim3_check = var_3d.get()
    if dim1_check == '1':
        s_n = ''
        for i in range(1, len(sample_data) + 1):
            str = f"Вибірка {i}"
            if sample_data[str]['var'].get() == 1:
                s_n = str
        createHist(sample_data[s_n]["data"], frame_1d)
        create_distribution_function(sample_data[s_n]["data"], frame_1d)
        outputData(sample_data[s_n]["data"], frame_1d)

    elif dim2_check == '1':
        s_n = []
        for i in range(1, len(sample_data) + 1):
            str = f"Вибірка {i}"
            if sample_data[str]['var'].get() == 1:
                s_n.append(str)

        SndDimDens(sample_data[s_n[0]]["data"], sample_data[s_n[1]]["data"], frame_2d)
        variation_series(sample_data[s_n[0]]["data"], sample_data[s_n[1]]["data"], frame_2d)
        outputDataScdDim(sample_data[s_n[0]]["data"], sample_data[s_n[1]]["data"], frame_2d)

    elif dim3_check == '1':
        s_n = []
        for i in range(1, len(sample_data) + 1):
            str = f"Вибірка {i}"
            if sample_data[str]['var'].get() == 1:
                s_n.append(str)
        # ic(regBound)
        e, Y = outputDataMlt(sample_data, s_n, frame_3d, y_sample, regBound)

        visualization(sample_data, s_n, frame_3d, e, Y)


def expParams():
    user_input1 = simpledialog.askstring("N", "Введіть розмірність вибірки")
    user_input2 = simpledialog.askstring("Lambd", "Введіть параметр \u03BB")
    n = 0
    lambd = 0
    if user_input1 != "":
        n = int(user_input1)
    if user_input1 != "":
        lambd = float(user_input2)

    sampleGenExp(n, lambd)


def sampleGenExp(n, lambd):
    alfa = np.random.default_rng().uniform(0, 1, n)
    sample = (1 / lambd) * np.log(1 / (1 - alfa))
    # global arr
    arr = np.sort(sample)
    avr = average(arr)
    lambd_teta = 1 / avr
    lambd_sq = (lambd ** 2) / n

    t = 0
    l = 0
    if int(n) < 60:
        t = 1.16
        l = 0.3
    elif int(n) >= 60 and int(n) <= 200:
        t = 1.345
        l = 0.1
    elif int(n) >= 200 and int(n) <= 400:
        t = 1.44
        l = 0.075
    elif int(n) >= 400 and int(n) <= 1000:
        t = 1.645
        l = 0.05
    elif int(n) > 1000:
        t = 2.33
        l = 0.01

    tTest = (lambd - lambd_teta) / lambd_sq

    sample_num = len(sample_data) + 1
    sample_name = f"Вибірка {sample_num}"
    sample_var = tkinter.IntVar()
    sample_data[sample_name] = {"data": arr, "var": sample_var}

    sample_menu.add_checkbutton(label=sample_name, variable=sample_var)

    createHist()
    create_distribution_function()
    outputDataScdDim()


def switch_to_tab(tab_index):
    notebook.select(tab_index)


def switch_to_1d():
    switch_to_tab(0)


def switch_to_2d():
    switch_to_tab(1)


def switch_to_3d():
    switch_to_tab(2)


def addSampleGroup(sample_data):
    s_n = []
    for i in range(1, len(sample_data) + 1):
        str = f"Вибірка {i}"
        if sample_data[str]['var'].get() == 1:
            s_n.append(str)

    buff = np.array([sample_data[s_n[i]]["data"] for i in range(len(s_n))])

    group_num = len(samplesGroup) + 1
    group_name = f"Набір {group_num}"
    group_var = tkinter.IntVar()

    samplesGroup[group_name] = {"samples_group": buff, "var": group_var}
    sample_group.add_checkbutton(label=group_name, variable=group_var)


root = Tk()

root.geometry("1400x800")

notebook = Notebook(root)
notebook.pack(fill='both', expand=True)

# Frame for 1D tab
frame_1d = Frame(notebook)
frame_1d.configure(bg="pink")
notebook.add(frame_1d)
#
# # Frame for 2D tab
frame_2d = Frame(notebook)
frame_2d.configure(bg="pink")
notebook.add(frame_2d)

# Frame for 3D tab
frame_3d = Frame(notebook)
frame_3d.configure(bg="pink")
notebook.add(frame_3d)

var_1d = StringVar(value="0")
var_2d = StringVar(value="0")
var_3d = StringVar(value="0")

ch1 = Checkbutton(root, text="1D", variable=var_1d, onvalue="1", offvalue="0")
ch1.place(x=10, y=10)
ch2 = Checkbutton(root, text="2D", variable=var_2d, onvalue="1", offvalue="0")
ch2.place(x=70, y=10)
ch3 = Checkbutton(root, text="3D", variable=var_3d, onvalue="1", offvalue="0")
ch3.place(x=130, y=10)

"""menu widgets"""
menubar = Menu(root)
filemenu = Menu(menubar, tearoff=0)
menubar.add_cascade(label="Меню", menu=filemenu)
filemenu.add_command(label="Відкрити файл", command=openFile)

# filemenu.add_command(label="Класи", command=classes)

sample_menu = Menu(menubar, tearoff=0)

dimensions_menu = Menu(menubar, tearoff=0)
filemenu.add_cascade(label="Виміри", menu=dimensions_menu)
dimensions_menu.add_command(label="1D", command=switch_to_1d)
dimensions_menu.add_command(label="2D", command=switch_to_2d)
dimensions_menu.add_command(label="3D", command=switch_to_3d)

filemenu.add_cascade(label="Вибірки", menu=sample_menu)

filemenu.add_command(label="Відобразити", command=showSample)
filemenu.add_command(label="Вийти", command=root.quit)

"""tools widgets"""
tools = Menu(menubar, tearoff=0)
menubar.add_cascade(label="Інструментарій", menu=tools)

transformation_menu = Menu(menubar, tearoff=0)
tools.add_cascade(label="Перетворення", menu=transformation_menu)
transformation_menu.add_command(label="Логарифмувати", command=lambda: logarimization(sample_data))
transformation_menu.add_command(label="Стандартизувати", command=lambda: standartization(sample_data))
transformation_menu.add_command(label="Вилучення аномальних значень", command=lambda: removeAnomals(sample_data))

"""PCA"""
pca_menu = Menu(menubar, tearoff=0)
transformation_menu.add_cascade(label='МГК', menu=pca_menu)
pca_menu.add_command(label="Перехід до незалежних значень 2n",
                     command=lambda: pca_two_frwd(sample_data))
pca_menu.add_command(label="Зворотній перехід 2n", command=lambda: pca_two_bck(sample_data))
pca_menu.add_command(label="Перехід по w-ознакам", command=lambda: pcaFnDim(sample_data, sample_menu, sample_checkbuttons))
pca_menu.add_command(label="Зворотній перехід по w-ознакам", command=lambda: pcaBnDim(sample_data, sample_menu, sample_checkbuttons))

"""criteria"""
criteria_menu = Menu(menubar, tearoff=0)
tools.add_cascade(label="Критерії однорідності", menu=criteria_menu)
# sampleMod_menu = Menu(menubar, tearoff=0)
# sampleMod_menu.add_command(label="Експоненціальний", command=expParams)
# tools.add_cascade(label="Mоделювання вибірок", menu=sampleMod_menu)
depend_menu = Menu(criteria_menu, tearoff=0)
independ_menu = Menu(criteria_menu, tearoff=0)

criteria_menu.add_cascade(label="Залежні", menu=depend_menu)
criteria_menu.add_cascade(label="Незалежні", menu=independ_menu)

depend_menu.add_command(label="Перевірка збігу середніх", command=lambda: avr_match_dep(sample_data, T3))
depend_menu.add_command(label="Критерій знаків", command=lambda: sign_test(sample_data, T3))
depend_menu.add_command(label="Q Критерій", command=lambda: q_test(sample_data, T3))
depend_menu.add_command(label="Критерій Аббе", command=lambda: abbe_test(sample_data, T3))

independ_menu.add_command(label="Перевірка збігу середніх", command=lambda: avr_match_indep(sample_data, T3))
independ_menu.add_command(label="Перевірка збігу дисперсій", command=lambda: f_test(sample_data, T3))
independ_menu.add_command(label="Критерій Бартлетта", command=lambda: Bartlett_test(sample_data, T3))
independ_menu.add_command(label="Однофакторний дисперсійний аналіз", command=lambda: anova_func(sample_data, T3))
independ_menu.add_command(label="Критерій Смирнова-Колмогорова",
                          command=lambda: Kolmogorov_Smirnov_test(sample_data, T3))
independ_menu.add_command(label="Критерій Вілкоксона", command=lambda: wilcoxon_test(sample_data, T3))
independ_menu.add_command(label="U Критерій", command=lambda: u_test(sample_data, T3))
independ_menu.add_command(label="Різниця середніх рангів", command=lambda: diff_mean_rank(sample_data, T3))
independ_menu.add_command(label="H Критерій", command=lambda: h_test(sample_data, T3))

sample_group = Menu(menubar, tearoff=0)
tools.add_cascade(label='Набір вибірок', menu=sample_group)
tools.add_command(label='Додати набір', command=lambda: addSampleGroup(sample_data))
tools.add_command(label='Перевірка на збіг', command=lambda: checkMatching_dc(samplesGroup))
tools.add_command(label='Ознака для Рег.аналіз', command=regrSample)
tools.add_command(label='Зріз для регресії', command=regrSlice)
tools.add_command(label='Матриця діаграм розкиду', command=lambda: scatterplotMatrix(sample_data))

root.config(menu=menubar)

root.configure(bg="pink")

root.mainloop()
