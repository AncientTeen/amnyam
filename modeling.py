from tkinter import *
from tkinter import filedialog as fd, simpledialog

def normSampParams():
    user_input1 = simpledialog.askstring("x mean", "Введіть x середнє")
    user_input2 = simpledialog.askstring("y mean", "Введіть y середнє")
    user_input3 = simpledialog.askstring("std err x", "Введіть середньоквадратичне відхилення х")
    user_input4 = simpledialog.askstring("std err y", "Введіть середньоквадратичне відхилення y")
    user_input5 = simpledialog.askstring("correlation coeeficient(r)", "Введіть коефіцієнт r")
    user_input6 = simpledialog.askstring("N", "Введіть обсяг вибірки")
    user_input7 = simpledialog.askstring("File name", "Введіть назву файлу для збереження")
    x_mean = 0
    y_mean = 0
    x_std_err = 0
    y_std_err = 0
    r = 0
    N = 0
    file_name = ''
    if user_input1 != "":
        x_mean = float(user_input1)
    if user_input2 != "":
        y_mean = float(user_input2)
    if user_input3 != "":
        x_std_err = float(user_input3)
    if user_input4 != "":
        y_std_err = float(user_input4)
    if user_input5 != "":
        r = float(user_input5)
    if user_input6 != "":
        N = int(user_input6)
    if user_input7 != "":
        file_name = str(user_input7)

    np.random.seed(42)  # Set seed for reproducibility

    # Covariance matrix
    cov_matrix = [[x_std_err ** 2, r * x_std_err * y_std_err],
                  [r * x_std_err * y_std_err, y_std_err ** 2]]

    # Generate 2D normal distributed data
    x, y = np.random.multivariate_normal([x_mean, y_mean], cov_matrix, N).T

    with open(file_name, 'w') as file:
        for xi, yi in zip(x, y):
            file.write(f"{xi}\t{yi}\n")


def nonLinearSampParams():
    user_input1 = simpledialog.askstring("a", "Введіть параметр a")
    user_input2 = simpledialog.askstring("b", "Введіть параметр b")
    user_input3 = simpledialog.askstring("c", "Введіть параметр c")
    user_input4 = simpledialog.askstring("q", "Введіть параметр q")
    user_input5 = simpledialog.askstring("N", "Введіть обсяг вибірки")
    user_input6 = simpledialog.askstring("File name", "Введіть назву файлу для збереження")

    a = 0
    b = 0
    c = 0
    q = 0
    N = 0
    file_name = ''

    if user_input1 != "":
        a = float(user_input1)
    if user_input2 != "":
        b = float(user_input2)
    if user_input3 != "":
        c = float(user_input3)
    if user_input4 != "":
        q = float(user_input4)
    if user_input5 != "":
        N = int(user_input5)
    if user_input6 != "":
        file_name = str(user_input6)

    np.random.seed(42)  # Set seed for reproducibility

    # Generate nonlinear data based on a quadratic polynomial
    x = np.linspace(-10, 10, N)
    y = c * x ** 2 + b * x + a + np.random.normal(0, q, N)

    with open(file_name, 'w') as file:
        for xi, yi in zip(x, y):
            file.write(f"{xi}\t{yi}\n")


def quasiLinearSampParams():
    user_input1 = simpledialog.askstring("a", "Введіть параметр a")
    user_input2 = simpledialog.askstring("b", "Введіть параметр b")
    user_input4 = simpledialog.askstring("q", "Введіть параметр q")
    user_input7 = simpledialog.askstring("X min", "Введіть мінімальне значення Х")
    user_input8 = simpledialog.askstring("X max", "Введіть максимальне значення Х")
    user_input5 = simpledialog.askstring("N", "Введіть обсяг вибірки")
    user_input6 = simpledialog.askstring("File name", "Введіть назву файлу для збереження")

    a = float(user_input1) if user_input1 else 0
    b = float(user_input2) if user_input2 else 0
    q = float(user_input4) if user_input4 else 0
    N = int(user_input5) if user_input5 else 0
    x_min = float(user_input7) if user_input7 else 0
    x_max = float(user_input8) if user_input8 else 0
    file_name = str(user_input6) if user_input6 else ''

    np.random.seed(42)  # Set seed for reproducibility

    # Generate quasi-linear data based on a linear relationship with some curvature
    x = np.linspace(x_min, x_max, N)

    noise = np.random.normal(0, q, N)
    y = [(a * np.exp(b * x[i]) + noise[i]) for i in range(N)]

    with open(file_name, 'w') as file:
        for xi, yi in zip(x, y):
            file.write(f"{xi}\t{yi}\n")

