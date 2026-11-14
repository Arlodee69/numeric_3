from pathlib import Path
from typing import Sequence

from lab_common import check_table, plot_arguments, polynomial_text, print_table, pyplot, save_figure


NODES = (-1.0, 0.0, 1.0, 2.0, 3.0, 4.0)
VALUES = (0.86603, 1.0, 0.86603, 0.5, 0.0, -0.5)


def gaussian_solve(matrix: Sequence[Sequence[float]], rhs: Sequence[float]) -> list[float]:
    """Метод Гаусса с выбором главного элемента по столбцу."""
    n = len(rhs)
    if n == 0 or len(matrix) != n or any(len(row) != n for row in matrix):
        raise ValueError("Нужна квадратная матрица и согласованный столбец")
    augmented = [list(map(float, row)) + [float(value)] for row, value in zip(matrix, rhs)]
    for k in range(n):
        pivot = max(range(k, n), key=lambda i: abs(augmented[i][k]))
        if augmented[pivot][k] == 0:
            raise ValueError("Система вырождена")
        augmented[k], augmented[pivot] = augmented[pivot], augmented[k]
        for i in range(k + 1, n):
            factor = augmented[i][k] / augmented[k][k]
            augmented[i][k] = 0.0
            for j in range(k + 1, n + 1):
                augmented[i][j] -= factor * augmented[k][j]
    result = [0.0] * n
    for i in range(n - 1, -1, -1):
        result[i] = (augmented[i][n] - sum(augmented[i][j] * result[j]
                                          for j in range(i + 1, n))) / augmented[i][i]
    return result


def normal_equations(nodes, values, degree):
    check_table(nodes, values)
    if not isinstance(degree, int) or not 0 <= degree < len(nodes):
        raise ValueError("Степень должна быть целой, от 0 до количества узлов минус 1")
    moments = [sum(x**k for x in nodes) for k in range(2 * degree + 1)]
    matrix = [[moments[j + k] for k in range(degree + 1)] for j in range(degree + 1)]
    rhs = [sum(y * x**j for x, y in zip(nodes, values)) for j in range(degree + 1)]
    return matrix, rhs


def least_squares(nodes, values, degree) -> list[float]:
    matrix, rhs = normal_equations(nodes, values, degree)
    return gaussian_solve(matrix, rhs)


def polynomial(coefficients: Sequence[float], x: float) -> float:
    result = 0.0
    for coefficient in reversed(coefficients):
        result = coefficient + x * result
    return result


def errors(nodes, values, coefficients) -> tuple[float, float]:
    check_table(nodes, values)
    sse = sum((polynomial(coefficients, x) - y)**2 for x, y in zip(nodes, values))
    return sse, sse / len(nodes)


def plot_solution(fits, output_dir: Path, show: bool = False) -> None:
    plt = pyplot(show)
    fig, ax = plt.subplots(figsize=(9, 6), layout="constrained")
    grid = [NODES[0] + (NODES[-1] - NODES[0]) * i / 500 for i in range(501)]
    ax.scatter(NODES, VALUES, color="black", s=45, zorder=5, label="Исходные данные")
    for degree, coefficients in fits.items():
        ax.plot(grid, [polynomial(coefficients, x) for x in grid],
                label=f"P_{degree}(x), MSE = {errors(NODES, VALUES, coefficients)[1]:.6g}")
    ax.set(title="3.3. Метод наименьших квадратов, вариант 2", xlabel="x", ylabel="y")
    ax.grid(alpha=0.3)
    ax.legend()
    save_figure(fig, plt, output_dir, "least_squares", show)


def main(argv=None) -> None:
    args = plot_arguments("3.3. Метод наименьших квадратов, вариант 2", argv)
    print("Лабораторная работа 3.3. Вариант 2")
    print_table(NODES, VALUES)
    print("Коэффициенты: P_m(x) = a_0 + a_1*x + ... + a_m*x^m.")
    fits = {}
    for degree in (1, 2, 3):
        matrix, rhs = normal_equations(NODES, VALUES, degree)
        fits[degree] = gaussian_solve(matrix, rhs)
        print(f"\n{degree}. Многочлен степени {degree}. Нормальная система A*a = b:")
        for row, value in zip(matrix, rhs):
            print("  " + " ".join(f"{a:12.6g}" for a in row) + f" | {value:12.6g}")
        for i, value in enumerate(fits[degree]):
            print(f"a_{i} = {value:.12f}")
    print("\n4. Многочлены в явном виде:")
    for degree, coefficients in fits.items():
        print(f"P_{degree}(x) = {polynomial_text(coefficients)}")
    print("\n5. Ошибки по всем 6 исходным точкам:")
    print("SSE = sum((P(x_i)-y_i)^2); MSE = SSE / N, N = 6.")
    for degree, coefficients in fits.items():
        sse, mse = errors(NODES, VALUES, coefficients)
        print(f"SSE{degree} = {sse:.12e}; MSE{degree} = {mse:.12e}")
    best = min(fits, key=lambda degree: errors(NODES, VALUES, fits[degree])[1])
    print(f"\n6. Минимальная MSE на исходных данных у P_{best}(x), степень {best}.")
    print("Это сравнение на исходной выборке, оно не гарантирует лучшего прогноза вне неё.")
    if not args.no_plot:
        plot_solution(fits, args.output_dir, args.show)


if __name__ == "__main__":
    main()
