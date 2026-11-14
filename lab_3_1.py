from math import cos, factorial, isfinite, pi, prod
from pathlib import Path
from typing import Sequence

from lab_common import plot_arguments, pyplot, save_figure


X_STAR = pi / 4
NODE_SETS = {
    "а": (0.0, pi / 6, pi / 3, pi / 2),
    "б": (0.0, pi / 6, 5 * pi / 12, pi / 2),
}


def validate_table(nodes: Sequence[float], values: Sequence[float]) -> None:
    if not nodes or len(nodes) != len(values):
        raise ValueError("Нужны непустые списки узлов и значений одинаковой длины")
    if not all(isfinite(v) for v in (*nodes, *values)):
        raise ValueError("Узлы и значения должны быть конечными числами")
    if len(set(nodes)) != len(nodes):
        raise ValueError("Узлы интерполяции должны быть различны")


def lagrange_weights(nodes: Sequence[float], values: Sequence[float]) -> list[float]:
    validate_table(nodes, values)
    return [
        y / prod(x_i - x_j for j, x_j in enumerate(nodes) if j != i)
        for i, (x_i, y) in enumerate(zip(nodes, values))
    ]


def lagrange(nodes: Sequence[float], values: Sequence[float], x: float) -> float:
    weights = lagrange_weights(nodes, values)
    return sum(
        weight * prod(x - x_j for j, x_j in enumerate(nodes) if j != i)
        for i, weight in enumerate(weights)
    )


def divided_differences(
    nodes: Sequence[float], values: Sequence[float]
) -> list[list[float]]:
    validate_table(nodes, values)
    n = len(nodes)
    table = [[float(value)] + [0.0] * (n - i - 1) for i, value in enumerate(values)]
    for order in range(1, n):
        for i in range(n - order):
            table[i][order] = (
                (table[i + 1][order - 1] - table[i][order - 1])
                / (nodes[i + order] - nodes[i])
            )
    return table


def newton(nodes: Sequence[float], coefficients: Sequence[float], x: float) -> float:
    validate_table(nodes, coefficients)
    result = coefficients[-1]
    for i in range(len(coefficients) - 2, -1, -1):
        result = coefficients[i] + (x - nodes[i]) * result
    return result

def factor(node: float) -> str:
    return "x" if node == 0 else f"(x - {node:.12g})"


def signed_sum(terms: Sequence[tuple[float, str]]) -> str:
    result = ""
    for coefficient, basis in terms:
        sign = "-" if coefficient < 0 else "+"
        term = f"{abs(coefficient):.12g}" + (f"*{basis}" if basis else "")
        result += (f" {sign} " if result else ("-" if coefficient < 0 else "")) + term
    return result


def print_solution(label: str, nodes: Sequence[float]) -> None:
    values = [cos(x) for x in nodes]
    table = divided_differences(nodes, values)
    coefficients = table[0]
    print(f"\nНабор {label})")
    print("1. Исходные узлы интерполяции и значения функции:")
    print(f"{'i':>3} {'x_i':>16} {'y_i = cos(x_i)':>16}")
    for i, (x_i, y_i) in enumerate(zip(nodes, values)):
        print(f"{i:3d} {x_i:16.12f} {y_i:16.12f}")

    print("\n2. Интерполяционный многочлен Лагранжа:")
    print("L_3(x) = " + signed_sum([
        (weight, "*".join(factor(x_j) for j, x_j in enumerate(nodes) if j != i))
        for i, weight in enumerate(lagrange_weights(nodes, values))
    ]))
    print("\n3. Таблица разделённых разностей:")
    print(f"{'i':>3} {'x_i':>16} {'f[x_i]':>16} {'порядок 1':>16} "
          f"{'порядок 2':>16} {'порядок 3':>16}")
    for i, row in enumerate(table):
        print(f"{i:3d} {nodes[i]:16.12f} " + " ".join(f"{v:16.12f}" for v in row))

    print("\n4. Интерполяционный многочлен Ньютона:")
    print("N_3(x) = " + signed_sum([
        (coefficient, "*".join(factor(nodes[j]) for j in range(i)))
        for i, coefficient in enumerate(coefficients)
    ]))

    exact = cos(X_STAR)
    l_value = lagrange(nodes, values, X_STAR)
    n_value = newton(nodes, coefficients, X_STAR)
    # f^(4)(x) = cos(x), max |cos(x)| = 1 на [0, pi/2].
    bound = abs(prod(X_STAR - x for x in nodes)) / factorial(4)
    print(f"\nx* = pi/4 = {X_STAR:.12f}")
    print(f"5. Значение Лагранжа L_3(x*) = {l_value:.12f}")
    print(f"6. Значение Ньютона N_3(x*)  = {n_value:.12f}")
    print(f"7. Точное значение f(x*) = cos(pi/4) = sqrt(2)/2 ≈ {exact:.12f}")
    print("8. Абсолютная погрешность интерполяции:")
    print(f"|cos(x*) - L_3(x*)|  = {abs(exact - l_value):.12e}")
    print(f"|cos(x*) - N_3(x*)|  = {abs(exact - n_value):.12e}")
    print("9. Разность между результатами методов:")
    print(f"L_3(x*) - N_3(x*)   = {l_value - n_value:.12e}")
    print(f"|L_3(x*) - N_3(x*)| = {abs(l_value - n_value):.12e}")


def plot_solutions(output_dir: Path, show: bool = False) -> None:
    """Общие графики для обоих наборов; сохраняются без открытия окна."""
    plt = pyplot(show)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), layout="constrained")
    for ax, (label, nodes) in zip(axes, NODE_SETS.items()):
        values = [cos(x) for x in nodes]
        coefficients = divided_differences(nodes, values)[0]
        grid = [nodes[0] + (nodes[-1] - nodes[0]) * i / 600 for i in range(601)]
        exact = [cos(x) for x in grid]
        l_values = [lagrange(nodes, values, x) for x in grid]
        n_values = [newton(nodes, coefficients, x) for x in grid]
        for target in (ax, ax.inset_axes([0.13, 0.12, 0.42, 0.32])):
            target.plot(grid, exact, color="black", linewidth=1.5, label="f(x) = cos(x)")
            target.plot(grid, l_values, color="tab:orange", linewidth=3,
                        label="Многочлен Лагранжа")
            target.plot(grid, n_values, color="tab:blue", linestyle="--", linewidth=1.5,
                        label="Многочлен Ньютона")
            target.axvline(X_STAR, color="tab:green", linestyle=":", label="x* = π/4")
            target.scatter([X_STAR], [cos(X_STAR)], color="tab:green", s=35, zorder=5,
                           label="(x*, f(x*))")
            target.scatter([X_STAR], [lagrange(nodes, values, X_STAR)], marker="x",
                           color="tab:red", s=45, zorder=6, label="(x*, P₃(x*))")
            target.grid(alpha=0.25)
            if target is not ax:
                target.set_xlim(X_STAR - 0.015, X_STAR + 0.015)
                target.set_ylim(cos(X_STAR) - 0.014, cos(X_STAR) + 0.012)
                target.set_title("Вблизи x*", fontsize=9)
                target.tick_params(labelsize=7)

        ax.scatter(nodes, values, color="black", s=40, zorder=5, label="Узлы (xᵢ, yᵢ)")
        ax.set_title(f"Набор {label})")
        ax.set_xlabel("x, радианы")
        ax.set_ylabel("y")
        ax.set_xticks([0, pi / 6, pi / 4, pi / 3, 5 * pi / 12, pi / 2],
                      ["0", "π/6", "π/4", "π/3", "5π/12", "π/2"])
        ax.legend(loc="upper right", fontsize=8)

    fig.suptitle("Интерполяция cos(x). Лабораторная 3.1, вариант 2")
    save_figure(fig, plt, output_dir, "interpolation", show)


def main(argv=None) -> None:
    args = plot_arguments("3.1. Интерполяция, вариант 2", argv)
    print("Лабораторная работа 3.1. Вариант 2")
    print("f(x) = cos(x); x* = pi/4. Все аргументы в радианах.")
    for label, nodes in NODE_SETS.items():
        print_solution(label, nodes)
    if not args.no_plot:
        plot_solutions(args.output_dir, args.show)


if __name__ == "__main__":
    main()
