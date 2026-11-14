import argparse
from math import isfinite

from lab_common import check_table, print_table


NODES = (-1.0, 0.0, 1.0, 2.0, 3.0)
VALUES = (-0.5, 0.0, 0.5, 0.86603, 1.0)
X_STAR = 1.0


def select_nodes(nodes, values, x):
    check_table(nodes, values, minimum=3)
    if not isfinite(x) or not nodes[0] <= x <= nodes[-1]:
        raise ValueError("Точка должна принадлежать отрезку таблицы")
    indices = sorted(sorted(range(len(nodes)), key=lambda i: (abs(nodes[i] - x), i))[:3])
    return indices, [nodes[i] for i in indices], [values[i] for i in indices]


def derivatives(nodes, values, x) -> tuple[float, float]:
    check_table(nodes, values, minimum=3)
    if len(nodes) != 3:
        raise ValueError("Квадратичная формула использует ровно три узла")
    if not isfinite(x) or not nodes[0] <= x <= nodes[-1]:
        raise ValueError("Точка должна лежать между выбранными узлами")
    x0, x1, x2 = nodes
    y0, y1, y2 = values
    first_difference = (y1 - y0) / (x1 - x0)
    next_difference = (y2 - y1) / (x2 - x1)
    second_difference = (next_difference - first_difference) / (x2 - x0)
    # P_2 = y0 + f[x0,x1]*(x-x0) + f[x0,x1,x2]*(x-x0)*(x-x1).
    first = first_difference + second_difference * (2 * x - x0 - x1)
    second = 2 * second_difference
    return first, second


def main(argv=None) -> None:
    argparse.ArgumentParser(description="3.4. Численное дифференцирование, вариант 2").parse_args(argv)
    print("Лабораторная работа 3.4. Вариант 2")
    print("1. Исходная таблица из P_3_4.doc (значения не заменяются аналитической функцией):")
    print_table(NODES, VALUES)
    print(f"\n2. Точка вычисления x* = {X_STAR}")
    indices, nodes, values = select_nodes(NODES, VALUES, X_STAR)
    print(f"3. Использованные узлы: индексы {indices} (от 0), x = {nodes}, y = {values}")
    print("Метод: аналитическое дифференцирование локального квадратичного интерполянта.")
    first, second = derivatives(nodes, values, X_STAR)
    print(f"4. Первая производная f'(x*) ≈ {first:.12f}")
    print(f"5. Вторая производная f''(x*) ≈ {second:.12f}")
    print("При равном шаге h и центральном узле:")
    print("f'(x*) ≈ (y_right-y_left)/(2h)")
    print("f''(x*) ≈ (y_left-2*y_center+y_right)/h^2")


if __name__ == "__main__":
    main()
