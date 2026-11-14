from bisect import bisect_right
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
from typing import Sequence

from lab_common import check_table, plot_arguments, print_table, pyplot, save_figure


NODES = (0.0, 1.0, 2.0, 3.0, 4.0)
VALUES = (1.0, 0.86603, 0.5, 0.0, -0.5)
X_STAR = 1.5


def thomas(lower, diagonal, upper, rhs) -> list[float]:
    n = len(diagonal)
    if n == 0 or any(len(row) != n for row in (lower, upper, rhs)):
        raise ValueError("Неверные размеры трёхдиагональной системы")
    alpha, beta = [0.0] * n, [0.0] * n
    for i in range(n):
        denominator = diagonal[i] + (lower[i] * alpha[i - 1] if i else 0.0)
        if denominator == 0:
            raise ValueError("Нулевой знаменатель в методе прогонки")
        alpha[i] = -upper[i] / denominator
        beta[i] = (rhs[i] - (lower[i] * beta[i - 1] if i else 0.0)) / denominator
    result = [0.0] * n
    result[-1] = beta[-1]
    for i in range(n - 2, -1, -1):
        result[i] = alpha[i] * result[i + 1] + beta[i]
    return result


@dataclass(frozen=True)
class Segment:
    left: float
    right: float
    a: float
    b: float
    c: float
    d: float

    def value(self, x: float) -> float:
        t = x - self.left
        return self.a + t * (self.b + t * (self.c + t * self.d))

    def first(self, x: float) -> float:
        t = x - self.left
        return self.b + 2 * self.c * t + 3 * self.d * t**2

    def second(self, x: float) -> float:
        return 2 * self.c + 6 * self.d * (x - self.left)


def natural_spline(nodes: Sequence[float], values: Sequence[float]) -> list[Segment]:
    """S_i(x) = a_i + b_i*t + c_i*t² + d_i*t³; t = x - x_i."""
    check_table(nodes, values)
    n = len(nodes)
    h = [nodes[i + 1] - nodes[i] for i in range(n - 1)]
    slopes = [(values[i + 1] - values[i]) / h[i] for i in range(n - 1)]
    lower, diagonal, upper, rhs = ([0.0] * n for _ in range(4))
    diagonal[0] = diagonal[-1] = 1.0
    for i in range(1, n - 1):
        lower[i] = h[i - 1]
        diagonal[i] = 2 * (h[i - 1] + h[i])
        upper[i] = h[i]
        rhs[i] = 3 * (slopes[i] - slopes[i - 1])
    c = thomas(lower, diagonal, upper, rhs)
    return [Segment(
        nodes[i], nodes[i + 1], values[i],
        slopes[i] - h[i] * (2 * c[i] + c[i + 1]) / 3,
        c[i], (c[i + 1] - c[i]) / (3 * h[i]),
    ) for i in range(n - 1)]


def interval_index(segments: Sequence[Segment], x: float) -> int:
    if not segments or not isfinite(x) or not segments[0].left <= x <= segments[-1].right:
        raise ValueError("Точка должна принадлежать отрезку сплайна")
    return min(bisect_right([part.left for part in segments], x) - 1, len(segments) - 1)


def spline_value(segments: Sequence[Segment], x: float) -> float:
    return segments[interval_index(segments, x)].value(x)


def plot_solution(segments: Sequence[Segment], output_dir: Path, show: bool = False) -> None:
    plt = pyplot(show)
    fig, ax = plt.subplots(figsize=(9, 6), layout="constrained")
    for i, part in enumerate(segments):
        grid = [part.left + (part.right - part.left) * j / 100 for j in range(101)]
        ax.plot(grid, [part.value(x) for x in grid], linewidth=2,
                label=f"S_{i + 1}: [{part.left:g}, {part.right:g}]")
    ax.scatter(NODES, VALUES, color="black", zorder=5, label="Исходные узлы")
    value = spline_value(segments, X_STAR)
    ax.scatter([X_STAR], [value], color="tab:red", marker="*", s=140, zorder=6,
               label=f"(x*, S(x*)) = ({X_STAR:g}, {value:.6f})")
    ax.axvline(X_STAR, color="tab:red", linestyle=":", alpha=0.5)
    ax.set(title="3.2. Естественный кубический сплайн, вариант 2", xlabel="x", ylabel="S(x)")
    ax.grid(alpha=0.3)
    ax.legend()
    save_figure(fig, plt, output_dir, "spline", show)


def main(argv=None) -> None:
    args = plot_arguments("3.2. Кубический сплайн, вариант 2", argv)
    print("Лабораторная работа 3.2. Вариант 2")
    print("1. Исходные узлы и значения (таблица из задания, без замены на cos):")
    print_table(NODES, VALUES)
    segments = natural_spline(NODES, VALUES)
    print("\n2–3. Коэффициенты и интервалы. t = x - левая граница участка.")
    print("S_i(x) = a_i + b_i*t + c_i*t^2 + d_i*t^3; номера участков — от 1.")
    print(f"{'№':>2} {'Интервал':>12} {'a_i':>16} {'b_i':>16} {'c_i':>16} {'d_i':>16}")
    for i, s in enumerate(segments, 1):
        print(f"{i:2d} {f'[{s.left:g}, {s.right:g}]':>12} "
              f"{s.a:16.10f} {s.b:16.10f} {s.c:16.10f} {s.d:16.10f}")
        print(f"   S_{i}(x) = {s.a:.12g} + ({s.b:.12g})*(x-{s.left:g}) "
              f"+ ({s.c:.12g})*(x-{s.left:g})^2 + ({s.d:.12g})*(x-{s.left:g})^3")
    index = interval_index(segments, X_STAR)
    print(f"\n4. S({X_STAR}) = {spline_value(segments, X_STAR):.12f}")
    print(f"5. Номер интервала: {index + 1}, [{segments[index].left:g}, {segments[index].right:g}]")
    print("\n6. Проверка прохождения через узлы (допуск 1e-10):")
    for x, y in zip(NODES, VALUES):
        actual = spline_value(segments, x)
        print(f"x = {x:g}: S(x) = {actual:.12f}, y = {y:.12f}, "
              f"|S-y| = {abs(actual-y):.3e}, {'OK' if abs(actual-y) < 1e-10 else 'ОШИБКА'}")
    print("\n7. Естественные граничные условия S''(x_0) = S''(x_n) = 0:")
    print(f"S''({NODES[0]:g}) = {segments[0].second(NODES[0]):.12e}")
    print(f"S''({NODES[-1]:g}) = {segments[-1].second(NODES[-1]):.12e}")
    print("\nДополнительно: скачки S, S', S'' на внутренних стыках:")
    for left, right in zip(segments, segments[1:]):
        x = right.left
        print(f"x = {x:g}: {left.value(x)-right.value(x):.3e}, "
              f"{left.first(x)-right.first(x):.3e}, {left.second(x)-right.second(x):.3e}")
    if not args.no_plot:
        plot_solution(segments, args.output_dir, args.show)


if __name__ == "__main__":
    main()
