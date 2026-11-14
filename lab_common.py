"""Общие операции вывода и построения графиков; численные методы — в файлах работ."""

import argparse
import os
from math import isfinite
from pathlib import Path
from typing import Sequence


ROOT = Path(__file__).resolve().parent


def check_table(nodes: Sequence[float], values: Sequence[float], minimum: int = 2) -> None:
    if len(nodes) != len(values) or len(nodes) < minimum:
        raise ValueError(f"Нужны не менее {minimum} пар (x_i, y_i)")
    if not all(isfinite(v) for v in (*nodes, *values)):
        raise ValueError("Все узлы и значения должны быть конечными")
    if any(right <= left for left, right in zip(nodes, nodes[1:])):
        raise ValueError("Узлы должны строго возрастать")


def print_table(nodes: Sequence[float], values: Sequence[float]) -> None:
    print(f"{'i':>3} {'x_i':>16} {'y_i':>16}")
    for i, (x, y) in enumerate(zip(nodes, values)):
        print(f"{i:3d} {x:16.10f} {y:16.10f}")


def polynomial_text(coefficients: Sequence[float]) -> str:
    """Коэффициенты передаются в порядке a_0, a_1, ... ."""
    terms = []
    for i, value in enumerate(coefficients):
        basis = "" if i == 0 else ("*x" if i == 1 else f"*x^{i}")
        magnitude = f"{abs(value):.12g}{basis}"
        if not terms:
            terms.append(("-" if value < 0 else "") + magnitude)
        else:
            terms.append((" - " if value < 0 else " + ") + magnitude)
    return "".join(terms)


def plot_arguments(description: str, argv=None):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--no-plot", action="store_true", help="Только расчёты, без matplotlib")
    parser.add_argument("--show", action="store_true", help="Открыть график в окне после сохранения")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "plots", help="Каталог графиков")
    args = parser.parse_args(argv)
    if args.no_plot and args.show:
        parser.error("--no-plot и --show нельзя использовать вместе")
    return args


def pyplot(show: bool = False):
    cache = ROOT / ".cache" / "matplotlib"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache))
    import matplotlib

    if not show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def save_figure(fig, plt, output_dir: Path, name: str, show: bool = False) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{name}.png"
    fig.savefig(path, dpi=180)
    print(f"График сохранён: {path}")
    if show:
        plt.show()
    plt.close(fig)
