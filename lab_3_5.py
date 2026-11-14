"""3.5. Составные квадратуры и метод Рунге–Ромберга, вариант 2."""

import argparse
from math import isclose, isfinite, log
from typing import Callable


A, B = 0.0, 4.0
H1, H2 = 1.0, 0.5
METHODS = {"Левые прямоугольники": 1, "Правые прямоугольники": 1,
           "Средние прямоугольники": 2, "Трапеции": 2, "Симпсон": 4}


def integrand(x: float) -> float:
    return x / (3 * x + 4)**2


def antiderivative(x: float) -> float:
    # u = 3x+4: integral = (ln|u| + 4/u) / 9.
    return (log(abs(3 * x + 4)) + 4 / (3 * x + 4)) / 9


def subinterval_count(a: float, b: float, h: float) -> int:
    if not all(isfinite(v) for v in (a, b, h)) or h <= 0 or b <= a:
        raise ValueError("Нужны конечные a < b и положительный шаг")
    quotient = (b - a) / h
    count = round(quotient)
    if count < 1 or not isclose(quotient, count, rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError("Шаг должен делить длину отрезка на целое число частей")
    return count


def integrate(f: Callable[[float], float], a: float, b: float, h: float, method: str) -> float:
    n = subinterval_count(a, b, h)
    if method not in METHODS:
        raise ValueError("Неизвестный метод интегрирования")
    if method == "Левые прямоугольники":
        return h * sum(f(a + i * h) for i in range(n))
    if method == "Правые прямоугольники":
        return h * sum(f(a + (i + 1) * h) for i in range(n))
    if method == "Средние прямоугольники":
        return h * sum(f(a + (i + 0.5) * h) for i in range(n))
    if method == "Трапеции":
        return h * ((f(a) + f(b)) / 2 + sum(f(a + i * h) for i in range(1, n)))
    if n % 2:
        raise ValueError("Для составного метода Симпсона число частей должно быть чётным")
    return h / 3 * (f(a) + f(b)
                    + 4 * sum(f(a + i * h) for i in range(1, n, 2))
                    + 2 * sum(f(a + i * h) for i in range(2, n, 2)))


def runge_romberg(coarse: float, fine: float, h1: float, h2: float, p: int):
    if not all(isfinite(v) for v in (coarse, fine, h1, h2, p)) or not h1 > h2 > 0 or p <= 0:
        raise ValueError("Нужны конечные значения, h1 > h2 > 0 и p > 0")
    correction = (fine - coarse) / ((h1 / h2)**p - 1)
    return fine + correction, abs(correction)


def main(argv=None) -> None:
    argparse.ArgumentParser(description="3.5. Численное интегрирование, вариант 2").parse_args(argv)
    print("Лабораторная работа 3.5. Вариант 2")
    print(f"I = integral x/(3x+4)^2 dx от {A:g} до {B:g}, h1 = {H1:g}, h2 = {H2:g}")
    exact = antiderivative(B) - antiderivative(A)
    print("Первообразная F(x) = (ln|3x+4| + 4/(3x+4))/9")
    print(f"Точное значение I = (ln(4)-3/4)/9 ≈ {exact:.12f}")
    print("\nI_RR = I_h2 + (I_h2-I_h1)/((h1/h2)^p-1)")
    print("Оценка ошибки h2 = |I_RR-I_h2|. Это асимптотическая оценка, не строгая граница.")
    print(f"\n{'Метод':24} {'I(h1)':>16} {'I(h2)':>16} {'p':>2} {'I_RR':>16} {'Оценка h2':>14}")
    results = {}
    for method, p in METHODS.items():
        coarse = integrate(integrand, A, B, H1, method)
        fine = integrate(integrand, A, B, H2, method)
        refined, estimate = runge_romberg(coarse, fine, H1, H2, p)
        results[method] = (coarse, fine, refined)
        print(f"{method:24} {coarse:16.12f} {fine:16.12f} {p:2d} "
              f"{refined:16.12f} {estimate:14.6e}")
    print(f"\n{'Метод':24} {'|I(h1)-I|':>16} {'|I(h2)-I|':>16} {'|I_RR-I|':>16}")
    for method, values in results.items():
        print(f"{method:24} " + " ".join(f"{abs(value-exact):16.8e}" for value in values))
    best = min(results, key=lambda method: abs(results[method][1] - exact))
    best_refined = min(results, key=lambda method: abs(results[method][2] - exact))
    print(f"\nНаиболее точный исходный результат при h2: {best}.")
    print(f"Наиболее точный результат после уточнения: {best_refined}.")


if __name__ == "__main__":
    main()
