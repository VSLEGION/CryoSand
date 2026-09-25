"""Thermal resistance network primitives. R in K/W, G = 1/R in W/K.

Analogy: temperature difference <-> voltage, heat rate <-> current.
"""
from __future__ import annotations

import math
from typing import Iterable


def R_slab(thickness: float, k: float, area: float) -> float:
    """Plane-wall conduction resistance R = t / (k A) [K/W]."""
    if thickness < 0 or k <= 0 or area <= 0:
        raise ValueError("need thickness >= 0, k > 0, area > 0")
    return thickness / (k * area)


def R_series(resistances: Iterable[float]) -> float:
    """Series resistances add: R = sum(R_i)."""
    return float(sum(resistances))


def R_parallel(resistances: Iterable[float]) -> float:
    """Parallel resistances: 1/R = sum(1/R_i). An infinite R is an open path."""
    G = sum(0.0 if math.isinf(R) else 1.0 / R for R in resistances)
    return math.inf if G == 0.0 else 1.0 / G


def heat_rate(T_hot: float, T_cold: float, R: float) -> float:
    """Q = (T_hot - T_cold) / R [W]."""
    return (T_hot - T_cold) / R


def k_face(k_P: float, k_N: float) -> float:
    """Harmonic-mean face conductivity between two equal half-cells.

    Two half-cell resistances in series: (h/2)/k_P + (h/2)/k_N = h / k_face
    => k_face = 2 k_P k_N / (k_P + k_N). The arithmetic mean overstates flux
    when k_P and k_N differ by orders of magnitude (e.g. MLI next to aluminium).
    """
    return 2.0 * k_P * k_N / (k_P + k_N)


def k_bar(k_of_T, T_cold: float, T_hot: float, n: int = 2001) -> float:
    """Temperature-averaged conductivity k_bar = (1/dT) * integral(k dT) [W/(m K)].

    k_of_T is a callable k(T). Trapezoidal rule on n points; exact for linear k.
    """
    if T_hot <= T_cold:
        raise ValueError("T_hot must exceed T_cold")
    h = (T_hot - T_cold) / (n - 1)
    s = 0.5 * (k_of_T(T_cold) + k_of_T(T_hot))
    s += sum(k_of_T(T_cold + i * h) for i in range(1, n - 1))
    return s * h / (T_hot - T_cold)
