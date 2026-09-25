"""Parameter sweeps and the passive/active crossover.

Design variables: number of MLI blankets n (integer) and cryocooler lift.
Because m_penalty is concave in Q_lift (linear boil-off/power/radiator terms
plus a concave Q^0.7 cooler mass), its minimum over 0 <= Q_lift <= Q_leak is
at an end point: fully passive or full zero-boil-off (ZBO). So the trade
compares min_n(passive) against min_n(ZBO).

Every function that depends on ullage physics takes an explicit closure and
the stored energy E_store (0 for VENTED) computed by model.stored_energy.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from cryosand.core.types import UllageClosure
from cryosand.layers import l1_insulation as L1
from cryosand.layers import l4_cooling as L4


@dataclass(frozen=True)
class TradeCase:
    h_fg: float            # J/kg
    T_cold: float          # K (propellant T_sat)
    T_hot: float           # K (outer surface)
    A: float               # m^2
    G_strut: float         # W/K
    t_blanket: float       # m
    k_eff: float           # W/(m K)
    rho_A_blanket: float   # kg/m^2 per blanket
    R_wall: float          # K/W
    eta_carnot: float
    T_reject: float        # K
    s_pow: float           # kg/W
    s_rad: float           # kg/m^2
    q_rad: float           # W/m^2
    n_max: int = 180


def Q_leak(c: TradeCase, n: int) -> float:
    R_b = n * c.t_blanket / (c.k_eff * c.A)
    return (c.T_hot - c.T_cold) / (R_b + c.R_wall) + c.G_strut * (c.T_hot - c.T_cold)


def passive_mass(c: TradeCase, n: int, t: float, closure: UllageClosure, E_store: float) -> float:
    Q = Q_leak(c, n)
    E_s = 0.0 if closure is UllageClosure.VENTED else E_store
    return n * c.rho_A_blanket * c.A + max(Q * t - E_s, 0.0) / c.h_fg


def zbo_mass(c: TradeCase, n: int) -> float:
    Q = Q_leak(c, n)
    phi = L4.specific_power(c.T_cold, c.T_reject, c.eta_carnot)
    return (n * c.rho_A_blanket * c.A + L4.cooler_mass_strobridge(Q, c.T_cold, c.T_reject)
            + c.s_pow * phi * Q + c.s_rad * (1.0 + phi) * Q / c.q_rad)


def best_passive(c, t, closure, E_store):
    return min((passive_mass(c, n, t, closure, E_store), n) for n in range(1, c.n_max + 1))


def best_zbo(c):
    return min((zbo_mass(c, n), n) for n in range(1, c.n_max + 1))


def t_crossover(c: TradeCase, closure: UllageClosure, E_store: float,
                t_lo: float = 3600.0, t_hi: float = 3.2e9) -> float:
    """Mission duration [s] at which ZBO and passive optimum masses are equal."""
    m_zbo = best_zbo(c)[0]
    f = lambda t: best_passive(c, t, closure, E_store)[0] - m_zbo  # noqa: E731
    if f(t_lo) > 0:
        return t_lo
    if f(t_hi) < 0:
        return math.inf
    for _ in range(200):
        tm = math.sqrt(t_lo * t_hi)
        if f(tm) > 0:
            t_hi = tm
        else:
            t_lo = tm
        if t_hi / t_lo - 1.0 < 1e-9:
            break
    return math.sqrt(t_lo * t_hi)
