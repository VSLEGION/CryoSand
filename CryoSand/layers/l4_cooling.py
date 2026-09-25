"""L4 - Active cooling and heat rejection.

COP_carnot = T_cold / (T_hot - T_cold)
W_in       = Q_lift / (eta_carnot * COP_carnot)
Q_reject   = Q_lift + W_in          <- the radiator rejects BOTH
A_rad      = Q_reject / (eps * sigma * (T_rad^4 - T_sink^4))

T_rad has an interior optimum (larger T_rad shrinks A_rad but raises W_in via
T_hot). Held fixed here; the optimum is noted, not exploited.
Thermodynamic vent system: deferred.
"""
from __future__ import annotations

from cryosand.core.constants import SIGMA


def cop_carnot(T_cold: float, T_hot: float) -> float:
    if not T_hot > T_cold > 0:
        raise ValueError("need T_hot > T_cold > 0")
    return T_cold / (T_hot - T_cold)


def specific_power(T_cold: float, T_hot: float, eta_carnot: float) -> float:
    """phi = W_in / Q_lift [W/W] = 1 / (eta * COP_carnot)."""
    return 1.0 / (eta_carnot * cop_carnot(T_cold, T_hot))


def work_input(Q_lift: float, T_cold: float, T_hot: float, eta_carnot: float) -> float:
    """Electrical input power W_in [W]."""
    return Q_lift * specific_power(T_cold, T_hot, eta_carnot)


def heat_rejected(Q_lift: float, W_in: float) -> float:
    """Q_reject = Q_lift + W_in [W] (first law on the cooler)."""
    return Q_lift + W_in


def radiator_flux(eps: float, T_rad: float, T_sink: float) -> float:
    """Net radiated flux per unit radiator area [W/m^2] (one-sided)."""
    return eps * SIGMA * (T_rad**4 - T_sink**4)


def radiator_area(Q_reject: float, eps: float, T_rad: float, T_sink: float) -> float:
    """A_rad [m^2]."""
    return Q_reject / radiator_flux(eps, T_rad, T_sink)


def cooler_mass_strobridge(Q_lift: float, T_cold: float, T_reject: float) -> float:
    """Cryocooler mass [kg] from the Strobridge correlation scaled by 0.2 for
    modern coolers: m = 0.2 * Q^0.7 * ((T_h - T_c)/T_c)^1.45, Q in W.
    Source: Plachta & Kittel (2003), NTRS 20030067928.
    Concave in Q (economy of scale), so specific mass falls as Q^-0.3."""
    if Q_lift <= 0.0:
        return 0.0
    return 0.2 * Q_lift**0.7 * ((T_reject - T_cold) / T_cold) ** 1.45
