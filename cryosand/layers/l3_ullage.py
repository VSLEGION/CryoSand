"""L3 - Ullage and liquid. Three closures, all first-class (PROJECT.md s.2).

Pure functions only: property values arrive as arguments. The CoolProp-backed
pressure histories are assembled in cryosand.model.

VENTED (constant pressure)
    mdot_boil = Q_net / h_fg

CLOSED, HOMOGENEOUS (slow bound) - whole tank at one state, V and m fixed
    du/dt = Q_net / m_total ,  rho = m_total / V ,  P = P(rho, u)

CLOSED, SURFACE (fast bound) - all energy into an interface layer m_surf
    dT_int/dt = Q_net / (m_surf * cp) ,  P = P_sat(T_int)

Energy stored before the first vent at P_max, and the hold time, follow
directly because Q_net is constant in the steady-leak problem:
    t_hold = E_store / Q_net
"""
from __future__ import annotations

import math


# ---------------------------------------------------------------- vented ---
def boiloff_rate_vented(Q_net: float, h_fg: float) -> float:
    """Vented boil-off [kg/s]. Q_net <= 0 (net cooling) gives zero: liquid
    subcooling under net refrigeration is not modelled at this fidelity."""
    if h_fg <= 0:
        raise ValueError("h_fg must be positive")
    return Q_net / h_fg if Q_net > 0.0 else 0.0


# ----------------------------------------------------------- homogeneous ---
def two_phase_mixture(fill_fraction: float, V: float, rho_l: float, rho_v: float,
                      u_l: float, u_v: float) -> tuple[float, float, float]:
    """Initial saturated two-phase state.
    Returns (m_total [kg], rho_mix [kg/m^3], u_mix [J/kg])."""
    m_l = rho_l * fill_fraction * V
    m_v = rho_v * (1.0 - fill_fraction) * V
    m = m_l + m_v
    return m, m / V, (m_l * u_l + m_v * u_v) / m


def du_dt_homogeneous(Q_net: float, m_total: float) -> float:
    """Specific internal energy rate [J/(kg s)] for the closed, well-mixed tank."""
    return Q_net / m_total


def stored_energy_homogeneous(m_total: float, u_0: float, u_at_Pmax: float) -> float:
    """Energy absorbed before P reaches P_max [J]. u_at_Pmax = u(rho_mix, P_max)."""
    return m_total * (u_at_Pmax - u_0)


# --------------------------------------------------------------- surface ---
def dT_dt_surface(Q_net: float, m_surf: float, cp: float) -> float:
    """Interface-layer heating rate [K/s]."""
    return Q_net / (m_surf * cp)


def stored_energy_surface(m_surf: float, cp_of_T, T_0: float, T_max: float,
                          n: int = 401) -> float:
    """E = m_surf * integral(cp dT) from T_0 to T_max = T_sat(P_max) [J].
    cp_of_T is a callable (CoolProp saturated-liquid cp in production)."""
    h = (T_max - T_0) / (n - 1)
    s = 0.5 * (cp_of_T(T_0) + cp_of_T(T_max))
    s += sum(cp_of_T(T_0 + i * h) for i in range(1, n - 1))
    return m_surf * s * h


def m_surf_conduction(rho_l: float, alpha_l: float, A_interface: float, t: float) -> float:
    """Physically anchored interface-layer mass [kg] for pure conduction.

    Constant flux q into a semi-infinite liquid (Carslaw & Jaeger 1959, s.2.9):
    surface rise dT_s = (2q/k) sqrt(alpha t / pi). Defining m_eff by
    E = m_eff cp dT_s with E = q A t gives
        m_eff = rho * A * sqrt(pi alpha t) / 2 .
    A stagnant, stably stratified layer; a reference point for the m_surf sweep,
    not a replacement for it.
    """
    return 0.5 * rho_l * A_interface * math.sqrt(math.pi * alpha_l * t)


# ------------------------------------------------------------------ both ---
def hold_time(E_store: float, Q_net: float) -> float:
    """Time to first vent [s] at constant Q_net > 0."""
    if Q_net <= 0.0:
        return math.inf
    return E_store / Q_net
