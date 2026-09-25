"""L5 - System mass closure.

m_penalty = m_insulation + m_cooler + m_power + m_radiator + m_lost

Linearity result (derived in the Fall report, s.4): with constant specific
masses, every active-system term is proportional to Q_lift, and m_lost is
proportional to the uncovered leak (Q_leak - Q_lift). m_penalty is therefore
linear in Q_lift, so its minimum over 0 <= Q_lift <= Q_leak sits at an end
point: fully passive or full zero-boil-off. Comparing per-watt costs,

    mu_passive = t_mission / h_fg                                   [kg/W]
    mu_active  = s_cool + s_pow*phi + s_rad*(1 + phi)/q_rad         [kg/W]

and because both multiply the same Q_leak(n), the passive design wins at
every insulation thickness iff mu_passive < mu_active, i.e. for

    t_mission < t* = h_fg * mu_active .

For a vented tank t* is independent of insulation, tank size and environment.
A closed tank that can store energy E_store before first vent shifts it to
t* + E_store/Q_leak = t* + t_hold (derivation in report).
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class MassBreakdown:
    m_insulation: float
    m_cooler: float
    m_power: float
    m_radiator: float
    m_lost: float

    @property
    def m_penalty(self) -> float:
        return (self.m_insulation + self.m_cooler + self.m_power
                + self.m_radiator + self.m_lost)


def mass_breakdown(m_insulation: float, Q_lift: float, s_cool: float, W_in: float,
                   s_pow: float, A_rad: float, s_rad: float, mdot_boil: float,
                   t_mission: float) -> MassBreakdown:
    return MassBreakdown(
        m_insulation=m_insulation,
        m_cooler=Q_lift * s_cool,
        m_power=W_in * s_pow,
        m_radiator=A_rad * s_rad,
        m_lost=mdot_boil * t_mission,
    )


def mu_passive(t_mission: float, h_fg: float) -> float:
    """Mass cost of one watt of leak left to boil for t_mission [kg/W]."""
    return t_mission / h_fg


def mu_active(s_cool: float, phi: float, s_pow: float, s_rad: float, q_rad: float) -> float:
    """Mass cost of removing one watt of leak with a cryocooler [kg/W].
    phi = W_in/Q_lift, q_rad = radiator net flux [W/m^2]."""
    return s_cool + s_pow * phi + s_rad * (1.0 + phi) / q_rad


def t_crossover_vented(h_fg: float, mu_act: float) -> float:
    """Passive/active crossover duration for a vented tank [s]."""
    return h_fg * mu_act
