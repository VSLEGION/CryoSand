"""L1 - Insulation stack.

Blankets are series resistances from the outer surface to the tank wall.
Struts and penetrations are a PARALLEL conductance: they bypass the blanket.

    G_total = 1 / (sum R_blanket + R_wall) + G_strut
    Q_leak  = G_total * (T_hot - T_cold)

k_eff is an effective blanket conductivity. Radiation between shields, solid
conduction through spacers, and residual gas conduction are all absorbed into
k_eff at this fidelity. k_eff is held constant; in reality it depends strongly
on the boundary temperatures (known trap, flagged in the report). The modified
Lockheed correlation is a spring upgrade.
"""
from __future__ import annotations

from dataclasses import dataclass

from cryosand.core.network import R_series, R_slab, heat_rate
from cryosand.core.types import InsulationStack


@dataclass(frozen=True)
class HeatLeak:
    Q_blanket: float  # W, through blankets + wall
    Q_strut: float    # W, parallel path
    R_blanket: float  # K/W, all blankets in series
    R_wall: float     # K/W

    @property
    def Q_total(self) -> float:
        return self.Q_blanket + self.Q_strut


def R_blankets(stack: InsulationStack, A: float) -> float:
    """Series resistance of all blankets [K/W]."""
    return R_series(R_slab(b.thickness, b.k_eff, A) for b in stack.blankets)


def heat_leak(stack: InsulationStack, A: float, T_hot: float, T_cold: float,
              R_wall: float) -> HeatLeak:
    """Steady heat leak into the tank [W], split by path."""
    R_b = R_blankets(stack, A)
    Q_b = heat_rate(T_hot, T_cold, R_b + R_wall)
    Q_s = stack.strut_conductance * (T_hot - T_cold)
    return HeatLeak(Q_blanket=Q_b, Q_strut=Q_s, R_blanket=R_b, R_wall=R_wall)


def insulation_mass(stack: InsulationStack, A: float) -> float:
    """m_insulation = sum(areal_density_i) * A [kg]."""
    return sum(b.areal_density for b in stack.blankets) * A
