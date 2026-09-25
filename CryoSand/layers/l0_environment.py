"""L0 - Orbital environment.

q_abs = alpha * (F_sun*G_sun + F_E*a*G_sun) + eps * F_E * G_IR     [W/m^2]

Area-averaged over a sphere. F_sun = 1/4 (projected area pi r^2 over total
area 4 pi r^2). F_E is the sphere-to-Earth view factor. Using F_E for albedo
is an upper bound (it assumes the subsolar point); orbit-averaging with albedo
phase is a spring upgrade together with attitude dependence.

Outer-surface radiative equilibrium (inward heat leak << absorbed flux):
eps * sigma * (T_s^4 - T_sink^4) = q_abs
"""
from __future__ import annotations

import math

from cryosand.core.constants import SIGMA

F_SUN_SPHERE = 0.25


def view_factor_sphere_to_earth(altitude: float, R_planet: float) -> float:
    """View factor from a small sphere to a planet [-].
    F = (1 - sqrt(1 - (R/(R+h))^2)) / 2   (Gilmore 2002, App. / standard result).
    h -> 0 gives 1/2 (half the sky is planet); h -> inf gives 0."""
    x = R_planet / (R_planet + altitude)
    return 0.5 * (1.0 - math.sqrt(1.0 - x * x))


def solar_flux(distance_AU: float, G_sun_1AU: float) -> float:
    """Inverse-square solar irradiance [W/m^2]."""
    return G_sun_1AU / distance_AU**2


def absorbed_flux(alpha: float, eps: float, G_sun: float, albedo: float,
                  G_IR: float, F_planet: float, F_sun: float = F_SUN_SPHERE) -> float:
    """Area-averaged absorbed flux [W/m^2]."""
    return alpha * (F_sun * G_sun + F_planet * albedo * G_sun) + eps * F_planet * G_IR


def surface_temperature(q_abs: float, eps: float, T_sink: float) -> float:
    """Radiative-equilibrium outer-surface temperature [K]."""
    return (q_abs / (eps * SIGMA) + T_sink**4) ** 0.25
