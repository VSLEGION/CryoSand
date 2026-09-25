"""Core data structures (PROJECT.md s.6). All frozen; all SI."""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum, auto


@dataclass(frozen=True)
class Propellant:
    name: str  # CoolProp name: 'Hydrogen', 'Oxygen', 'Methane'


@dataclass(frozen=True)
class SaturatedState:
    """Saturated fluid properties at the storage pressure. Filled by
    model.saturated_state() from CoolProp; passed explicitly to physics."""
    fluid: str
    P: float        # Pa
    T_sat: float    # K
    h_fg: float     # J/kg
    rho_l: float    # kg/m^3
    rho_v: float    # kg/m^3
    source: str = "CoolProp"


@dataclass(frozen=True)
class TankGeometry:
    """Cylinder with hemispherical end caps."""
    radius: float           # m
    length: float           # m, cylindrical section only
    wall_thickness: float   # m
    wall_material: str
    fill_fraction: float    # liquid volume / total volume

    @property
    def volume(self) -> float:
        r, L = self.radius, self.length
        return math.pi * r**2 * L + 4.0 / 3.0 * math.pi * r**3

    @property
    def surface_area(self) -> float:
        r, L = self.radius, self.length
        return 2.0 * math.pi * r * L + 4.0 * math.pi * r**2

    @property
    def wetted_area(self) -> float:
        """First-order estimate: fill_fraction * surface_area. Orientation and
        microgravity liquid position are not modelled at this fidelity."""
        return self.fill_fraction * self.surface_area


@dataclass(frozen=True)
class Blanket:
    thickness: float       # m
    k_eff: float           # W/(m K); residual gas conduction absorbed here
    areal_density: float   # kg/m^2


@dataclass(frozen=True)
class InsulationStack:
    blankets: tuple[Blanket, ...]
    strut_conductance: float  # W/K, parallel path

    @property
    def n_blankets(self) -> int:
        return len(self.blankets)


@dataclass(frozen=True)
class Environment:
    q_absorbed: float     # W/m^2, area-averaged absorbed flux on outer surface
    T_sink: float         # K
    eps_surface: float    # outer-surface IR emissivity (needed for T_surface)


@dataclass(frozen=True)
class Cryocooler:
    Q_lift: float         # W; 0 means fully passive
    T_cold: float         # K
    T_reject: float       # K
    eta_carnot: float     # fraction of Carnot COP achieved
    specific_mass: float  # kg per W of lift


@dataclass(frozen=True)
class SystemMassParams:
    """Specific masses for L5. PROVISIONAL until the Stage 1 register is done."""
    specific_mass_power: float     # kg per W of electrical input
    areal_density_radiator: float  # kg/m^2
    eps_radiator: float            # radiator IR emissivity
    T_radiator: float              # K (held fixed; an interior optimum exists)


@dataclass(frozen=True)
class Mission:
    name: str
    duration: float        # s
    environment: Environment


class UllageClosure(Enum):
    VENTED = auto()
    HOMOGENEOUS = auto()   # slow bound
    SURFACE = auto()       # fast bound
