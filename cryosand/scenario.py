"""Scenario builder: turns data/parameters.yaml into model inputs.

The parameter register is the ONLY place parameter values live. Both
scripts/run_case.py and scripts/make_figures.py build their cases here, so a
change to the register reaches every result.

Not a physics module: it reads configuration and calls CoolProp via
cryosand.core.properties, like model.py and trade.py.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml

from cryosand.core import constants as C
from cryosand.core import properties as props
from cryosand.core.types import TankGeometry
from cryosand.layers import l0_environment as L0
from cryosand.layers import l2_wall as L2
from cryosand.layers import l4_cooling as L4
from cryosand.trade import TradeCase

REGISTER = Path(__file__).resolve().parents[1] / "data" / "parameters.yaml"
FLUIDS = ("Hydrogen", "Methane", "Oxygen")


def _v(entry) -> float:
    """Register entries are either bare numbers or {value: ..., ...}. Coerced
    to float: PyYAML reads e.g. 3.0e5 as a string, and a non-numeric typo
    should fail here rather than deep inside the physics."""
    return float(entry["value"] if isinstance(entry, dict) else entry)


@dataclass(frozen=True)
class Params:
    """Flat, typed view of the register. SI units."""
    eta_carnot: dict          # fluid -> fraction of Carnot
    T_reject: float           # K
    s_pow: float              # kg/W
    s_rad: float              # kg/m^2
    eps_rad: float
    T_rad: float              # K
    mli_layer_areal_density: float   # kg/m^2 per layer
    mli_blanket_layers: int
    mli_blanket_thickness: float     # m
    k_eff: float              # W/(m K)
    alpha_out: float
    eps_out: float
    G_strut: float            # W/K
    wall_thickness: float     # m
    k_wall: float             # W/(m K)
    P_storage: float          # Pa
    fill: float
    P_vent: float             # Pa
    radius: float             # m
    length_over_radius: float
    unit_layers: int
    environments: dict        # key -> {label, distance_AU, altitude, planet}
    source: str               # path of the register used

    @property
    def unit_thickness(self) -> float:
        """Thickness of one design-unit sub-blanket [m]."""
        return self.mli_blanket_thickness * self.unit_layers / self.mli_blanket_layers

    @property
    def unit_areal_density(self) -> float:
        """Areal density of one design-unit sub-blanket [kg/m^2]."""
        return self.unit_layers * self.mli_layer_areal_density


def load_params(path: str | Path | None = None) -> Params:
    path = Path(path) if path else REGISTER
    d = yaml.safe_load(path.read_text())
    cc, mli, b = d["cryocooler"], d["mli"], d["baseline"]
    return Params(
        eta_carnot={f: _v(v) for f, v in cc["eta_carnot"].items()},
        T_reject=_v(cc["T_reject"]),
        s_pow=_v(d["power_system"]["specific_mass"]),
        s_rad=_v(d["radiator"]["areal_density"]),
        eps_rad=_v(d["radiator"]["emissivity"]),
        T_rad=_v(d["radiator"]["T_radiator"]),
        mli_layer_areal_density=_v(mli["per_layer_areal_density"]),
        mli_blanket_layers=int(mli["blanket"]["layers"]),
        mli_blanket_thickness=_v(mli["blanket"]["thickness"]),
        k_eff=_v(mli["blanket"]["k_eff"]),
        alpha_out=_v(mli["outer_cover"]["alpha"]),
        eps_out=_v(mli["outer_cover"]["eps"]),
        G_strut=_v(d["struts"]["conductance"]),
        wall_thickness=_v(d["wall"]["thickness"]),
        k_wall=_v(d["wall"]["k"]),
        P_storage=_v(b["P_storage"]),
        fill=_v(b["fill_fraction"]),
        P_vent=_v(b["P_vent"]),
        radius=_v(b["radius"]),
        length_over_radius=_v(b["length_over_radius"]),
        unit_layers=int(_v(b["design_unit_layers"])),
        environments=d["environments"],
        source=str(path),
    )


# ------------------------------------------------------------ environment ---
def environment(p: Params, key: str) -> dict:
    """Absorbed flux [W/m^2] and outer-surface temperature [K] for a preset."""
    if key not in p.environments:
        raise KeyError(f"unknown environment {key!r}; choose from {sorted(p.environments)}")
    e = p.environments[key]
    G_sun = L0.solar_flux(float(e["distance_AU"]), C.G_SUN_1AU)
    if e.get("planet") == "earth":
        F = L0.view_factor_sphere_to_earth(float(e["altitude"]), C.R_EARTH)
        q = L0.absorbed_flux(p.alpha_out, p.eps_out, G_sun, C.ALBEDO_EARTH, C.G_IR_EARTH, F)
    elif e.get("planet") is None:
        F = 0.0
        q = L0.absorbed_flux(p.alpha_out, p.eps_out, G_sun, 0.0, 0.0, 0.0)
    else:
        raise ValueError(f"planet {e['planet']!r} not modelled")
    return dict(label=e["label"], q_abs=q, F_planet=F,
                T_s=L0.surface_temperature(q, p.eps_out, C.T_DEEP_SPACE))


# ------------------------------------------------------------------ tank ---
def tank(p: Params, radius: float | None = None, fill: float | None = None) -> TankGeometry:
    r = p.radius if radius is None else radius
    return TankGeometry(radius=r, length=p.length_over_radius * r,
                        wall_thickness=p.wall_thickness, wall_material="Al6061",
                        fill_fraction=p.fill if fill is None else fill)


def saturated(fluid: str, p: Params) -> dict:
    """Saturated properties at storage pressure, from CoolProp."""
    P = p.P_storage
    return dict(T_sat=props.T_sat(fluid, P), rho_l=props.rho_l_sat(fluid, P),
                h_fg=props.h_fg(fluid, P))


def trade_case(p: Params, fluid: str, T_s: float, radius: float | None = None) -> TradeCase:
    g = tank(p, radius)
    A = g.surface_area
    s = saturated(fluid, p)
    return TradeCase(h_fg=s["h_fg"], T_cold=s["T_sat"], T_hot=T_s, A=A,
                     G_strut=p.G_strut,  # held constant across scale (provisional)
                     t_blanket=p.unit_thickness, k_eff=p.k_eff,
                     rho_A_blanket=p.unit_areal_density,
                     R_wall=L2.R_wall_plane(p.wall_thickness, p.k_wall, A),
                     eta_carnot=p.eta_carnot[fluid], T_reject=p.T_reject,
                     s_pow=p.s_pow, s_rad=p.s_rad,
                     q_rad=L4.radiator_flux(p.eps_rad, p.T_rad, C.T_DEEP_SPACE))

