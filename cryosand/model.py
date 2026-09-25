"""Assembles L0-L5. Steady driver plus closed-tank pressure histories.

This is the only layer (with trade.py) that talks to CoolProp. Every entry
point that depends on ullage physics takes an explicit UllageClosure.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from cryosand.core import properties as props
from cryosand.core.types import (Cryocooler, Environment, InsulationStack, SaturatedState,
                                 SystemMassParams, TankGeometry, UllageClosure)
from cryosand.layers import l0_environment as L0
from cryosand.layers import l1_insulation as L1
from cryosand.layers import l2_wall as L2
from cryosand.layers import l3_ullage as L3
from cryosand.layers import l4_cooling as L4
from cryosand.layers import l5_closure as L5


# ------------------------------------------------------------ fluid state ---
def saturated_state(fluid: str, P: float) -> SaturatedState:
    """Saturated properties at storage pressure P [Pa], from CoolProp."""
    return SaturatedState(fluid=fluid, P=P, T_sat=props.T_sat(fluid, P),
                          h_fg=props.h_fg(fluid, P), rho_l=props.rho_l_sat(fluid, P),
                          rho_v=props.rho_v_sat(fluid, P), source="CoolProp")


# ------------------------------------------------------ closed-tank bounds ---
@dataclass(frozen=True)
class ClosedTankSpec:
    P_max: float                  # Pa, vent set point
    surf_fraction: float = 1.0    # SURFACE only: m_surf / m_liquid (swept)


class LiquidFullError(ValueError):
    """The (warm zone of the) tank fills with liquid before reaching P_max."""


def homogeneous_initial(fluid: str, geom: TankGeometry, P0: float):
    """(m_total, rho_mix, u0) for a saturated tank at P0 and the given fill."""
    return L3.two_phase_mixture(geom.fill_fraction, geom.volume,
                                props.rho_l_sat(fluid, P0), props.rho_v_sat(fluid, P0),
                                props.u_l_sat(fluid, P0), props.u_v_sat(fluid, P0))


def warm_zone_initial(closure: UllageClosure, fluid: str, geom: TankGeometry, P0: float,
                      spec: ClosedTankSpec):
    """(m, rho, u0) of the zone that absorbs Q_net: the whole tank for
    HOMOGENEOUS, the interface layer plus ullage vapour for SURFACE."""
    if closure is UllageClosure.HOMOGENEOUS:
        return homogeneous_initial(fluid, geom, P0)
    if closure is UllageClosure.SURFACE:
        return L3.two_zone_warm_state(geom.fill_fraction, geom.volume, spec.surf_fraction,
                                      props.rho_l_sat(fluid, P0), props.rho_v_sat(fluid, P0),
                                      props.u_l_sat(fluid, P0), props.u_v_sat(fluid, P0))
    raise ValueError("defined for the closed closures only")


def stored_energy(closure: UllageClosure, fluid: str, geom: TankGeometry, P0: float,
                  spec: ClosedTankSpec | None) -> float:
    """Energy a closed tank absorbs between P0 and spec.P_max [J].

    Raises LiquidFullError if the absorbing zone becomes liquid-full first:
    then the vent pressure is reached hydraulically and neither bound applies.
    """
    if closure is UllageClosure.VENTED:
        return 0.0
    if spec is None:
        raise ValueError("closed closures need a ClosedTankSpec")
    m, rho, u0 = warm_zone_initial(closure, fluid, geom, P0, spec)
    if rho > props.rho_l_sat(fluid, spec.P_max):
        raise LiquidFullError(
            f"{fluid}, fill {geom.fill_fraction}: liquid-full before "
            f"{spec.P_max:.0f} Pa (rho {rho:.2f} > rho_l_sat {props.rho_l_sat(fluid, spec.P_max):.2f})")
    u_max = props.u_from_rho_P(fluid, rho, spec.P_max)
    return L3.stored_energy_homogeneous(m, u0, u_max)


def pressure_history(closure: UllageClosure, fluid: str, geom: TankGeometry, P0: float,
                     Q_net: float, t: np.ndarray, spec: ClosedTankSpec) -> np.ndarray:
    """Closed-tank pressure P(t) [Pa] under either bound (no venting). Exact:
    u of the absorbing zone is linear in t, and P follows from the EOS."""
    m, rho, u0 = warm_zone_initial(closure, fluid, geom, P0, spec)
    du = L3.du_dt_homogeneous(Q_net, m)
    return np.array([props.P_from_rho_u(fluid, rho, u0 + du * ti) for ti in t])


# ------------------------------------------------------------ steady run ---
@dataclass(frozen=True)
class SteadyResult:
    T_surface: float
    Q_leak: float
    Q_strut: float
    R_wall: float
    R_blanket: float
    Q_net: float
    W_in: float
    Q_reject: float
    A_rad: float
    mdot_boil: float
    E_store: float
    t_hold: float
    mass: L5.MassBreakdown


def run_steady(sat: SaturatedState, geom: TankGeometry, stack: InsulationStack,
               env: Environment, cooler: Cryocooler, sysp: SystemMassParams,
               k_wall: float, t_mission: float, closure: UllageClosure,
               E_store: float = 0.0) -> SteadyResult:
    """End-to-end steady evaluation. For closed closures pass E_store from
    stored_energy(); boil-off begins after t_hold = E_store / Q_net."""
    A = geom.surface_area
    T_s = L0.surface_temperature(env.q_absorbed, env.eps_surface, env.T_sink)
    R_w = L2.R_wall_plane(geom.wall_thickness, k_wall, A)
    leak = L1.heat_leak(stack, A, T_s, sat.T_sat, R_w)
    Q_net = leak.Q_total - cooler.Q_lift
    W = L4.work_input(cooler.Q_lift, cooler.T_cold, cooler.T_reject, cooler.eta_carnot) \
        if cooler.Q_lift > 0 else 0.0
    Q_rej = L4.heat_rejected(cooler.Q_lift, W)
    A_rad = L4.radiator_area(Q_rej, sysp.eps_radiator, sysp.T_radiator, env.T_sink)
    mdot = L3.boiloff_rate_vented(Q_net, sat.h_fg)
    if closure is UllageClosure.VENTED:
        E_store = 0.0
    t_hold = L3.hold_time(E_store, Q_net)
    t_venting = max(t_mission - t_hold, 0.0)
    mass = L5.mass_breakdown(L1.insulation_mass(stack, A), cooler.Q_lift, cooler.specific_mass,
                             W, sysp.specific_mass_power, A_rad, sysp.areal_density_radiator,
                             mdot, t_venting)
    return SteadyResult(T_s, leak.Q_total, leak.Q_strut, R_w, leak.R_blanket, Q_net, W, Q_rej,
                        A_rad, mdot, E_store, t_hold, mass)
