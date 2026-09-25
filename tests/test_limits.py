"""Analytic limits (PROJECT.md s.7). Property-free: runs without CoolProp."""
import math

from cryosand.core.network import R_parallel, R_series, R_slab, heat_rate, k_bar, k_face
from cryosand.core.types import (Blanket, Cryocooler, Environment, InsulationStack,
                                 SaturatedState, SystemMassParams, TankGeometry, UllageClosure)
from cryosand.layers import l0_environment as L0
from cryosand.layers import l1_insulation as L1
from cryosand.layers import l2_wall as L2
from cryosand.layers import l3_ullage as L3
from cryosand.layers import l4_cooling as L4
from cryosand.layers import l5_closure as L5
from cryosand.model import run_steady

# Test fixture only: arbitrary but plausible saturated state, NOT a property source.
SAT = SaturatedState("TestFluid", 101325.0, 20.0, 4.0e5, 70.0, 1.3, source="test fixture")
GEOM = TankGeometry(2.0, 6.0, 0.002, "Al6061", 0.9)
SYSP = SystemMassParams(0.025, 5.0, 0.9, 300.0)
ENV = Environment(100.0, 4.0, 0.8)


def _stack(n, G_strut=5e-3):
    return InsulationStack(tuple(Blanket(0.025, 3e-5, 0.6) for _ in range(n)), G_strut)


def _passive():
    return Cryocooler(0.0, 20.0, 300.0, 0.075, 5.0)


# ------------------------------------------------------------- network ---
def test_harmonic_mean():
    assert k_face(2.0, 6.0) == 3.0


def test_harmonic_mean_dominated_by_poor_conductor():
    # MLI next to aluminium: face conductance ~ 2x the poor side, not ~half the good side
    assert k_face(3e-5, 150.0) < 6.001e-5


def test_series_chain():
    A, dT = 7.0, 180.0
    layers = [(0.01, 2.0), (0.002, 150.0), (0.03, 3e-5)]
    R = R_series(R_slab(t, k, A) for t, k in layers)
    R_hand = sum(t / (k * A) for t, k in layers)
    assert abs(heat_rate(dT, 0.0, R) - dT / R_hand) <= 1e-10 * dT / R_hand


def test_parallel():
    assert math.isclose(R_parallel([4.0, 4.0]), 2.0)
    assert R_parallel([math.inf, math.inf]) == math.inf
    assert math.isclose(R_parallel([3.0, math.inf]), 3.0)


def test_k_bar_linear_exact():
    a, b = 0.1, 0.002
    kb = k_bar(lambda T: a + b * T, 20.0, 300.0)
    assert math.isclose(kb, a + b * 160.0, rel_tol=1e-12)


# ------------------------------------------------------------------ L1 ---
def test_blankets_in_series_scale_as_one_over_n():
    A, Th, Tc = 100.0, 250.0, 20.0
    Q1 = L1.heat_leak(_stack(1, 0.0), A, Th, Tc, 0.0).Q_blanket
    Q4 = L1.heat_leak(_stack(4, 0.0), A, Th, Tc, 0.0).Q_blanket
    assert math.isclose(Q1 / Q4, 4.0, rel_tol=1e-12)


def test_struts_are_parallel():
    A, Th, Tc, G = 100.0, 250.0, 20.0, 5e-3
    no = L1.heat_leak(_stack(2, 0.0), A, Th, Tc, 0.0)
    yes = L1.heat_leak(_stack(2, G), A, Th, Tc, 0.0)
    assert yes.Q_blanket == no.Q_blanket            # blanket path untouched
    assert math.isclose(yes.Q_total - no.Q_total, G * (Th - Tc), rel_tol=1e-12)


# ------------------------------------------------------------------ L0 ---
def test_view_factor_limits():
    assert math.isclose(L0.view_factor_sphere_to_earth(0.0, 6.371e6), 0.5)
    assert L0.view_factor_sphere_to_earth(1e12, 6.371e6) < 1e-10


def test_surface_temperature_zero_flux_is_sink():
    assert math.isclose(L0.surface_temperature(0.0, 0.8, 4.0), 4.0)


# ------------------------------------------------------------------ L2 ---
def test_wall_cylinder_tends_to_plane():
    r, t, L, k = 2.0, 1e-5, 6.0, 150.0
    plane = L2.R_wall_plane(t, k, 2 * math.pi * r * L)
    assert math.isclose(L2.R_wall_cylinder(r, t, L, k), plane, rel_tol=1e-5)


def test_wall_is_negligible_but_computed():
    A = GEOM.surface_area
    R_w = L2.R_wall_plane(GEOM.wall_thickness, 10.0, A)  # k=10 is pessimistic for Al at 20 K
    R_b = L1.R_blankets(_stack(1), A)
    assert R_w / R_b < 1e-6


# ------------------------------------------------------------------ L3 ---
def test_zero_leak():
    assert L3.boiloff_rate_vented(0.0, 4.46e5) == 0.0


def test_m_surf_conduction_consistent_with_semi_infinite_solution():
    rho, cp, k, A, q, t = 70.0, 9700.0, 0.1, 10.0, 0.3, 86400.0
    alpha = k / (rho * cp)
    dTs = 2 * q / k * math.sqrt(alpha * t / math.pi)       # Carslaw & Jaeger
    m = L3.m_surf_conduction(rho, alpha, A, t)
    assert math.isclose(m * cp * dTs, q * A * t, rel_tol=1e-12)


def test_hold_time():
    assert L3.hold_time(1e6, 10.0) == 1e5
    assert L3.hold_time(1e6, 0.0) == math.inf


# ------------------------------------------------------------------ L4 ---
def test_carnot_limits():
    assert math.isclose(L4.specific_power(20.0, 300.0, 1.0), 14.0)
    assert math.isclose(L4.specific_power(90.0, 300.0, 1.0), 210.0 / 90.0)


def test_radiator_rejects_lift_plus_work():
    W = L4.work_input(10.0, 20.0, 300.0, 0.1)
    assert math.isclose(L4.heat_rejected(10.0, W), 10.0 + 140.0 * 10.0)


# --------------------------------------------------------- end to end ---
def test_zbo_balance():
    base = run_steady(SAT, GEOM, _stack(2), ENV, _passive(), SYSP, 150.0, 1e7,
                      UllageClosure.VENTED)
    cooler = Cryocooler(base.Q_leak, 20.0, 300.0, 0.075, 5.0)
    zbo = run_steady(SAT, GEOM, _stack(2), ENV, cooler, SYSP, 150.0, 1e7, UllageClosure.VENTED)
    assert abs(zbo.mdot_boil) <= 1e-6
    assert zbo.mass.m_lost == 0.0


def _best(t, mode, n_max=12):
    best = math.inf
    for n in range(1, n_max + 1):
        st = _stack(n)
        r0 = run_steady(SAT, GEOM, st, ENV, _passive(), SYSP, 150.0, t, UllageClosure.VENTED)
        if mode == "active":
            c = Cryocooler(r0.Q_leak, 20.0, 300.0, 0.075, 5.0)
            r0 = run_steady(SAT, GEOM, st, ENV, c, SYSP, 150.0, t, UllageClosure.VENTED)
        best = min(best, r0.mass.m_penalty)
    return best


def test_crossover_matches_linearity_result():
    """Brute-force optimisation over insulation reproduces t* = h_fg * mu_active."""
    phi = L4.specific_power(20.0, 300.0, 0.075)
    q_rad = L4.radiator_flux(0.9, 300.0, 4.0)
    t_star = L5.t_crossover_vented(SAT.h_fg, L5.mu_active(5.0, phi, 0.025, 5.0, q_rad))
    assert _best(0.98 * t_star, "passive") < _best(0.98 * t_star, "active")
    assert _best(1.02 * t_star, "passive") > _best(1.02 * t_star, "active")


def test_two_zone_collapses_to_homogeneous():
    """surf_fraction = 1: the warm zone is the whole tank (property-free)."""
    args = (0.9, 70.0, 70.8, 1.3, -3.0e3, 4.5e5)   # fill, V, rho_l, rho_v, u_l, u_v
    V_, fill = args[1], args[0]
    hom = L3.two_phase_mixture(fill, V_, *args[2:])
    two = L3.two_zone_warm_state(fill, V_, 1.0, *args[2:])
    assert all(math.isclose(a, b, rel_tol=1e-12) for a, b in zip(hom, two))


def test_two_zone_warm_zone_shrinks_with_layer():
    small = L3.two_zone_warm_state(0.9, 70.0, 1e-3, 70.8, 1.3, -3.0e3, 4.5e5)
    big = L3.two_zone_warm_state(0.9, 70.0, 0.5, 70.8, 1.3, -3.0e3, 4.5e5)
    assert small[0] < big[0]


def test_trade_leak_matches_layer_physics():
    """trade.Q_leak re-states L1 for speed; it must agree with L1.heat_leak."""
    from cryosand.trade import Q_leak, TradeCase
    A, Th, Tc, n, R_w = GEOM.surface_area, 259.0, 20.3, 7, 2e-6
    c = TradeCase(h_fg=4.46e5, T_cold=Tc, T_hot=Th, A=A, G_strut=5e-3, t_blanket=0.025,
                  k_eff=3e-5, rho_A_blanket=0.6, R_wall=R_w, eta_carnot=0.075,
                  T_reject=300.0, s_pow=0.025, s_rad=5.0, q_rad=400.0)
    ref = L1.heat_leak(_stack(n), A, Th, Tc, R_w).Q_total
    assert math.isclose(Q_leak(c, n), ref, rel_tol=1e-12)


# ------------------------------------------------------ register / trade ---
def test_register_loads_report_baseline():
    """The register reproduces the report's baseline inputs (property-free)."""
    from cryosand import scenario as S
    p = S.load_params()
    assert math.isclose(p.unit_thickness, 0.025 / 3, rel_tol=1e-12)
    assert math.isclose(p.unit_areal_density, 0.2, rel_tol=1e-12)
    assert isinstance(p.P_vent, float) and p.P_vent == 3.0e5
    assert math.isclose(S.environment(p, "leo")["T_s"], 259.0, abs_tol=0.05)
    assert S.tank(p).length == 1.5 * p.radius


def test_mass_breakdowns_sum_to_totals():
    from cryosand.trade import (TradeCase, passive_breakdown, passive_mass, zbo_breakdown,
                                zbo_mass)
    c = TradeCase(h_fg=4.46e5, T_cold=20.3, T_hot=259.0, A=88.0, G_strut=5e-3,
                  t_blanket=0.025 / 3, k_eff=3e-5, rho_A_blanket=0.2, R_wall=2e-6,
                  eta_carnot=0.075, T_reject=300.0, s_pow=0.025, s_rad=5.0, q_rad=410.0)
    for n in (1, 6, 40):
        assert math.isclose(sum(zbo_breakdown(c, n).values()), zbo_mass(c, n), rel_tol=1e-12)
        for cl, E in ((UllageClosure.VENTED, 0.0), (UllageClosure.HOMOGENEOUS, 1e8)):
            assert math.isclose(sum(passive_breakdown(c, n, 3e6, cl, E).values()),
                                passive_mass(c, n, 3e6, cl, E), rel_tol=1e-12)
