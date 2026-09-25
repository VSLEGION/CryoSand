"""Tests that need CoolProp. Skipped automatically where it is not installed.

test_bounds_ordered protects the entire research claim (PROJECT.md s.7).
"""
import math

import numpy as np
import pytest

pytest.importorskip("CoolProp")

from cryosand.core import properties as props  # noqa: E402
from cryosand.core.types import TankGeometry, UllageClosure  # noqa: E402
from cryosand.layers import l3_ullage as L3  # noqa: E402
from cryosand.model import (ClosedTankSpec, LiquidFullError, homogeneous_initial,  # noqa: E402
                            pressure_history, stored_energy, warm_zone_initial)

P0 = 101325.0
P_MAX = 3.0e5
FLUIDS = ["Hydrogen", "Oxygen", "Methane"]
FILLS = [0.1, 0.3, 0.6, 0.9]
FRACTIONS = [1e-4, 1e-3, 1e-2, 1e-1, 0.5, 0.9, 1.0]   # m_surf / m_liquid, ascending


def _geom(fill):
    return TankGeometry(2.0, 6.0, 0.002, "Al6061", fill)


def test_hfg_LH2_sane():
    # Reference (PROJECT.md s.9): ~446 kJ/kg at 1 atm
    assert 440e3 < props.h_fg("Hydrogen", P0) < 455e3


@pytest.mark.parametrize("fluid", FLUIDS)
@pytest.mark.parametrize("fill", FILLS)
def test_bounds_ordered(fluid, fill):
    """Surface bound stores no more energy than homogeneous (reaches P_max no
    later) at every m_surf, stores more as m_surf grows, and pressurises at
    least as fast."""
    g = _geom(fill)
    E_h = stored_energy(UllageClosure.HOMOGENEOUS, fluid, g, P0, ClosedTankSpec(P_MAX))
    E_s = [stored_energy(UllageClosure.SURFACE, fluid, g, P0, ClosedTankSpec(P_MAX, f))
           for f in FRACTIONS]
    for f, E in zip(FRACTIONS, E_s):
        assert E <= E_h * (1.0 + 1e-12), f"{fluid} fill={fill} f={f}: E_s/E_h={E / E_h:.6f}"
    assert all(a < b for a, b in zip(E_s, E_s[1:])), "E_surf must increase with m_surf"
    t = np.array([0.0, 3600.0])
    Q = 50.0
    Ph = pressure_history(UllageClosure.HOMOGENEOUS, fluid, g, P0, Q, t, ClosedTankSpec(P_MAX))
    for f in (1e-2, 0.5):
        Ps = pressure_history(UllageClosure.SURFACE, fluid, g, P0, Q, t, ClosedTankSpec(P_MAX, f))
        assert Ps[1] - Ps[0] >= Ph[1] - Ph[0]


@pytest.mark.parametrize("fluid", FLUIDS)
@pytest.mark.parametrize("fill", FILLS)
def test_bounds_collapse(fluid, fill):
    """m_surf -> m_liquid recovers the homogeneous bound exactly (two-zone form)."""
    g = _geom(fill)
    spec = ClosedTankSpec(P_MAX, 1.0)
    assert warm_zone_initial(UllageClosure.SURFACE, fluid, g, P0, spec) == \
        pytest.approx(homogeneous_initial(fluid, g, P0), rel=1e-12)
    E_h = stored_energy(UllageClosure.HOMOGENEOUS, fluid, g, P0, spec)
    E_s = stored_energy(UllageClosure.SURFACE, fluid, g, P0, spec)
    assert math.isclose(E_s, E_h, rel_tol=1e-10)


@pytest.mark.parametrize("fluid", FLUIDS)
def test_cp_form_violates_ordering_at_high_fill(fluid):
    """Regression record: the textbook c_p surface closure over-stores energy
    at 90 % fill, which is why the two-zone form replaced it (report s.5)."""
    g = _geom(0.9)
    m, _, _ = homogeneous_initial(fluid, g, P0)
    E_h = stored_energy(UllageClosure.HOMOGENEOUS, fluid, g, P0, ClosedTankSpec(P_MAX))
    E_cp = L3.stored_energy_surface(m, lambda T: props.cp_l_sat(fluid, T),
                                    props.T_sat(fluid, P0), props.T_sat(fluid, P_MAX))
    assert E_cp > E_h


@pytest.mark.parametrize("fluid", FLUIDS)
def test_liquid_full_guard(fluid):
    """95 % full: the liquid swells to fill the tank before 3 bar."""
    with pytest.raises(LiquidFullError):
        stored_energy(UllageClosure.HOMOGENEOUS, fluid, _geom(0.95), P0, ClosedTankSpec(P_MAX))


@pytest.mark.parametrize("closure,f", [(UllageClosure.HOMOGENEOUS, 1.0),
                                       (UllageClosure.SURFACE, 0.05)])
def test_energy_balance(closure, f):
    """Transient energy closes to better than 0.5 % (PROJECT.md s.7)."""
    g, spec = _geom(0.9), ClosedTankSpec(P_MAX, f)
    m, rho, u0 = warm_zone_initial(closure, "Hydrogen", g, P0, spec)
    Q, t = 40.0, 2 * 86400.0
    P_end = pressure_history(closure, "Hydrogen", g, P0, Q, np.array([t]), spec)[0]
    u_end = props.u_from_rho_P("Hydrogen", rho, P_end)
    assert abs(m * (u_end - u0) - Q * t) / (Q * t) < 5e-3


def test_boiloff_published():
    pytest.skip("Awaiting MHTB / GFSSP demonstration-tank data in the Stage 1 register")
