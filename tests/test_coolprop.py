"""Tests that need CoolProp. Skipped automatically where it is not installed.

test_bounds_ordered protects the entire research claim (PROJECT.md s.7).
"""
import math

import numpy as np
import pytest

pytest.importorskip("CoolProp")

from cryosand.core import properties as props  # noqa: E402
from cryosand.core.types import TankGeometry, UllageClosure  # noqa: E402
from cryosand.model import ClosedTankSpec, homogeneous_initial, pressure_history, stored_energy  # noqa: E402

P0 = 101325.0
P_MAX = 3.0e5
FLUIDS = ["Hydrogen", "Oxygen", "Methane"]
FILLS = [0.3, 0.6, 0.9]
FRACTIONS = [1e-4, 1e-3, 1e-2, 1e-1, 0.5, 1.0]


def _geom(fill):
    return TankGeometry(2.0, 6.0, 0.002, "Al6061", fill)


def test_hfg_LH2_sane():
    # Reference (PROJECT.md s.9): ~446 kJ/kg at 1 atm
    assert 440e3 < props.h_fg("Hydrogen", P0) < 455e3


@pytest.mark.parametrize("fluid", FLUIDS)
@pytest.mark.parametrize("fill", FILLS)
def test_bounds_ordered(fluid, fill):
    """Surface bound must store no more energy than homogeneous => reaches P_max
    no later, at every m_surf <= m_total. Also compares initial dP/dt."""
    g = _geom(fill)
    E_h = stored_energy(UllageClosure.HOMOGENEOUS, fluid, g, P0, ClosedTankSpec(P_MAX))
    for f in FRACTIONS:
        E_s = stored_energy(UllageClosure.SURFACE, fluid, g, P0, ClosedTankSpec(P_MAX, f))
        assert E_s <= E_h, f"{fluid} fill={fill} f={f}: E_s/E_h={E_s / E_h:.4f}"
    t = np.array([0.0, 3600.0])
    Q = 50.0
    Ph = pressure_history(UllageClosure.HOMOGENEOUS, fluid, g, P0, Q, t, ClosedTankSpec(P_MAX))
    Ps = pressure_history(UllageClosure.SURFACE, fluid, g, P0, Q, t, ClosedTankSpec(P_MAX, 1.0))
    assert (Ps[1] - P0) >= (Ph[1] - Ph[0])


@pytest.mark.parametrize("fluid", FLUIDS)
def test_bounds_collapse(fluid):
    """m_surf -> m_total must recover the homogeneous bound. EXPECTED TO
    EXPOSE the cp-vs-constant-volume-path inconsistency (report s.5.2); the
    tolerance is provisional pending the first CoolProp run."""
    g = _geom(0.9)
    E_h = stored_energy(UllageClosure.HOMOGENEOUS, fluid, g, P0, ClosedTankSpec(P_MAX))
    E_s = stored_energy(UllageClosure.SURFACE, fluid, g, P0, ClosedTankSpec(P_MAX, 1.0))
    print(f"{fluid}: E_surf(m_total)/E_hom = {E_s / E_h:.4f}")
    assert math.isclose(E_s, E_h, rel_tol=0.05)


def test_energy_balance_homogeneous():
    g = _geom(0.9)
    m, rho, u0 = homogeneous_initial("Hydrogen", g, P0)
    Q, t = 40.0, 30 * 86400.0
    P_end = pressure_history(UllageClosure.HOMOGENEOUS, "Hydrogen", g, P0, Q, np.array([t]),
                             ClosedTankSpec(P_MAX))[0]
    u_end = props._props("U", "D", rho, "P", P_end, "Hydrogen")
    assert abs(m * (u_end - u0) - Q * t) / (Q * t) < 5e-3


def test_boiloff_published():
    pytest.skip("Awaiting MHTB / GFSSP demonstration-tank data in the Stage 1 register")
