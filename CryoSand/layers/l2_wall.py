"""L2 - Tank wall.

Lumped plane-wall resistance for now. The 2-D ThermoSand FDM solver
(harmonic-mean faces, Dirichlet/Neumann/Robin BCs) is ported here in build
step 8; this lumped form is its 1-D limit and the check it must reproduce.

Expected result: negligible. For 2 mm aluminium over ~100 m^2 the wall
resistance is ~1e-7 of the blanket resistance. We compute it rather than
assume it (PROJECT.md s.4, L2).
"""
from __future__ import annotations

import math

from cryosand.core.network import R_slab


def R_wall_plane(thickness: float, k_wall: float, A: float) -> float:
    """Thin-wall plane approximation R = t/(k A) [K/W]."""
    return R_slab(thickness, k_wall, A)


def R_wall_cylinder(r_inner: float, thickness: float, length: float, k_wall: float) -> float:
    """Exact cylindrical-shell resistance R = ln(r_o/r_i) / (2 pi k L) [K/W].
    Tends to the plane form as thickness/r_inner -> 0 (tested)."""
    r_o = r_inner + thickness
    return math.log(r_o / r_inner) / (2.0 * math.pi * k_wall * length)
