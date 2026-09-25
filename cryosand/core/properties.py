"""Single entry point for fluid properties. Wraps CoolProp (Bell et al. 2014,
Ind. Eng. Chem. Res. 53(6), 2498-2508).

Rules (PROJECT.md s.3):
  * CoolProp is the only source of fluid properties. Nothing is hard-coded.
  * CoolProp returns SI (Pa, K, J/kg, kg/m^3), so no conversion happens here.
  * Calls are memoised; all arguments are hashable floats/strings.

Physics functions in cryosand.layers never import this module. They take
property values as arguments, which keeps them pure and testable without
CoolProp. model.py and trade.py are the only callers.
"""
from __future__ import annotations

from functools import lru_cache

try:
    from CoolProp.CoolProp import PropsSI as _PropsSI
except ImportError as exc:  # pragma: no cover - environment dependent
    _PropsSI = None
    _IMPORT_ERROR = exc
else:
    _IMPORT_ERROR = None


def coolprop_available() -> bool:
    return _PropsSI is not None


def _props(output: str, k1: str, v1: float, k2: str, v2: float, fluid: str) -> float:
    if _PropsSI is None:
        raise ImportError(
            "CoolProp is required for fluid properties: pip install CoolProp"
        ) from _IMPORT_ERROR
    return float(_PropsSI(output, k1, v1, k2, v2, fluid))


@lru_cache(maxsize=4096)
def T_sat(fluid: str, P: float) -> float:
    """Saturation temperature [K] at pressure P [Pa]."""
    return _props("T", "P", P, "Q", 0.0, fluid)


@lru_cache(maxsize=4096)
def P_sat(fluid: str, T: float) -> float:
    """Saturation pressure [Pa] at temperature T [K]."""
    return _props("P", "T", T, "Q", 0.0, fluid)


@lru_cache(maxsize=4096)
def h_fg(fluid: str, P: float) -> float:
    """Latent heat of vaporisation [J/kg] at pressure P [Pa]."""
    return _props("H", "P", P, "Q", 1.0, fluid) - _props("H", "P", P, "Q", 0.0, fluid)


@lru_cache(maxsize=4096)
def rho_l_sat(fluid: str, P: float) -> float:
    """Saturated-liquid density [kg/m^3]."""
    return _props("D", "P", P, "Q", 0.0, fluid)


@lru_cache(maxsize=4096)
def rho_v_sat(fluid: str, P: float) -> float:
    """Saturated-vapour density [kg/m^3]."""
    return _props("D", "P", P, "Q", 1.0, fluid)


@lru_cache(maxsize=4096)
def u_l_sat(fluid: str, P: float) -> float:
    """Saturated-liquid specific internal energy [J/kg]."""
    return _props("U", "P", P, "Q", 0.0, fluid)


@lru_cache(maxsize=4096)
def u_v_sat(fluid: str, P: float) -> float:
    """Saturated-vapour specific internal energy [J/kg]."""
    return _props("U", "P", P, "Q", 1.0, fluid)


@lru_cache(maxsize=4096)
def cp_l_sat(fluid: str, T: float) -> float:
    """Saturated-liquid isobaric specific heat [J/(kg K)] at temperature T."""
    return _props("C", "T", T, "Q", 0.0, fluid)


def P_from_rho_u(fluid: str, rho: float, u: float) -> float:
    """Pressure [Pa] from density [kg/m^3] and specific internal energy [J/kg].
    Valid in the two-phase dome (CoolProp flashes D-U). Not cached: continuous
    state variable in the ODE."""
    return _props("P", "D", rho, "U", u, fluid)


def T_from_rho_u(fluid: str, rho: float, u: float) -> float:
    """Temperature [K] from density and specific internal energy."""
    return _props("T", "D", rho, "U", u, fluid)


@lru_cache(maxsize=4096)
def k_l_sat(fluid: str, P: float) -> float:
    """Saturated-liquid thermal conductivity [W/(m K)]."""
    return _props("L", "P", P, "Q", 0.0, fluid)
