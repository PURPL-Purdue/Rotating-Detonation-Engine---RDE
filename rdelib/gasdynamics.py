"""Compressible-flow and pipe-flow relations (SI units throughout).

Perfect-gas relations used by the HADES sizing model, the feed-system
sizing and DEIMOS. Every function works on scalars and on numpy arrays
unless it solves for a root.

Several ports keep their own copies of algebraically identical relations
(e.g. the choked-flow factor in ``HADES.feed_sizing.choked_flow``) because
their operation order differs and the ports are pinned bit-for-bit to the
MATLAB/Octave results. Only expressions with the same operation order live
here; see the PORTING_NOTES.md of each package.

Conventions
-----------
* ``gamma``  ratio of specific heats [-]
* ``R``      specific gas constant [J/(kg K)]
* ``M``      Mach number [-]
* pressures in Pa, temperatures in K, areas in m^2, mass flow in kg/s
* friction factors are **Darcy** factors (f = 4 * Fanning)
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq

__all__ = [
    "T0_over_T",
    "P0_over_P",
    "rho0_over_rho",
    "area_ratio",
    "critical_pressure_ratio",
    "choked_flow_factor",
    "choked_throat_area",
    "choked_mass_flow",
    "subcritical_orifice_mass_flow",
    "incompressible_orifice_area",
    "fanno_fLD",
    "fanno_p_over_pstar",
    "fanno_T_over_Tstar",
    "fanno_fLD_max",
    "fanno_mach_from_fLD",
    "colebrook",
    "colebrook_darcy",
    "tkachenko_friction_factor",
    "sutherland_viscosity",
    "reynolds_number",
]


# --- Isentropic relations ---------------------------------------------------
def T0_over_T(M, gamma):
    """Stagnation-to-static temperature ratio, 1 + (gamma-1)/2 M^2."""
    return 1.0 + 0.5 * (gamma - 1.0) * np.asarray(M) ** 2


def P0_over_P(M, gamma):
    """Stagnation-to-static pressure ratio."""
    return T0_over_T(M, gamma) ** (gamma / (gamma - 1.0))


def rho0_over_rho(M, gamma):
    """Stagnation-to-static density ratio."""
    return T0_over_T(M, gamma) ** (1.0 / (gamma - 1.0))


def area_ratio(M, gamma):
    """Isentropic area ratio A/A* for Mach number ``M``."""
    M = np.asarray(M, dtype=float)
    return (1.0 / M) * ((2.0 / (gamma + 1.0)) * T0_over_T(M, gamma)) ** (
        (gamma + 1.0) / (2.0 * (gamma - 1.0))
    )


def critical_pressure_ratio(gamma):
    """Critical (sonic) pressure ratio p*/p0 = (2/(gamma+1))^(gamma/(gamma-1))."""
    return (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))


# --- Choked / orifice flow --------------------------------------------------
def choked_flow_factor(gamma, R):
    """sqrt(gamma/R * (2/(gamma+1))^((gamma+1)/(gamma-1))) [s sqrt(K)/m].

    Choked mass flux is ``mdot/A = P0 * factor / sqrt(T0)``.
    """
    return math.sqrt(gamma / R * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (gamma - 1.0)))


def choked_throat_area(mdot, P0, T0, gamma, R):
    """Minimum (sonic) throat area [m^2] that passes ``mdot`` [kg/s].

    Same formula as ``HADES_size_throat``: ``A* = mdot sqrt(T0) / (P0 * F)``
    with ``F = choked_flow_factor(gamma, R)``. P0 in Pa, T0 in K.
    """
    return mdot * np.sqrt(T0) / (P0 * choked_flow_factor(gamma, R))


def choked_mass_flow(A, P0, T0, gamma, R, Cd=1.0):
    """Choked mass flow [kg/s] through area ``A`` [m^2] with discharge coefficient ``Cd``."""
    return Cd * A * P0 * choked_flow_factor(gamma, R) / np.sqrt(T0)


def subcritical_orifice_mass_flow(Cd, A, P0, P, T, gamma, R):
    """Compressible (unchoked) orifice mass flow [kg/s].

    ``mdot = Cd A P0 sqrt( 2 gamma/(R T (gamma-1)) (P/P0)^(2/gamma) (1-(P/P0)^((gamma-1)/gamma)) )``

    This is the form used in ``HADES_size_P_plenum_fanno`` (``T`` is the
    temperature used in that script, i.e. the static inlet temperature).
    """
    r = P / P0
    return Cd * A * P0 * np.sqrt(
        (2.0 * gamma / (R * T * (gamma - 1.0))) * r ** (2.0 / gamma) * (1.0 - r ** ((gamma - 1.0) / gamma))
    )


def incompressible_orifice_area(mdot, Cd, rho, dP):
    """Orifice area [m^2] from the incompressible orifice equation.

    ``A = mdot / (Cd sqrt(2 rho dP))`` with mdot in kg/s, rho in kg/m^3, dP in Pa.
    """
    return mdot / (Cd * np.sqrt(2.0 * rho * dP))


# --- Fanno flow -------------------------------------------------------------
def fanno_fLD(M, gamma):
    """Fanno parameter f L*/D (Darcy f) from Mach ``M`` to the choking point.

    ``(1-M^2)/(gamma M^2) + (gamma+1)/(2 gamma) ln( (gamma+1) M^2 / (2 + (gamma-1) M^2) )``
    """
    M2 = np.asarray(M, dtype=float) ** 2
    return (1.0 - M2) / (gamma * M2) + (gamma + 1.0) / (2.0 * gamma) * np.log(
        (gamma + 1.0) * M2 / (2.0 + (gamma - 1.0) * M2)
    )


def fanno_p_over_pstar(M, gamma):
    """Fanno static pressure ratio p/p* = (1/M) sqrt((gamma+1)/(2+(gamma-1)M^2))."""
    M = np.asarray(M, dtype=float)
    return (1.0 / M) * np.sqrt((gamma + 1.0) / (2.0 + (gamma - 1.0) * M**2))


def fanno_T_over_Tstar(M, gamma):
    """Fanno static-temperature ratio T/T* = (gamma+1) / (2 + (gamma-1) M^2)."""
    return (gamma + 1) / (2 + (gamma - 1) * M**2)


def fanno_fLD_max(gamma: float) -> float:
    """Limit of :func:`fanno_fLD` as M -> infinity (largest value on the supersonic branch)."""
    return ((gamma + 1) / (2 * gamma)) * math.log((gamma + 1) / (gamma - 1)) - 1 / gamma


def fanno_mach_from_fLD(fLD, gamma, branch: str = "subsonic"):
    """Invert :func:`fanno_fLD` for the inlet Mach number.

    ``branch='subsonic'`` searches (1e-6, 1); ``'supersonic'`` searches (1, 50).
    """
    g = lambda M: fanno_fLD(M, gamma) - fLD  # noqa: E731
    if branch == "subsonic":
        return brentq(g, 1e-6, 1.0 - 1e-12, xtol=1e-15, rtol=4 * np.finfo(float).eps)
    if branch == "supersonic":
        return brentq(g, 1.0 + 1e-12, 50.0, xtol=1e-15, rtol=4 * np.finfo(float).eps)
    raise ValueError("branch must be 'subsonic' or 'supersonic'")


# --- Friction factor / viscosity -----------------------------------------
_BRENT_RTOL = 4 * np.finfo(float).eps  # brentq minimum


def colebrook_darcy(Re: float, epsilon: float, D: float) -> float:
    """Darcy friction factor from the Colebrook-White equation.

    Solves ``1/sqrt(f) + 2 log10( epsilon/(3.7 D) + 2.51/(Re sqrt(f)) ) = 0``.
    ``epsilon`` and ``D`` must share units. Solved in x = 1/sqrt(f), where the
    residual x + 2 log10(eps/(3.7 D) + 2.51 x/Re) is strictly increasing, so the
    root is unique and bracketed by [1e-6, 1e3] for any eps/D < 3.7. The
    bracketed solve needs no start point and, unlike ``fsolve`` from f = 0.02,
    never stops early with a "not making good progress" warning.
    """
    a = epsilon / (3.7 * D)
    b = 2.51 / Re
    x = brentq(lambda x: x + 2 * math.log10(a + b * x), 1e-6, 1e3, xtol=1e-15, rtol=_BRENT_RTOL, maxiter=500)
    return 1 / x**2


def colebrook(Re: float, rel_rough: float, f0: float | None = None) -> float:
    """Darcy friction factor from Colebrook-White for relative roughness ``rel_rough`` = epsilon/D.

    Same bracketed solve as :func:`colebrook_darcy` (``colebrook_darcy(Re, rel_rough, 1.0)``).
    ``f0`` (the MATLAB ``fsolve`` start point, 0.02) is accepted for compatibility
    and ignored: the root is unique.
    """
    return colebrook_darcy(Re, rel_rough, 1.0)


def tkachenko_friction_factor(Re, rel_rough):
    """Explicit Darcy friction factor (Tkachenko & Mileikovskyi 2020 form).

    Exactly the expression used in ``HADES_size_pdrop_fanno``.
    """
    A0 = -0.79638 * np.log(rel_rough / 8.208 + 7.3357 / Re)
    A1 = Re * rel_rough + 9.3120665 * A0
    num = 8.128943 + A1
    den = 8.128943 * A0 - 0.86859209 * A1 * np.log(A1 / (3.7099535 * Re))
    return (num / den) ** 2


def sutherland_viscosity(T, mu_ref, T_ref, S):
    """Sutherland's law, mu = mu_ref (T_ref + S)/(T + S) (T/T_ref)^1.5 [Pa s]."""
    return mu_ref * (T_ref + S) / (T + S) * (T / T_ref) ** 1.5


def reynolds_number(rho, U, D, mu):
    """Reynolds number rho U D / mu."""
    return rho * U * D / mu
