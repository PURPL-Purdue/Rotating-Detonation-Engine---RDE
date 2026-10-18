"""Ideal-gas choked (sonic) flow relations used by the HADES feed-sizing tools.

These are the equations that every ``SonicNozzleSizingV*.m`` script and the
blowdown Cd scripts wrote out inline.  They are kept local to this package
until the shared ``rdelib`` compressible-flow helpers exist; every function
here is a candidate for consolidation into ``rdelib``.

The expressions keep the operation order of the MATLAB originals so the port
reproduces MATLAB results to round-off.

Choked mass flow through an orifice/throat of area ``A`` fed from stagnation
conditions ``P0, T0``::

    mdot = Cd * A * P0 / sqrt(T0) * K
    K    = sqrt(gamma / R) * (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1)))

The flow is choked when ``P_down / P0 <= (2 / (gamma + 1)) ** (gamma / (gamma - 1))``.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "critical_pressure_ratio",
    "choked_flow_parameter",
    "choked_mass_flow",
    "choked_area",
    "choked_upstream_pressure",
    "area_to_diameter",
]


def critical_pressure_ratio(gamma):
    """Critical (sonic) pressure ratio ``P*/P0 = (2/(gamma+1))**(gamma/(gamma-1))``."""
    return (2 / (gamma + 1)) ** (gamma / (gamma - 1))


def choked_flow_parameter(gamma, R):
    """Choked mass-flow parameter ``K`` [sqrt(K)*s/m] (``mdot = Cd*A*P0/sqrt(T0)*K``).

    Parameters
    ----------
    gamma : float
        Specific heat ratio [-].
    R : float
        Specific gas constant [J/(kg*K)].
    """
    return np.sqrt(gamma / R) * (2 / (gamma + 1)) ** ((gamma + 1) / (2 * (gamma - 1)))


def choked_mass_flow(Cd, A, P0, T0, gamma, R):
    """Choked mass flow [kg/s] through area ``A`` [m^2] from ``P0`` [Pa], ``T0`` [K]."""
    return Cd * A * P0 / np.sqrt(T0) * choked_flow_parameter(gamma, R)


def choked_area(m_dot, Cd, P0, T0, K):
    """Flow area [m^2] that passes ``m_dot`` [kg/s] choked from ``P0`` [Pa], ``T0`` [K].

    ``A = m_dot * sqrt(T0) / (Cd * P0 * K)`` with ``K`` from
    :func:`choked_flow_parameter`.
    """
    return (m_dot * np.sqrt(T0)) / (Cd * P0 * K)


def choked_upstream_pressure(m_dot, Cd, A, T0, K):
    """Stagnation pressure [Pa] needed to pass ``m_dot`` [kg/s] choked through ``A`` [m^2].

    ``P0 = m_dot * sqrt(T0) / (Cd * A * K)`` with ``K`` from
    :func:`choked_flow_parameter`.  Works element-wise on arrays.
    """
    return (m_dot * np.sqrt(T0)) / (Cd * A * K)


def area_to_diameter(A):
    """Equivalent circular diameter ``sqrt(4*A/pi)`` for area ``A``."""
    return np.sqrt(4 * A / np.pi)
