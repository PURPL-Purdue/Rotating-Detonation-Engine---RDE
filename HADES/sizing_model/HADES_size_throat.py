"""Port of ``HADES_size_throat.m`` (minimum choked throat area)."""

from __future__ import annotations

import math

__all__ = ["HADES_size_throat"]


def HADES_size_throat(mdot, P0, T0, gamma, R):
    """Minimum throat area [m^2] for a choked RDE nozzle.

    Parameters
    ----------
    mdot : mass flow rate [kg/s]
    P0 : stagnation pressure [Pa]
    T0 : stagnation temperature [K]
    gamma : ratio of specific heats [-]
    R : specific gas constant [J/(kg K)]
    """
    term1 = gamma / R
    term2 = (2 / (gamma + 1)) ** ((gamma + 1) / (gamma - 1))
    flow_factor = math.sqrt(term1 * term2)
    At = (mdot * math.sqrt(T0)) / (P0 * flow_factor)
    return At
