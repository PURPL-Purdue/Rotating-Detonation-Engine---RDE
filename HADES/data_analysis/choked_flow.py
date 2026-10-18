"""Choked (sonic) orifice / nozzle relations shared by the data-analysis scripts.

All scripts in this folder use the same isentropic choked mass-flow equation in
one of three equivalent forms::

    mdot = Cd * A * sqrt(gamma * rho0 * P0 * K)                 (PURPLE_MassFlow_..., cD_analysis)
         = Cd * A * P0 * MFF,      MFF = sqrt(gamma * K / (R*T0))
         = Cd * A * P0 * MFP / sqrt(R*T0)                       (CDRConditions, UpstreamConditions)

    K   = (2/(gamma+1))^((gamma+1)/(gamma-1))
    MFP = sqrt(gamma) * (2/(gamma+1))^((gamma+1)/(2*(gamma-1)))

Candidates for later consolidation into ``rdelib`` (compressible-flow helpers).
"""

from __future__ import annotations

import numpy as np


def choked_flow_constant(gamma):
    """``K = (2/(gamma+1))^((gamma+1)/(gamma-1))``."""
    return (2 / (gamma + 1)) ** ((gamma + 1) / (gamma - 1))


def mass_flow_factor(gamma, R, T0):
    """``MFF = sqrt(gamma*K/(R*T0))`` so that ``mdot = Cd*A*P0*MFF`` [s/m]."""
    return np.sqrt(gamma * choked_flow_constant(gamma) / (R * T0))


def mass_flow_parameter(gamma):
    """Dimensionless ``MFP = sqrt(gamma)*(2/(gamma+1))^((gamma+1)/(2(gamma-1)))``."""
    return np.sqrt(gamma) * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))


def choked_mass_flow(Cd, A, gamma, rho0, P0):
    """``mdot = Cd*A*sqrt(gamma*rho0*P0*K)`` [kg/s] (SI inputs)."""
    return Cd * A * np.sqrt(gamma * rho0 * P0 * choked_flow_constant(gamma))


def stagnation_pressure_for_mass_flow(mdot, Cd, A, gamma, R, T0):
    """Upstream stagnation pressure ``P0 = mdot/(Cd*A*MFF)`` [Pa] that chokes ``mdot``."""
    return mdot / (Cd * A * mass_flow_factor(gamma, R, T0))
