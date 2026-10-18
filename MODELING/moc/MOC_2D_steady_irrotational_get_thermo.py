"""Thermodynamic state from the velocity components.

Python translation of ``MOC_2D_steady_irrotational_get_thermo.m`` from
xavierdechamps/Nozzle_MOC_supersonic (GPL-3.0); see MODELING/moc/README.md.
"""

import numpy as np

from ._matlab import msqrt


def MOC_2D_steady_irrotational_get_thermo(u, v, params):
    """Return ``(Mach, pressure, density, temperature, sound)`` from ``u``, ``v`` [m/s]."""
    u = np.asarray(u)
    v = np.asarray(v)
    with np.errstate(divide="ignore", invalid="ignore"):
        # Square of velocity amplitude
        Q = u**2 + v**2

        # Stagnation acoustic speed
        sound_stag = np.sqrt(params.gamma * params.R * params.T)

        # Local acoustic speed
        sound = msqrt(sound_stag**2 - 0.5 * (params.gamma - 1) * Q)

        # Mach number
        Mach = msqrt(Q) / sound

        # Static temperature
        temperature = params.T / (1 + 0.5 * (params.gamma - 1) * Mach**2)

        # Static pressure
        pressure = params.P / (1 + 0.5 * (params.gamma - 1) * Mach**2) ** (
            params.gamma / (params.gamma - 1)
        )

        # Density
        density = pressure / (params.R * temperature)
    return Mach, pressure, density, temperature, sound
