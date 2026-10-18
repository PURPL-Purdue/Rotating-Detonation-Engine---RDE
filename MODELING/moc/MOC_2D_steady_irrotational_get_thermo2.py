"""Thermodynamic state from the Mach number.

Python translation of ``MOC_2D_steady_irrotational_get_thermo2.m`` from
xavierdechamps/Nozzle_MOC_supersonic (GPL-3.0); see MODELING/moc/README.md.
"""

import numpy as np

from ._matlab import msqrt


def MOC_2D_steady_irrotational_get_thermo2(Mach, params):
    """Return ``(pressure, density, temperature, sound)`` for isentropic flow at ``Mach``."""
    # Static temperature
    temperature = params.T / (1 + 0.5 * (params.gamma - 1) * Mach**2)

    # Local acoustic speed
    sound = msqrt(params.gamma * params.R * temperature)

    # Static pressure
    pressure = params.P / (1 + 0.5 * (params.gamma - 1) * Mach**2) ** (
        params.gamma / (params.gamma - 1)
    )

    # Density
    density = pressure / (params.R * temperature)
    if np.ndim(Mach) == 0:
        return float(pressure), float(density), float(temperature), sound
    return pressure, density, temperature, sound
