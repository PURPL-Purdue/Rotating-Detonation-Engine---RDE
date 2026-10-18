"""Speed of sound from the energy equation.

Python translation of ``get_speed_sound.m`` from xavierdechamps/Nozzle_MOC_supersonic
(GPL-3.0); see MODELING/moc/README.md.
"""

from ._matlab import msqrt


def get_speed_sound(params, V):
    """Speed of sound a^2 = a0^2 - 0.5*(gamma-1)*V^2 = gamma*R*T0 - 0.5*(gamma-1)*V^2."""
    return msqrt(params.gamma * params.R * params.T - 0.5 * (params.gamma - 1.0) * V**2)
