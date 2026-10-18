"""Pressure mismatch after a Prandtl-Meyer expansion of a given deflection.

Python translation of ``FinalWaveExpansion.m`` from xavierdechamps/Nozzle_MOC_supersonic
(GPL-3.0); see MODELING/moc/README.md.
"""

from .get_Mach_from_nu import get_Mach_from_nu


def FinalWaveExpansion(deflection, gamma, nu1, Pstag, pambiant):
    """Return ``pambiant - p2static`` after expanding from ``nu1`` by ``deflection`` [deg]."""
    M2 = get_Mach_from_nu(nu1 + deflection, gamma)

    p2static = Pstag / (1.0 + 0.5 * (gamma - 1) * M2**2) ** (gamma / (gamma - 1))

    deltaP = pambiant - p2static
    return deltaP
