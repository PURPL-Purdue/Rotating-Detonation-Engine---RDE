"""Deflection of a Prandtl-Meyer expansion matching a given ambient pressure.

Python translation of the script ``expansion.m`` (src2/) from
xavierdechamps/Nozzle_MOC_supersonic (GPL-3.0); see MODELING/moc/README.md.
"""

from __future__ import annotations

from dataclasses import dataclass

from ._matlab import fzero
from .FinalWaveExpansion import FinalWaveExpansion
from .get_Mach_from_nu import get_Mach_from_nu
from .get_prandtl_meyer_function import get_prandtl_meyer_function


@dataclass
class ExpansionResult:
    """Workspace of ``expansion.m``."""

    p1static: float  # [Pa] static pressure before the expansion
    nu1: float  # [deg] Prandtl-Meyer function before the expansion
    mu1: float  # [deg] Mach angle before the expansion
    M2: float  # Mach number after the *initial guess* deflection (20 deg)
    p2static: float  # [Pa] static pressure after the initial guess deflection
    delta: float  # [deg] deflection matching Pstatic


def expansion(M1=3.518326326701526, theta1=15, gamma=1.2, Pstag=150.0e5, Pstatic=101325):
    """Solve ``FinalWaveExpansion(delta) = 0`` on [0, 20] deg (defaults = MATLAB script values)."""
    p1static = Pstag / (1.0 + 0.5 * (gamma - 1) * M1**2) ** (gamma / (gamma - 1))

    nu1, mu1 = get_prandtl_meyer_function(M1, gamma)

    # Initial values, to launch the iterative process
    delta = 20  # deflection angle, nu2 = nu1 + delta
    M2 = get_Mach_from_nu(nu1 + delta, gamma)
    p2static = Pstag / (1.0 + 0.5 * (gamma - 1) * M2**2) ** (gamma / (gamma - 1))

    FinalWaveExpansion(delta, gamma, nu1, Pstag, Pstatic)

    delta = fzero(
        lambda delta12: FinalWaveExpansion(delta12, gamma, nu1, Pstag, Pstatic), (0.0, float(delta)), tolx=1e-12
    )
    return ExpansionResult(float(p1static), float(nu1), float(mu1), float(M2), float(p2static), float(delta))


def main():
    """Equivalent of running ``expansion.m`` (prints ``delta``)."""
    res = expansion()
    print(f"delta = {res.delta:.4f}")
    return res


if __name__ == "__main__":
    main()
