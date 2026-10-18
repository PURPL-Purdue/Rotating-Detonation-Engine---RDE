"""Prandtl-Meyer function and Mach angle.

Python translation of ``get_prandtl_meyer_function.m`` from
xavierdechamps/Nozzle_MOC_supersonic (GPL-3.0); see MODELING/moc/README.md.
"""

import numpy as np

from ._matlab import asind, atand, msqrt


def get_prandtl_meyer_function(Mach, gamma):
    """Return ``(nu, mu)`` [deg]: Prandtl-Meyer function and Mach angle of ``Mach``."""
    if not np.all(np.real(np.asarray(Mach)) > 0):
        raise AssertionError("assert(Mach>0) failed")
    mu = asind(1.0 / np.asarray(Mach, dtype=float)) if np.ndim(Mach) else asind(1.0 / float(Mach))

    sqroot = np.sqrt((gamma + 1) / (gamma - 1))
    if np.ndim(Mach) == 0:
        root = msqrt(float(Mach) * float(Mach) - 1)
    else:
        M = np.asarray(Mach, dtype=float)
        root = msqrt(M * M - 1)
    nu = sqroot * atand(root / sqroot) - atand(root)
    return nu, mu
