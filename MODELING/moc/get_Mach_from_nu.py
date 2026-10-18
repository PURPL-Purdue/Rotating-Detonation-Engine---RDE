"""Inverse Prandtl-Meyer function.

Python translation of ``get_Mach_from_nu.m`` from xavierdechamps/Nozzle_MOC_supersonic
(GPL-3.0); see MODELING/moc/README.md.
"""

import math

import numpy as np

from ._matlab import fzero

_R2D = 180.0 / math.pi


def get_Mach_from_nu(nu, gamma):
    """Return the Mach number whose Prandtl-Meyer function equals ``nu`` [deg].

    Solved with MATLAB's bracketed ``fzero`` on [1.00001, 50].  Scalar input
    returns a float; array input returns a 1-D array (MATLAB row vector).
    """
    sqroot = math.sqrt((gamma + 1) / (gamma - 1))

    def equation(Mach, nu_k):
        root = math.sqrt(Mach * Mach - 1)
        return sqroot * (_R2D * math.atan(root / sqroot)) - _R2D * math.atan(root) - nu_k

    guessMach = (1.00001, 50.0)
    nu_arr = np.ravel(np.asarray(nu, dtype=float))
    ma = np.zeros(nu_arr.size)
    for k in range(nu_arr.size):
        if not nu_arr[k] > 0:
            raise AssertionError("assert(nu(k)>0) failed")
        if not nu_arr[k] <= 125.0:
            raise AssertionError("assert(nu(k)<=125.) failed")
        nu_k = float(nu_arr[k])
        ma[k] = fzero(lambda Mach: equation(Mach, nu_k), guessMach)

    if np.ndim(nu) == 0:
        return float(ma[0])
    return ma
