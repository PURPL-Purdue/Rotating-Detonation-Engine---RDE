"""1-D isentropic (choked de Laval nozzle) Mach number and pressure.

Python translation of ``get_Laval_theory.m`` from xavierdechamps/Nozzle_MOC_supersonic
(GPL-3.0); see MODELING/moc/README.md.

Porting fix (see PORTING_NOTES.md): the area ratio is ``(Y/Y(1))^(1+delta)``.
The MATLAB always used ``Y.^2/Y(1)^2`` (axisymmetric), which is wrong for the
planar nozzle (``geom.delta = 0``) that the team's driver uses.
"""

import numpy as np

from ._matlab import fzero
from .structs import Chocked


def get_Laval_theory(geom, params, X, Y):
    """Return a :class:`Chocked` with the isentropic ``mach`` and ``pressure`` (p/P0).

    ``X``, ``Y`` are the wall-node coordinates downstream of the throat; the area
    ratio is referenced to the first node ``Y[0]``.
    """
    gam = params.gamma
    gamm1 = gam - 1
    gamp1 = gam + 1

    Y = np.ravel(np.asarray(Y, dtype=float))
    if geom.delta == 1:
        area = Y**2 / (Y[0] ** 2)  # axisymmetric (unchanged MATLAB expression)
    else:
        area = Y / Y[0]  # planar: A/A* = y/y*  (MATLAB bug fix)

    # Equation linking the area ratio and the Mach number ratio
    # m = M/M* and a = A/A* where the star denotes conditions at the throat
    def the_eq(m, a):
        return (2 * (1 + 0.5 * gamm1 * m**2) / gamp1) ** (0.5 * gamp1 / gamm1) / m - a

    chocked = Chocked()
    chocked.mach = np.zeros(area.size)
    for i in range(area.size):
        a_i = float(area[i])
        chocked.mach[i] = fzero(lambda m: the_eq(m, a_i), (1.0, 25.0))

    # Pressure is the ratio p/P_stag
    chocked.pressure = (1 + 0.5 * gamm1 * chocked.mach**2) ** (-gam / gamm1)
    return chocked
