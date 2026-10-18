"""Nozzle wall contour: downstream circular arc + 2nd-order polynomial.

Python translation of ``MOC_2D_steady_irrotational_get_geometry.m`` from
xavierdechamps/Nozzle_MOC_supersonic (GPL-3.0); see MODELING/moc/README.md.
"""

import numpy as np

from ._matlab import cosd, mge, sind, tand


def MOC_2D_steady_irrotational_get_geometry(x, step, geom):
    """Return ``(abc, y, tangent)`` of the wall contour (Gas Dynamics V2, p.128).

    The divergent region has a circular arc of radius ``rhod`` up to ``ta`` deg;
    at the attachment point a parabola ``y = a + b x + c x^2`` (coefficients
    ``abc``) is tangent to the arc and reaches ``xe`` with exit angle ``te``.
    With ``step == 2`` the wall ordinate ``y`` and slope ``tangent`` at ``x``
    (which must lie downstream of the attachment point) are also returned;
    otherwise both are 0.
    """
    radtd = geom.rhod  # Radius of curvature of downstream circular arc
    yt = geom.yt  # Radius of nozzle throat
    thetaa = geom.ta  # Attachment angle between downstream arc and 2nd order poly
    thetae = geom.te  # Nozzle exit lip angle
    xe = geom.xe  # Nozzle length

    # Coordinates of attachment point between downstream arc and 2nd order poly
    xa = radtd * sind(thetaa)
    ya = yt + (1 - cosd(thetaa)) * radtd

    # 1st line: (x,y) at A known; 2nd line: slope at A = tan(ta); 3rd: slope at xe = tan(te)
    matA = np.array([[1.0, xa, xa**2], [0.0, 1.0, 2 * xa], [0.0, 1.0, 2 * xe]])
    RHS = np.array([ya, tand(thetaa), tand(thetae)])
    abc = np.linalg.solve(matA, RHS)

    y = 0.0
    tangent = 0.0
    if step == 2:  # geometry after the downstream circular arc
        if not mge(x, xa):
            raise AssertionError("assert(x >= xa) failed: point upstream of the attachment point")
        y = abc[0] + abc[1] * x + abc[2] * x**2
        tangent = abc[1] + 2 * abc[2] * x
    return abc, y, tangent
