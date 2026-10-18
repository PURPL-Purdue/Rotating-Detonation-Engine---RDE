"""Transonic initial-value line (Sauer/Hall small-perturbation solution at the throat).

Python translation of ``MOC_2D_steady_irrotational_IVLINE.m`` from
xavierdechamps/Nozzle_MOC_supersonic (GPL-3.0); see MODELING/moc/README.md.
"""

import numpy as np


def MOC_2D_steady_irrotational_IVLINE(geom, params, y):
    """Return ``(xsonic, xvnull, u)`` along the throat ordinates ``y`` [m].

    ``xsonic``: abscissae of the sonic line (M = 1); ``xvnull``: abscissae of the
    v = 0 line used as the initial-value line; ``u``: axial velocity on it.
    """
    y = np.asarray(y, dtype=float)

    # Origin of the coordinate system relative to the nozzle throat [m]
    eps = (
        -geom.yt
        * np.sqrt((params.gamma + 1) * (geom.delta + 1) * geom.yt / geom.rhou)
        * 0.5
        / (3 + geom.delta)
    )

    # Sonic speed a*
    astar = np.sqrt(2 * params.gamma * params.R * params.T / (params.gamma + 1))

    # Coefficient of the linear axial perturbation velocity [1/m]
    alpha = np.sqrt((1 + geom.delta) / (geom.yt * geom.rhou * (params.gamma + 1)))

    # Equation for the sonic line x = coeff_sonic * y^2, where coeff_sonic is [1/m]
    coeff = -(params.gamma + 1) * alpha * 0.5 / (1 + geom.delta)
    xsonic = coeff * y**2 - eps

    # Equation for the v=0 line: x = coeff_vnull * y^2, where coeff_vnull is [1/m]
    coeff_vnull = -(params.gamma + 1) * alpha * 0.5 / (3 + geom.delta)
    xvnull = coeff_vnull * y**2 - eps

    # Perturbation velocity field on the v=0 line
    u = astar * (
        1
        + alpha * (xvnull + eps)
        + (params.gamma + 1) * (alpha**2) * (y**2) * 0.5 / (1 + geom.delta)
    )
    return xsonic, xvnull, u
