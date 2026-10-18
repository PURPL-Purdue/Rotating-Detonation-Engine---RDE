"""Direct wall point: intersection of a C+ characteristic with the nozzle wall.

Python translation of ``MOC_2D_steady_irrotational_wall.m`` from
xavierdechamps/Nozzle_MOC_supersonic (GPL-3.0); see MODELING/moc/README.md.
"""

import numpy as np

from ._errors import MocConvergenceError
from ._matlab import narrow, asind, atand, mmax, msqrt, tand
from .get_speed_sound import get_speed_sound
from .MOC_2D_steady_irrotational_get_geometry import MOC_2D_steady_irrotational_get_geometry


def MOC_2D_steady_irrotational_wall(xp, yp, up, vp, geom, params):
    """Return ``(xo, yo, uo, vo)`` where the C+ characteristic from point p meets the wall."""
    xp, yp, up, vp = narrow(xp), narrow(yp), narrow(up), narrow(vp)
    x_orig = xp
    y_orig = yp
    u_orig = up
    v_orig = vp

    step_current = 0
    step_max = 50
    eps_pos = 1.0e-8
    eps_vel = 1.0e-5

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        # Predictor step
        xo, yo, uo, vo = _solve_wall(xp, yp, up, vp, x_orig, y_orig, u_orig, v_orig, geom, params)
        while True:
            step_current += 1
            # Corrector step
            xcp = xp
            ycp = 0.5 * (yp + yo)
            ucp = 0.5 * (up + uo)
            vcp = 0.5 * (vp + vo)
            xn, yn, un, vn = _solve_wall(xcp, ycp, ucp, vcp, x_orig, y_orig, u_orig, v_orig, geom, params)
            error_pos = mmax(xo - xn, yo - yn)
            error_vel = mmax(uo - un, vo - vn)
            xo = xn
            yo = yn
            uo = un
            vo = vn
            if abs(error_pos) < eps_pos and abs(error_vel) < eps_vel:
                break
            if step_current > step_max:
                raise MocConvergenceError(
                    "The maximum of iterations for the predictor-corrector algorithm has been "
                    "reached. Stopping the execution..."
                )
    return xo, yo, uo, vo


def _solve_wall(xp, yp, up, vp, x_orig, y_orig, u_orig, v_orig, geom, params):
    Vmag = msqrt(up**2 + vp**2)
    a = get_speed_sound(params, Vmag)
    theta = atand(vp / up)
    alpha = asind(a / Vmag)
    lam = tand(theta + alpha)  # lambda+
    Q = up**2 - a**2
    R = 2 * up * vp - Q * lam
    S = geom.delta * a**2 * vp / yp

    # Coefficients a, b, c of the 2nd order poly for the wall
    abcpoly, _, _ = MOC_2D_steady_irrotational_get_geometry(0.0, 1, geom)

    # C+ characteristic: y4 - lambda+ * x4 = y_orig - lambda+ * x_orig
    # -> second order equation a*x4^2 + b*x4 + c = 0
    a4 = abcpoly[0] - (y_orig - lam * x_orig)
    b4 = abcpoly[1] - lam
    c4 = abcpoly[2]
    if c4 == 0:
        xn = -a4 / b4  # Intersection between 2 lines
    else:
        xn = (-b4 - msqrt(b4**2 - 4 * a4 * c4)) / (2 * c4)  # line / parabola intersection

    _, yn, tangentwall = MOC_2D_steady_irrotational_get_geometry(xn, 2, geom)

    T = S * (xn - x_orig) + Q * u_orig + R * v_orig
    un = T / (Q + R * tangentwall)
    vn = un * tangentwall
    return xn, yn, un, vn
