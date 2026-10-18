"""Free-pressure (jet boundary) point.

Python translation of ``MOC_2D_steady_irrotational_free_pressure.m`` from
xavierdechamps/Nozzle_MOC_supersonic (GPL-3.0); see MODELING/moc/README.md.
"""

import numpy as np

from ._errors import MocConvergenceError
from ._matlab import narrow, asind, atand, mldivide2, mmax, msqrt, real_part, tand
from .get_speed_sound import get_speed_sound


def MOC_2D_steady_irrotational_free_pressure(xw, yw, uw, vw, xp, yp, up, vp, geom, params):
    """Intersection of a C+ characteristic with the plume boundary (p = ``params.Pstatic``).

    ``xw..vw``: point 3 on the jet boundary (previous line);
    ``xp..vp``: point 2, internal point of the previous line.
    Returns ``(xo, yo, uo, vo)`` of the new jet-boundary point.
    """
    xw, yw, uw, vw = narrow(xw), narrow(yw), narrow(uw), narrow(vw)
    xp, yp, up, vp = narrow(xp), narrow(yp), narrow(up), narrow(vp)
    x_orig = [xp, xw, xw]
    y_orig = [yp, yw, yw]
    u_orig = [up, uw, uw]
    v_orig = [vp, vw, vw]

    step_current = 0
    step_max = 50
    eps_pos = 1.0e-8
    eps_vel = 1.0e-5

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        # Predictor step
        xo, yo, uo, vo = _solve_free_pressure(xp, yp, up, vp, x_orig, y_orig, u_orig, v_orig, geom, params)
        # The first estimation can lead to complex numbers, keep only the real part
        xo = real_part(xo)
        yo = real_part(yo)
        uo = real_part(uo)
        vo = real_part(vo)
        while True:
            step_current += 1
            # Corrector step
            xcp = xp
            ycp = 0.5 * (yp + yo)
            ucp = 0.5 * (up + uo)
            vcp = 0.5 * (vp + vo)
            xn, yn, un, vn = _solve_free_pressure(
                xcp, ycp, ucp, vcp, x_orig, y_orig, u_orig, v_orig, geom, params
            )
            error_pos = mmax(xo - xn, yo - yn)
            error_vel = mmax(uo - un, vo - vn)
            xo = xn
            yo = yn
            uo = un
            vo = vn
            x_orig[2] = xo  # data for point 4, reinjected in the solver
            y_orig[2] = yo
            u_orig[2] = uo
            v_orig[2] = vo

            if abs(error_pos) < eps_pos and abs(error_vel) < eps_vel:
                break
            if step_current > step_max:
                raise MocConvergenceError(
                    "The maximum of iterations for the predictor-corrector algorithm has been "
                    "reached. Stopping the execution..."
                )
    return xo, yo, uo, vo


def _solve_free_pressure(xp, yp, up, vp, x_orig, y_orig, u_orig, v_orig, geom, params):
    # Velocity at solution point, flow is isentropic (free pressure)
    V4 = msqrt(
        2 * params.gamma * params.R * params.T
        * (1 - (params.Pstatic / params.P) ** ((params.gamma - 1) / params.gamma))
        / (params.gamma - 1)
    )
    Vmag = msqrt(up**2 + vp**2)
    a = get_speed_sound(params, Vmag)
    theta = atand(vp / up)
    alpha = asind(a / Vmag)
    lam = tand(theta + alpha)  # lambda+
    Q = up**2 - a**2
    R = 2 * up * vp - Q * lam
    S = geom.delta * a**2 * vp / yp
    lambda0 = (v_orig[1] + v_orig[2]) / (u_orig[1] + u_orig[2])  # slope of the jet boundary

    # C+ :  y4 - lambda  * x4 = y_orig(1) - lambda  * x_orig(1)
    # Jet:  y4 - lambda0 * x4 = y_orig(2) - lambda0 * x_orig(2)
    xn, yn = mldivide2(
        -lam, 1.0, -lambda0, 1.0,
        y_orig[0] - lam * x_orig[0],
        y_orig[1] - lambda0 * x_orig[1],
    )

    T = S * (xn - x_orig[0]) + Q * u_orig[0] + R * v_orig[0]
    un = (Q * T - R * msqrt(V4**2 * (Q**2 + R**2) - T**2)) / (Q**2 + R**2)
    vn = msqrt(V4**2 - un**2)
    return xn, yn, un, vn
