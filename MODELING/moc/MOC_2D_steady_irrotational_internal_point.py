"""Interior point: intersection of a C+ and a C- characteristic.

Python translation of ``MOC_2D_steady_irrotational_internal_point.m`` from
xavierdechamps/Nozzle_MOC_supersonic (GPL-3.0); see MODELING/moc/README.md.
Modified Euler predictor-corrector of Zucrow & Hoffman, Gas Dynamics Vol. II.
"""

import numpy as np

from ._errors import MocConvergenceError
from ._matlab import narrow, asind, atand, mldivide2, mmax, msqrt, tand
from .get_speed_sound import get_speed_sound


def MOC_2D_steady_irrotational_internal_point(xp, yp, up, vp, xm, ym, um, vm, geom, params):
    """Return ``(xo, yo, uo, vo)`` at the intersection of two characteristics.

    ``xp, yp, up, vp``: conditions on the left-running C+ characteristic;
    ``xm, ym, um, vm``: conditions on the right-running C- characteristic.
    Raises :class:`MocConvergenceError` after 50 corrector steps.
    """
    xp, yp, up, vp = narrow(xp), narrow(yp), narrow(up), narrow(vp)
    xm, ym, um, vm = narrow(xm), narrow(ym), narrow(um), narrow(vm)
    x_orig = (xp, xm)
    y_orig = (yp, ym)
    u_orig = (up, um)
    v_orig = (vp, vm)

    step_current = 0
    step_max = 50
    eps_pos = 1.0e-8
    eps_vel = 1.0e-5

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        # Predictor step: first estimation of the intersection of the 2 characteristics
        xo, yo, uo, vo = _solve_internal_point(
            xp, yp, up, vp, xm, ym, um, vm, x_orig, y_orig, u_orig, v_orig, geom, params
        )
        while True:
            step_current += 1
            # Corrector step
            xcp = xp
            ycp = 0.5 * (yp + yo)
            ucp = 0.5 * (up + uo)
            vcp = 0.5 * (vp + vo)
            xcm = xm
            ycm = 0.5 * (ym + yo)
            ucm = 0.5 * (um + uo)
            vcm = 0.5 * (vm + vo)
            xn, yn, un, vn = _solve_internal_point(
                xcp, ycp, ucp, vcp, xcm, ycm, ucm, vcm, x_orig, y_orig, u_orig, v_orig, geom, params
            )
            error_pos = mmax(xo - xn, yo - yn)
            error_vel = mmax(uo - un, vo - vn)
            xo = xn
            yo = yn
            uo = un
            vo = vn
            # Check if we converged on the position and on the velocity components
            if abs(error_pos) < eps_pos and abs(error_vel) < eps_vel:
                break
            if step_current > step_max:
                raise MocConvergenceError(
                    "The maximum of iterations for the predictor-corrector algorithm has been "
                    "reached. Stopping the execution..."
                )
    return xo, yo, uo, vo


def _solve_internal_point(xp, yp, up, vp, xm, ym, um, vm, x_orig, y_orig, u_orig, v_orig, geom, params):
    # Index 1 (MATLAB) -> C+ data (p), index 2 -> C- data (m)
    V1 = msqrt(up**2 + vp**2)
    V2 = msqrt(um**2 + vm**2)
    a1 = get_speed_sound(params, V1)
    a2 = get_speed_sound(params, V2)
    theta1 = atand(vp / up)
    theta2 = atand(vm / um)
    alpha1 = asind(a1 / V1)
    alpha2 = asind(a2 / V2)
    lambda1 = tand(theta1 + alpha1)  # lambda+
    lambda2 = tand(theta2 - alpha2)  # lambda-
    Q1 = up**2 - a1**2
    Q2 = um**2 - a2**2
    R1 = 2 * up * vp - Q1 * lambda1
    R2 = 2 * um * vm - Q2 * lambda2
    S1 = geom.delta * a1**2 * vp / yp
    S2 = geom.delta * a2**2 * vm / ym

    # Solve the system to get the position of the intersection of the C+ and C- characteristics
    xn, yn = mldivide2(
        -lambda1, 1.0, -lambda2, 1.0,
        y_orig[0] - lambda1 * x_orig[0],
        y_orig[1] - lambda2 * x_orig[1],
    )

    T1 = S1 * (xn - x_orig[0]) + Q1 * u_orig[0] + R1 * v_orig[0]
    T2 = S2 * (xn - x_orig[1]) + Q2 * u_orig[1] + R2 * v_orig[1]
    un, vn = mldivide2(Q1, R1, Q2, R2, T1, T2)
    return xn, yn, un, vn
