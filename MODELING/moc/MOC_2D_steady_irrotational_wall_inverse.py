"""Inverse wall point: flow data at a prescribed wall point.

Python translation of ``MOC_2D_steady_irrotational_wall_inverse.m`` from
xavierdechamps/Nozzle_MOC_supersonic (GPL-3.0); see MODELING/moc/README.md.
"""

import numpy as np

from ._errors import MocConvergenceError
from ._matlab import narrow, asind, atand, mldivide2, mmax, msqrt, tand
from .get_speed_sound import get_speed_sound


def MOC_2D_steady_irrotational_wall_inverse(xin, yin, uin, vin, theta4, geom, params):
    """Inverse method for the (u, v) data at a known wall point.

    ``xin..vin`` are 4-vectors: node 1 internal (known), node 2 internal
    (unknown, placeholder), node 3 upstream wall point (known), node 4 imposed
    wall point (known position, unknown data); ``theta4`` [deg] is the wall angle
    at node 4.  Returns ``(u4, v4, x2, y2, u2, v2)``: data at the imposed wall
    point and the position/data of node 2, the intersection of the C+
    characteristic through node 4 with the C- line 1-3.
    """
    step_current = 0
    step_max = 10
    eps_vel = 1.0e-5

    x = [narrow(val) for val in np.ravel(xin)]
    y = [narrow(val) for val in np.ravel(yin)]
    u = [narrow(val) for val in np.ravel(uin)]
    v = [narrow(val) for val in np.ravel(vin)]

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        # Predictor step
        x[1] = x[2]
        y[1] = y[2]
        u[1] = u[2]
        v[1] = v[2]
        x, y, u, v = _solve_wall_inverse(x, y, u, v, theta4, geom, params)
        while True:
            step_current += 1
            # Corrector step
            x[1] = x[2]
            y[1] = y[2]
            u[1] = 0.5 * (u[1] + u[3])
            v[1] = 0.5 * (v[1] + v[3])
            xn, yn, un, vn = _solve_wall_inverse(x, y, u, v, theta4, geom, params)

            error_vel = mmax(u[3] - un[3], v[3] - vn[3])
            x = xn
            y = yn
            u = un
            v = vn
            if abs(error_vel) < eps_vel:
                break
            if step_current > step_max:
                raise MocConvergenceError(
                    "The maximum of iterations for the predictor-corrector algorithm has been "
                    "reached. Stopping the execution..."
                )

    u4 = u[3]
    v4 = v[3]  # Data at the imposed wall point
    x2 = x[1]
    y2 = y[1]
    u2 = u[1]
    v2 = v[1]
    return u4, v4, x2, y2, u2, v2


def _solve_wall_inverse(x, y, u, v, theta4, geom, params):
    # Nodes (MATLAB 1..4 -> Python 0..3):
    #   1 internal (known position and data)  2 internal (unknown position and data)
    #   3 wall (known position and data)      4 wall (known position, unknown data)
    xn = list(x)
    yn = list(y)
    un = list(u)
    vn = list(v)
    # First step: position of the unknown internal point 2
    eps_vel = 1.0e-5
    steps_max = 10
    step_current = 1
    while True:
        Vmag = msqrt(un[1] ** 2 + vn[1] ** 2)
        a = get_speed_sound(params, Vmag)
        theta = atand(vn[1] / un[1])
        alpha = asind(a / Vmag)
        lam = tand(theta + alpha)  # lambda+
        # Slope of line 1-3
        slope13 = (y[2] - y[0]) / (x[2] - x[0])
        xn[1], yn[1] = mldivide2(
            -slope13, 1.0, -lam, 1.0,
            y[0] - slope13 * x[0],
            y[3] - lam * x[3],
        )
        # Interpolation of data at points 1 and 3 to get data at point 2
        unn = u[0] + (xn[1] - x[0]) / (x[2] - x[0]) * (u[2] - u[0])
        vnn = v[0] + (xn[1] - x[0]) / (x[2] - x[0]) * (v[2] - v[0])
        if abs(unn - un[1]) < eps_vel and abs(vnn - vn[1]) < eps_vel:
            break
        if step_current > steps_max:
            raise MocConvergenceError("Could not converge to a position for the new internal point")
        un[1] = unn
        vn[1] = vnn
        step_current += 1

    # Data at point 4 on the wall
    Q = un[1] ** 2 - a**2
    R = 2 * un[1] * vn[1] - Q * lam
    S = geom.delta * a**2 * vn[1] / yn[1]
    T = S * (x[3] - xn[1]) + Q * un[1] + R * vn[1]
    un[3], vn[3] = mldivide2(tand(theta4), -1.0, Q, R, 0.0, T)
    return xn, yn, un, vn
