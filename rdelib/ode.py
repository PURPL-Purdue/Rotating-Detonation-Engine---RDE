"""ODE / PDE integrators that reproduce MATLAB's ``ode45`` and ``pdepe``.

``ode45``
    A line-by-line reimplementation of MATLAB's ``ode45`` algorithm
    (Dormand-Prince 5(4) pair, MATLAB's initial-step heuristic, mixed
    relative/absolute *infinity-norm* error test, MATLAB step-size update
    and ``Refine = 4`` continuous-extension output). It exists for MATLAB
    *sample-grid compatibility*: use it when a port depends on the output
    time grid (e.g. "first sample below a target", which is quantised by the
    solver's steps). For new calculations prefer
    ``scipy.integrate.solve_ivp`` (with ``events`` for threshold crossings).

``pdepe``
    The Skeel-Berzins (1990) spatial discretisation used by MATLAB ``pdepe``
    for ``c(x,t,u,ux) u_t = x^-m (x^m f(x,t,u,ux))_x + s(x,t,u,ux)``
    with boundary conditions ``p + q f = 0``, integrated in time with a stiff
    BDF solver. Scalar PDEs (one component), ``m = 0`` or ``m = 1`` on a mesh
    with ``x[0] > 0`` when ``m = 1``, and ``q != 0`` at both boundaries
    (flux/Robin conditions) are supported, which covers the HADES use.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

__all__ = ["ode45", "ODE45Result", "pdepe"]


# ---------------------------------------------------------------------------
# ode45 (MATLAB algorithm)
# ---------------------------------------------------------------------------
_A = np.array([1 / 5, 3 / 10, 4 / 5, 8 / 9, 1.0, 1.0])
_B = np.array(
    [
        [1 / 5, 3 / 40, 44 / 45, 19372 / 6561, 9017 / 3168, 35 / 384],
        [0, 9 / 40, -56 / 15, -25360 / 2187, -355 / 33, 0],
        [0, 0, 32 / 9, 64448 / 6561, 46732 / 5247, 500 / 1113],
        [0, 0, 0, -212 / 729, 49 / 176, 125 / 192],
        [0, 0, 0, 0, -5103 / 18656, -2187 / 6784],
        [0, 0, 0, 0, 0, 11 / 84],
        [0, 0, 0, 0, 0, 0],
    ]
)
_E = np.array([71 / 57600, 0, -71 / 16695, 71 / 1920, -17253 / 339200, 22 / 525, -1 / 40])
_BI = np.array(
    [
        [1, -183 / 64, 37 / 12, -145 / 128],
        [0, 0, 0, 0],
        [0, 1500 / 371, -1000 / 159, 1000 / 371],
        [0, -125 / 32, 125 / 12, -375 / 64],
        [0, 9477 / 3392, -729 / 106, 25515 / 6784],
        [0, -11 / 7, 11 / 3, -55 / 28],
        [0, 3 / 2, -4, 5 / 2],
    ]
)


@dataclass
class ODE45Result:
    """Output of :func:`ode45`: ``t`` (n,), ``y`` (n, neq) like MATLAB ``[t, y]``."""

    t: np.ndarray
    y: np.ndarray
    nsteps: int
    nfailed: int
    nfevals: int


def _finite_state(v, t):
    """Reject overflowed / NaN stage, accepted or interpolated states."""
    if not np.all(np.isfinite(v)):
        raise FloatingPointError(f"ode45: non-finite solution state in the step from t = {t!r}")
    return v


def _ntrp45(tinterp, t, y, h, f):
    s = (np.atleast_1d(tinterp) - t) / h
    S = np.vstack([s, s**2, s**3, s**4])
    return y[:, None] + f @ (h * _BI) @ S


def ode45(fun, tspan, y0, *, rtol: float = 1e-3, atol: float = 1e-6, refine: int = 4,
          max_step: float | None = None, initial_step: float | None = None) -> ODE45Result:
    """MATLAB ``[t, y] = ode45(fun, tspan, y0)`` with default options.

    ``fun(t, y)`` returns dy/dt (1-D array). ``tspan`` is ``[t0, tfinal]``
    (output at refined solver steps) or a vector of >2 output times
    (output only at those times, like MATLAB). ``rtol``/``atol`` default to
    MATLAB's 1e-3 / 1e-6; ``atol`` must be scalar.
    """
    tspan = np.asarray(tspan, dtype=float).ravel()
    if tspan.size < 2:
        raise ValueError("tspan must have at least two entries")
    if not np.all(np.isfinite(tspan)):
        raise ValueError("tspan must be finite")
    t0, tfinal = float(tspan[0]), float(tspan[-1])
    if t0 == tfinal:
        raise ValueError("tspan start and end must differ")
    tdir = np.sign(tfinal - t0)
    if tspan.size > 2 and not np.all(tdir * np.diff(tspan) > 0):
        raise ValueError("tspan must be strictly monotonic")
    y = np.asarray(y0, dtype=float).ravel().copy()
    if y.size == 0 or not np.all(np.isfinite(y)):
        raise ValueError("y0 must be a non-empty finite vector")
    neq = y.size
    for name, v in (("rtol", rtol), ("atol", atol)):
        if not (np.isscalar(v) and np.isfinite(v) and v > 0):
            raise ValueError(f"{name} must be a positive finite scalar")
    if max_step is not None and not (np.isfinite(max_step) and max_step > 0):
        raise ValueError("max_step must be positive and finite")
    if initial_step is not None and not (np.isfinite(initial_step) and initial_step > 0):
        raise ValueError("initial_step must be positive and finite")
    if int(refine) != refine or refine < 1:
        raise ValueError("refine must be a positive integer")
    refine = int(refine)
    rtol = max(rtol, 100 * np.finfo(float).eps)
    threshold = atol / rtol
    hmax = abs(max_step) if max_step is not None else 0.1 * abs(tfinal - t0)
    htspan = abs(tspan[1] - tspan[0])
    ntspan = tspan.size
    pow_ = 1.0 / 5.0

    def F(t, yy):
        d = np.asarray(fun(t, yy), dtype=float).ravel()
        if d.shape != (neq,):
            raise ValueError(f"ode45: fun returned {d.size} values for {neq} equations")
        if not np.all(np.isfinite(d)):
            raise FloatingPointError(f"ode45: non-finite derivative at t = {t!r}")
        return d

    def eps_t(t):  # MATLAB eps(t): spacing of |t| (always positive)
        return np.spacing(abs(t))

    f0 = F(t0, y)
    nfevals = 1
    f = np.zeros((neq, 7))
    t = t0
    hmin = 16 * eps_t(t)
    if initial_step is None:
        absh = min(hmax, htspan)
        rh = np.max(np.abs(f0 / np.maximum(np.abs(y), threshold))) / (0.8 * rtol**pow_)
        if absh * rh > 1:
            absh = 1.0 / rh
        absh = max(absh, hmin)
    else:
        absh = min(hmax, max(hmin, initial_step))
    f[:, 0] = f0

    if ntspan > 2:
        output_at = "requested"
    elif refine <= 1:
        output_at = "steps"
    else:
        output_at = "refined"
        S = np.arange(1, refine) / refine

    tout = [t0]
    yout = [y.copy()]
    next_idx = 1  # next requested tspan index
    nsteps = nfailed = 0
    done = False
    while not done:
        hmin = 16 * eps_t(t)
        absh = min(hmax, max(hmin, absh))
        h = tdir * absh
        if 1.1 * absh >= abs(tfinal - t):
            h = tfinal - t
            absh = abs(h)
            done = True
        nofailed = True
        while True:
            hA = h * _A
            hB = h * _B
            f[:, 1] = F(t + hA[0], _finite_state(y + f @ hB[:, 0], t))
            f[:, 2] = F(t + hA[1], _finite_state(y + f @ hB[:, 1], t))
            f[:, 3] = F(t + hA[2], _finite_state(y + f @ hB[:, 2], t))
            f[:, 4] = F(t + hA[3], _finite_state(y + f @ hB[:, 3], t))
            f[:, 5] = F(t + hA[4], _finite_state(y + f @ hB[:, 4], t))
            tnew = t + hA[5]
            if done:
                tnew = tfinal
            h = tnew - t
            if h == 0:
                raise RuntimeError(f"ode45: step size underflow at t = {t!r}")
            ynew = _finite_state(y + f @ hB[:, 5], t)
            f[:, 6] = F(tnew, ynew)
            nfevals += 6
            err = absh * np.max(np.abs((f @ _E) / np.maximum(np.maximum(np.abs(y), np.abs(ynew)), threshold)))
            if not np.isfinite(err):
                raise FloatingPointError(f"ode45: non-finite error estimate at t = {t!r}")
            if err > rtol:
                nfailed += 1
                if absh <= hmin:
                    raise RuntimeError(f"ode45: unable to meet integration tolerances at t={t}")
                if nofailed:
                    nofailed = False
                    absh = max(hmin, absh * max(0.1, 0.8 * (rtol / err) ** pow_))
                else:
                    absh = max(hmin, 0.5 * absh)
                h = tdir * absh
                done = False
            else:
                break
        nsteps += 1

        if output_at == "refined":
            tref = t + (tnew - t) * S
            yref = _finite_state(_ntrp45(tref, t, y, h, f), t)
            tout.extend(tref.tolist())
            yout.extend(list(yref.T))
            tout.append(tnew)
            yout.append(ynew.copy())
        elif output_at == "steps":
            tout.append(tnew)
            yout.append(ynew.copy())
        else:
            while next_idx < ntspan and tdir * (tnew - tspan[next_idx]) >= 0:
                tn = tspan[next_idx]
                if tn == tnew:
                    yn = ynew.copy()
                else:
                    yn = _finite_state(_ntrp45(np.array([tn]), t, y, h, f)[:, 0], t)
                tout.append(float(tn))
                yout.append(yn)
                next_idx += 1

        if done:
            break
        if nofailed:
            temp = 1.25 * (err / rtol) ** pow_
            if temp > 0.2:
                absh = absh / temp
            else:
                absh = 5.0 * absh
        t = tnew
        y = ynew
        f[:, 0] = f[:, 6]

    return ODE45Result(t=np.asarray(tout), y=np.vstack(yout), nsteps=nsteps, nfailed=nfailed, nfevals=nfevals)


# ---------------------------------------------------------------------------
# pdepe (Skeel-Berzins discretisation)
# ---------------------------------------------------------------------------
def pdepe(m, pdefun, icfun, bcfun, xmesh, tspan, *, rtol: float = 1e-6, atol: float = 1e-8,
          method: str = "BDF"):
    """Solve a scalar 1-D parabolic PDE like MATLAB ``sol = pdepe(...)``.

    Signatures (same as MATLAB, scalar ``u``)::

        c, f, s = pdefun(x, t, u, dudx)
        u0      = icfun(x)
        pl, ql, pr, qr = bcfun(xl, ul, xr, ur, t)

    Returns ``sol`` with shape ``(len(tspan), len(xmesh))`` (MATLAB
    ``sol(:,:,1)``).

    ``rtol``/``atol`` default to 1e-6/1e-8, tighter than MATLAB's pdepe
    defaults (1e-3/1e-6): scipy's BDF at 1e-3 is measurably less accurate
    than MATLAB/Octave ``ode15s`` at 1e-3 (time-to-limit in
    HADES_size_thermal_v1 off by 0.07 % vs 0.02 % for Octave ode15s), so the
    tighter default lands closer to what MATLAB actually returns.
    """
    x = np.asarray(xmesh, dtype=float).ravel()
    tspan = np.asarray(tspan, dtype=float).ravel()
    nx = x.size
    if nx < 3:
        raise ValueError("xmesh must have at least 3 points")
    xL, xR = x[:-1], x[1:]
    if m == 0:
        xi = 0.5 * (xL + xR)
        zeta = xi.copy()
    elif m == 1:
        if x[0] <= 0:
            raise NotImplementedError("pdepe: m=1 with x[0] == 0 (singular) is not implemented")
        xi = (xR - xL) / np.log(xR / xL)
        zeta = np.sqrt(xi * (xR + xL) / 2.0)
    else:
        raise NotImplementedError("pdepe: only m = 0 and m = 1 are implemented")
    mp1 = m + 1
    wR = (zeta**mp1 - xL**mp1) / mp1  # weight of interval i seen from its left node i
    wL = (xR**mp1 - zeta**mp1) / mp1  # weight of interval i seen from its right node i+1
    xim = xi**m

    def interp(uL, uR):
        du = uR - uL
        if m == 0:
            return uL + du * (xi - xL) / (xR - xL), du / (xR - xL)
        lg = np.log(xR / xL)
        return uL + du * np.log(xi / xL) / lg, du / (xi * lg)

    def rhs(t, u):
        U, Ux = interp(u[:-1], u[1:])
        c = np.empty(nx - 1)
        fl = np.empty(nx - 1)
        s = np.empty(nx - 1)
        for i in range(nx - 1):
            ci, fi, si = pdefun(xi[i], t, U[i], Ux[i])
            c[i], fl[i], s[i] = ci, fi, si
        pl, ql, pr, qr = bcfun(x[0], u[0], x[-1], u[-1], t)
        if ql == 0 or qr == 0:
            raise NotImplementedError("pdepe: Dirichlet (q = 0) boundaries are not implemented")
        flux = xim * fl
        dudt = np.empty(nx)
        dudt[0] = (flux[0] + x[0] ** m * pl / ql + wR[0] * s[0]) / (wR[0] * c[0])
        dudt[1:-1] = (flux[1:] - flux[:-1] + wL[:-1] * s[:-1] + wR[1:] * s[1:]) / (
            wL[:-1] * c[:-1] + wR[1:] * c[1:]
        )
        dudt[-1] = (-x[-1] ** m * pr / qr - flux[-1] + wL[-1] * s[-1]) / (wL[-1] * c[-1])
        return dudt

    u0 = np.array([float(icfun(xx)) for xx in x])
    sol = solve_ivp(rhs, (tspan[0], tspan[-1]), u0, method=method, t_eval=tspan, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(f"pdepe: time integration failed: {sol.message}")
    return sol.y.T
