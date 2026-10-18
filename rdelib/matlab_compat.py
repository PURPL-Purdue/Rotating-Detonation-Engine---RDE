"""Small MATLAB-semantics helpers needed for faithful ports.

* :func:`colon`   - MATLAB ``a:step:b`` (inclusive end point, MATLAB's
  tolerance on the element count and symmetric element generation)
* :func:`num2str` - MATLAB ``num2str``/``string`` formatting of a scalar double
* :func:`fzero`   - MATLAB ``fzero`` with a scalar start point: same
  bracket-search sequence as MATLAB, then Brent's method on the bracket
"""

from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq

__all__ = ["colon", "num2str", "fzero"]

_EPS = np.finfo(float).eps


def colon(a: float, step: float, b: float) -> np.ndarray:
    """MATLAB ``a:step:b`` as a 1-D float array (inclusive end point).

    Uses MATLAB's algorithm: the number of intervals is rounded and accepted
    if ``a + n*step`` does not overshoot ``b`` by more than a few ulps, and the
    elements are generated symmetrically from both ends.
    """
    if step == 0 or (b - a) * step < 0:
        return np.zeros(0)
    tol = 2.0 * _EPS * max(abs(a), abs(b))
    sig = math.copysign(1.0, step)
    n = round((b - a) / step)
    if sig * (a + n * step - b) > tol:
        n -= 1
    last = a + n * step
    if sig * (last - b) > -tol and sig * (last - b) > 0:
        last = b
    out = np.empty(n + 1)
    k = np.arange(0, n // 2 + 1)
    out[k] = a + k * step
    out[n - k] = last - k * step
    if n % 2 == 0:
        out[n // 2] = (a + last) / 2.0
    return out


def num2str(x: float) -> str:
    """MATLAB ``num2str(x)`` for a real scalar (also what ``"txt" + x`` prints).

    Integers print without decimals; otherwise ``max(floor(log10|x|)+5, 5)``
    significant digits in ``%g`` style.
    """
    x = float(x)
    if math.isnan(x):
        return "NaN"
    if math.isinf(x):
        return "Inf" if x > 0 else "-Inf"
    if x == math.floor(x):
        return f"{x:.0f}"
    ndgt = math.floor(math.log10(abs(x))) if x != 0 else 0
    ndgt = max(ndgt + 5, 5)
    return f"{x:.{ndgt}g}"


def fzero(fun, x0, *, xtol: float | None = None):
    """MATLAB ``fzero`` for a scalar initial guess or a 2-element bracket.

    For a scalar ``x0`` the bracket is searched exactly as MATLAB does
    (``dx = x0/50`` - or 1/50 if ``x0 == 0`` - grown by sqrt(2) each try,
    testing ``x0-dx`` then ``x0+dx``), which determines *which* root is found
    when several exist. The root is then polished with Brent's method to
    MATLAB's default tolerance (TolX = eps).
    """
    x0 = np.atleast_1d(np.asarray(x0, dtype=float))
    if x0.size == 2:
        a, b = float(x0[0]), float(x0[1])
        fa, fb = fun(a), fun(b)
        if fa == 0:
            return a
        if fb == 0:
            return b
        if (fa > 0) == (fb > 0):
            raise ValueError("fzero: the function values at the interval endpoints must differ in sign.")
    else:
        x = float(x0[0])
        fx = fun(x)
        if fx == 0:
            return x
        dx = x / 50.0 if x != 0 else 1.0 / 50.0
        twosqrt = math.sqrt(2.0)
        a, fa, b, fb = x, fx, x, fx
        while (fa > 0) == (fb > 0):
            dx = twosqrt * dx
            a = x - dx
            fa = fun(a)
            if not (np.isfinite(fa) and np.isfinite(a)):
                raise ValueError("fzero: bracket search reached a non-finite value.")
            if (fa > 0) != (fb > 0):
                break
            b = x + dx
            fb = fun(b)
            if not (np.isfinite(fb) and np.isfinite(b)):
                raise ValueError("fzero: bracket search reached a non-finite value.")
        if fa == 0:
            return a
        if fb == 0:
            return b
    lo, hi = min(a, b), max(a, b)
    tol = xtol if xtol is not None else 2.0 * _EPS * max(abs(lo), abs(hi), 1.0)
    return brentq(fun, lo, hi, xtol=tol, rtol=4 * _EPS, maxiter=500)
