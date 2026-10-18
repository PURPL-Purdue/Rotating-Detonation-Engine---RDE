"""Small MATLAB-semantics helpers used by the faithful MoC port.

These exist so the translated solver behaves like the MATLAB original in the
places where Python/NumPy semantics differ:

* ``colon``          -- MATLAB's ``a:d:b`` (symmetric construction, inclusive end)
* ``sind``/``cosd``/``tand``/``asind``/``atand`` -- degree trig, exact zeros at
  multiples of 180 deg, complex results for ``asind`` outside [-1, 1]
* ``msqrt``          -- ``sqrt`` that returns a complex value for negative input
* ``mmax``           -- ``max([...])`` that ignores NaN and orders complex by abs
* ``mldivide2``      -- ``A\\b`` for 2x2 systems (LAPACK getrf/getrs order of
  operations; singular systems give Inf/NaN like MATLAB instead of raising)
* ``fzero``          -- MATLAB's bracketed ``fzero`` (Brent/Dekker zeroin)
* ``GrowArray``      -- auto-growing, zero-padded 2-D array addressed with
  1-based (row, column) indices, mirroring MATLAB's ``X(J,I) = value``
* ``num2str``        -- MATLAB ``num2str`` for real scalars
"""

from __future__ import annotations

import cmath
import math

import numpy as np

EPS = float(np.finfo(float).eps)
_D2R = math.pi / 180.0
_R2D = 180.0 / math.pi


# --------------------------------------------------------------------------
# Ranges
# --------------------------------------------------------------------------
def colon(a: float, d: float, b: float | None = None) -> np.ndarray:
    """MATLAB ``a:d:b`` (or ``a:b`` when called as ``colon(a, b)``).

    Follows the published MATLAB algorithm (C. Moler, "colonop"): the vector is
    built symmetrically from both ends so the last element equals ``b`` when the
    range lands on it.
    """
    if b is None:
        a, d, b = a, 1.0, d
    a = float(a)
    d = float(d)
    b = float(b)
    tol = 2.0 * EPS * max(abs(a), abs(b))
    sig = math.copysign(1.0, d) if d != 0 else 0.0
    if not (math.isfinite(a) and math.isfinite(d) and math.isfinite(b)):
        return np.array([np.nan])
    if d == 0 or (a < b and d < 0) or (b < a and d > 0):
        return np.zeros(0)
    if a == math.floor(a) and d == 1:
        n = math.floor(b) - a
    elif a == math.floor(a) and d == math.floor(d):
        q = math.floor(a / d)
        r = a - q * d
        n = math.floor((b - r) / d) - q
    else:
        n = round((b - a) / d)
        if sig * (a + n * d - b) > tol:
            n = n - 1
    n = int(n)
    c = a + n * d
    if sig * (c - b) > -tol:
        c = b
    v = np.zeros(n + 1)
    k = np.arange(0, n // 2 + 1)
    v[k] = a + k * d
    v[n - k] = c - k * d
    if n % 2 == 0:
        v[n // 2] = (a + c) / 2.0
    return v


# --------------------------------------------------------------------------
# Degree trigonometry
# --------------------------------------------------------------------------
def _is_complex(x) -> bool:
    return isinstance(x, complex) or np.iscomplexobj(x)


def narrow(x):
    """Scalar as MATLAB holds it: real ``float64`` unless the imaginary part is nonzero.

    MATLAB/Octave drop an all-zero imaginary part (e.g. when indexing a complex
    array), which matters for relational operators, ``max`` and branch cuts.
    """
    if isinstance(x, (complex, np.complexfloating)):
        if x.imag == 0:
            return np.float64(x.real)
        return np.complex128(x)
    return np.float64(x)


def _complex_trig_arg(xa):
    """Split a complex degree argument: real part reduced exactly mod 360 (sign kept)."""
    r = np.fmod(np.real(xa), 360.0)
    return r, (r + 1j * np.imag(xa)) * _D2R


def sind(x):
    """MATLAB ``sind`` (scalar or array).

    Real arguments: exactly 0 at integer multiples of 180.  Complex arguments:
    ``sin(z*pi/180)`` with the real part reduced exactly mod 360; only the real
    part is set to exactly 0 at multiples of 180 (where it vanishes
    mathematically), as Octave does.
    """
    if np.ndim(x) == 0:
        x = narrow(x)
        if _is_complex(x):
            return np.complex128(sind(np.asarray([x]))[0])
        x = float(x)
        if not math.isfinite(x):
            return np.float64(np.nan)
        r = math.fmod(x, 360.0)  # exact
        if r == 0.0 or abs(r) == 180.0:
            return np.float64(0.0)
        return np.float64(math.sin(r * _D2R))
    xa = np.asarray(x)
    if _is_complex(xa):
        r, z = _complex_trig_arg(xa)
        y = np.sin(z)
        on_axis = (r == 0.0) | (np.abs(r) == 180.0)
        return np.where(on_axis, 1j * y.imag, y)
    with np.errstate(invalid="ignore"):
        r = np.fmod(xa.astype(float), 360.0)
        y = np.sin(r * _D2R)
    return np.where((r == 0.0) | (np.abs(r) == 180.0), 0.0, y)


def cosd(x):
    """MATLAB ``cosd`` (scalar or array).

    Real arguments: exactly 0 at odd multiples of 90.  Complex arguments:
    ``cos(z*pi/180)`` with the sign of the real part preserved; only the real
    part is set to exactly 0 at odd multiples of 90.
    """
    if np.ndim(x) == 0:
        x = narrow(x)
        if _is_complex(x):
            return np.complex128(cosd(np.asarray([x]))[0])
        x = float(x)
        if not math.isfinite(x):
            return np.float64(np.nan)
        r = math.fmod(abs(x), 360.0)
        if r == 90.0 or r == 270.0:
            return np.float64(0.0)
        return np.float64(math.cos(r * _D2R))
    xa = np.asarray(x)
    if _is_complex(xa):
        r, z = _complex_trig_arg(xa)
        y = np.cos(z)
        odd90 = (np.abs(r) == 90.0) | (np.abs(r) == 270.0)
        return np.where(odd90, 1j * y.imag, y)
    with np.errstate(invalid="ignore"):
        r = np.fmod(np.abs(xa.astype(float)), 360.0)
        y = np.cos(r * _D2R)
    return np.where((r == 90.0) | (r == 270.0), 0.0, y)


def tand(x):
    """MATLAB/Octave ``tand``: 0 at multiples of 180, Inf at odd multiples of 90 (real input)."""
    if np.ndim(x) == 0:
        x = narrow(x)
    if _is_complex(x):
        if np.ndim(x) == 0:
            return np.complex128(np.tan(np.complex128(x) / 180.0 * math.pi))
        return np.tan(np.asarray(x) / 180.0 * math.pi)
    if np.ndim(x) == 0:
        x = float(x)
        i0 = x / 180.0
        i90 = (x - 90.0) / 180.0
        if math.isfinite(i0) and i0 == math.trunc(i0):
            return np.float64(0.0)
        if math.isfinite(i90) and i90 == math.trunc(i90):
            return np.float64(np.inf)
        return np.float64(math.tan(i0 * math.pi))
    xa = np.asarray(x, dtype=float)
    i0 = xa / 180.0
    i90 = (xa - 90.0) / 180.0
    y = np.tan(i0 * math.pi)
    y = np.where(np.isfinite(i0) & (i0 == np.trunc(i0)), 0.0, y)
    y = np.where(np.isfinite(i90) & (i90 == np.trunc(i90)), np.inf, y)
    return y


def asind(x):
    """MATLAB ``asind`` (scalar or array).

    Real |x| > 1 gives a complex result on MATLAB's branch
    (``asind(2) = 90 - 75.4561i``, ``asind(-3) = -90 + 100.998i``) for scalars
    and arrays alike; an array containing such an element becomes complex.
    """
    if np.ndim(x) == 0:
        x = narrow(x)
        if _is_complex(x):
            return np.complex128(cmath.asin(complex(x)) * 180.0 / math.pi)
        x = float(x)
        if x > 1.0:
            return np.complex128(cmath.asin(complex(x, -0.0)) * 180.0 / math.pi)
        if x < -1.0:
            return np.complex128(cmath.asin(complex(x, 0.0)) * 180.0 / math.pi)
        if math.isnan(x):
            return np.float64(np.nan)
        return np.float64(math.asin(x) * 180.0 / math.pi)
    xa = np.asarray(x)
    re = np.real(xa)
    if not _is_complex(xa) and not np.any(np.abs(re) > 1.0):
        with np.errstate(invalid="ignore"):
            return np.arcsin(xa.astype(float)) * 180.0 / math.pi
    im = np.imag(xa).astype(float) if _is_complex(xa) else np.zeros(re.shape)
    # elements with a zero imaginary part follow the real-argument (MATLAB) branch
    im = np.where((im == 0) & (re > 1.0), -0.0, im)
    im = np.where((im == 0) & (re < -1.0), 0.0, im)
    z = re.astype(float) + 0j
    z.imag = im
    out = np.arcsin(z) * 180.0 / math.pi
    # real elements inside [-1, 1] keep the real-argument value (as for scalars)
    inside = (im == 0) & (np.abs(re) <= 1.0)
    out[inside] = np.arcsin(re[inside].astype(float)) * 180.0 / math.pi
    return out


def atand(x):
    """MATLAB ``atand`` (``180/pi*atan(x)``)."""
    if np.ndim(x) == 0:
        x = narrow(x)
    if np.ndim(x) == 0 and not _is_complex(x):
        return np.float64(_R2D * math.atan(float(x)))
    if np.ndim(x) == 0:
        return np.complex128(_R2D * np.arctan(np.complex128(x)))
    return _R2D * np.arctan(np.asarray(x))


def msqrt(x):
    """MATLAB ``sqrt``: complex principal root for negative (or complex) input."""
    if np.ndim(x) == 0:
        x = narrow(x)
        if _is_complex(x):
            return np.complex128(cmath.sqrt(complex(x)))
        x = float(x)
        if x < 0.0:
            return np.complex128(cmath.sqrt(complex(x, 0.0)))
        return np.float64(math.sqrt(x)) if not math.isnan(x) else np.float64(np.nan)
    return np.emath.sqrt(np.asarray(x))


def real_part(x):
    """Real part of a scalar (``real(x)``)."""
    return x.real if _is_complex(x) else x


# Relational operators on complex scalars.  MATLAB compares only the real
# parts; GNU Octave orders complex numbers by modulus, then argument.  The port
# follows MATLAB.  ``complex_ordering("octave")`` exists only so the regression
# tests can reproduce Octave reference runs node-by-node where the (unphysical)
# complex values of the original free-pressure point propagate.
_COMPLEX_ORDER = ["matlab"]


class complex_ordering:
    """Context manager selecting ``"matlab"`` (default) or ``"octave"`` complex ordering."""

    def __init__(self, mode: str):
        if mode not in ("matlab", "octave"):
            raise ValueError(mode)
        self.mode = mode

    def __enter__(self):
        self._old = _COMPLEX_ORDER[0]
        _COMPLEX_ORDER[0] = self.mode
        return self

    def __exit__(self, *exc):
        _COMPLEX_ORDER[0] = self._old
        return False


def _cmp_keys(a, b):
    a = narrow(a)
    b = narrow(b)
    if _COMPLEX_ORDER[0] == "octave" and (_is_complex(a) or _is_complex(b)):
        za, zb = complex(a), complex(b)
        return (abs(za), cmath.phase(za)), (abs(zb), cmath.phase(zb))
    return real_part(a), real_part(b)


def mlt(a, b) -> bool:
    """MATLAB ``a < b`` for scalars (complex: real parts)."""
    ka, kb = _cmp_keys(a, b)
    return bool(ka < kb)


def mle(a, b) -> bool:
    """MATLAB ``a <= b`` for scalars."""
    ka, kb = _cmp_keys(a, b)
    return bool(ka <= kb)


def mgt(a, b) -> bool:
    """MATLAB ``a > b`` for scalars."""
    ka, kb = _cmp_keys(a, b)
    return bool(ka > kb)


def mge(a, b) -> bool:
    """MATLAB ``a >= b`` for scalars."""
    ka, kb = _cmp_keys(a, b)
    return bool(ka >= kb)


def mmax(*vals):
    """MATLAB ``max([v1, v2, ...])`` for scalars: NaN ignored, complex by abs."""
    vals = [narrow(v) for v in vals]
    vals = [v for v in vals if not (np.isnan(v.real) or np.isnan(v.imag) if _is_complex(v) else np.isnan(v))]
    if not vals:
        return np.float64(np.nan)
    if any(_is_complex(v) for v in vals):
        # MATLAB orders complex numbers by abs, then by angle
        return max(vals, key=lambda v: (abs(v), np.angle(v)))
    best = vals[0]
    for v in vals[1:]:
        if v > best:
            best = v
    return best


# --------------------------------------------------------------------------
# Linear algebra
# --------------------------------------------------------------------------
def mldivide2(a11, a12, a21, a22, b1, b2):
    """Solve ``[a11 a12; a21 a22] \\ [b1; b2]`` like MATLAB/LAPACK (dgesv).

    LU with partial pivoting (first max |a_i1| pivot), multiplier formed with the
    reciprocal of the pivot as in LAPACK's ``dgetrf2``, then forward/back
    substitution as in ``dgetrs``. Singular systems yield Inf/NaN (MATLAB only
    warns) rather than raising.
    """
    cplx = any(_is_complex(z) for z in (a11, a12, a21, a22, b1, b2))
    if cplx:
        a11, a12, a21, a22, b1, b2 = (np.complex128(z) for z in (a11, a12, a21, a22, b1, b2))

        def mag(z):  # LAPACK izamax uses |Re| + |Im|
            return abs(z.real) + abs(z.imag)
    else:
        mag = abs
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        if mag(a21) > mag(a11):  # row swap
            a11, a12, a21, a22 = a21, a22, a11, a12
            b1, b2 = b2, b1
        if a11 != 0:
            if abs(a11) >= np.finfo(float).tiny:
                l21 = a21 * (1.0 / a11 if cplx else np.float64(1.0) / a11)
            else:
                l21 = a21 / a11
        else:
            l21 = a21  # LAPACK leaves the column unscaled (INFO > 0)
        u22 = a22 - l21 * a12
        y2 = b2 - l21 * b1
        x2 = np.float64(y2) / u22 if not _is_complex(y2) else y2 / u22
        x1 = (b1 - x2 * a12) / np.float64(a11) if not _is_complex(a11) else (b1 - x2 * a12) / a11
    return x1, x2


# --------------------------------------------------------------------------
# Root finding
# --------------------------------------------------------------------------
def fzero(fun, x0, tolx: float = EPS):
    """MATLAB ``fzero(fun, [a b])`` / ``fzero(fun, [a b], optimset('TolX',tolx))``.

    Port of the bracketed algorithm used by MATLAB's ``fzero`` (Dekker/Brent
    ``zeroin``: bisection, secant and inverse quadratic interpolation) with the
    MATLAB convergence test ``abs(m) <= 2*tol*max(abs(b),1)``.
    """
    a, b = float(x0[0]), float(x0[1])
    fa = float(fun(a))
    fb = float(fun(b))
    if not (math.isfinite(fa) and math.isfinite(fb)):
        raise ValueError("fzero: function values at interval endpoints must be finite and real.")
    if fa == 0:
        return a
    if fb == 0:
        return b
    if (fa > 0) == (fb > 0):
        raise ValueError("fzero: the function values at the interval endpoints must differ in sign.")
    tol = tolx
    fc = fb
    c = b
    d = e = b - a
    while fb != 0 and a != b:
        # Ensure b is the best result so far, a the previous value of b and
        # c on the opposite side of the zero from b.
        if (fb > 0) == (fc > 0):
            c = a
            fc = fa
            d = b - a
            e = d
        if abs(fc) < abs(fb):
            a = b
            b = c
            c = a
            fa = fb
            fb = fc
            fc = fa
        # Convergence test and possible exit
        m = 0.5 * (c - b)
        toler = 2.0 * tol * max(abs(b), 1.0)
        if abs(m) <= toler or fb == 0.0:
            break
        # Choose bisection or interpolation
        if abs(e) < toler or abs(fa) <= abs(fb):
            d = m
            e = m
        else:
            s = fb / fa
            if a == c:
                # Linear interpolation
                p = 2.0 * m * s
                q = 1.0 - s
            else:
                # Inverse quadratic interpolation
                q = fa / fc
                r = fb / fc
                p = s * (2.0 * m * q * (q - r) - (b - a) * (r - 1.0))
                q = (q - 1.0) * (r - 1.0) * (s - 1.0)
            if p > 0:
                q = -q
            else:
                p = -p
            # Is the interpolated point acceptable?
            if (2.0 * p < 3.0 * m * q - abs(toler * q)) and (p < abs(0.5 * e * q)):
                e = d
                d = p / q
            else:
                d = m
                e = m
        # Next point
        a = b
        fa = fb
        if abs(d) > toler:
            b = b + d
        elif b > c:
            b = b - toler
        else:
            b = b + toler
        fb = float(fun(b))
    return b


# --------------------------------------------------------------------------
# Growing arrays
# --------------------------------------------------------------------------
class GrowArray:
    """2-D array with MATLAB growth semantics and 1-based (row, col) indexing.

    ``A[J, I] = v`` grows the array (zero padded) when (J, I) is out of range,
    exactly like ``A(J,I) = v`` in MATLAB.  Reading outside the current size
    raises ``IndexError`` ("Index exceeds matrix dimensions").  If a complex
    value is stored, the whole array becomes complex (as in MATLAB).
    """

    __slots__ = ("_buf", "nrows", "ncols")

    def __init__(self, init=None):
        if init is None:
            arr = np.zeros((0, 0))
        else:
            arr = np.asarray(init)
            if arr.ndim == 0:
                arr = arr.reshape(1, 1)
            elif arr.ndim == 1:
                arr = arr.reshape(1, -1)
        dtype = np.complex128 if np.iscomplexobj(arr) else np.float64
        self.nrows, self.ncols = arr.shape
        self._buf = np.zeros((max(16, self.nrows), max(16, self.ncols)), dtype=dtype)
        self._buf[: self.nrows, : self.ncols] = arr

    @property
    def shape(self):
        return (self.nrows, self.ncols)

    @property
    def capacity(self):
        """Allocated (rows, cols) of the backing buffer."""
        return self._buf.shape

    def _ensure(self, r: int, c: int) -> None:
        # Amortised doubling per axis: only an axis whose capacity is exceeded grows.
        cap_r, cap_c = self._buf.shape
        new_r = cap_r if r <= cap_r else max(r, 2 * cap_r)
        new_c = cap_c if c <= cap_c else max(c, 2 * cap_c)
        if (new_r, new_c) != (cap_r, cap_c):
            new = np.zeros((new_r, new_c), dtype=self._buf.dtype)
            new[:cap_r, :cap_c] = self._buf
            self._buf = new
        self.nrows = max(self.nrows, r)
        self.ncols = max(self.ncols, c)

    def __getitem__(self, idx):
        r, c = idx
        if not (1 <= r <= self.nrows and 1 <= c <= self.ncols):
            raise IndexError(
                f"Index ({r},{c}) exceeds matrix dimensions ({self.nrows}x{self.ncols})"
            )
        val = self._buf[r - 1, c - 1]
        if val.dtype.kind == "c" and val.imag == 0:
            return np.float64(val.real)  # MATLAB narrows on indexing
        return val

    def __setitem__(self, idx, val):
        r, c = idx
        if r < 1 or c < 1:
            raise IndexError("MATLAB indices must be positive integers")
        if _is_complex(val) and not np.iscomplexobj(self._buf):
            self._buf = self._buf.astype(np.complex128)
        self._ensure(r, c)
        self._buf[r - 1, c - 1] = val

    def to_numpy(self) -> np.ndarray:
        out = self._buf[: self.nrows, : self.ncols].copy()
        if np.iscomplexobj(out) and not np.any(out.imag):
            out = out.real.copy()
        return out


class GrowVector:
    """MATLAB row vector with 1-based indexing and zero-padded growth."""

    __slots__ = ("_data",)

    def __init__(self, init=None):
        self._data = [] if init is None else [v for v in np.ravel(np.asarray(init)).tolist()]

    def __len__(self):
        return len(self._data)

    def __getitem__(self, i: int):
        if not 1 <= i <= len(self._data):
            raise IndexError(f"Index {i} exceeds vector length {len(self._data)}")
        return self._data[i - 1]

    def __setitem__(self, i: int, val):
        if i < 1:
            raise IndexError("MATLAB indices must be positive integers")
        if i > len(self._data):
            self._data.extend([0] * (i - len(self._data)))
        self._data[i - 1] = val

    def to_numpy(self, dtype=int) -> np.ndarray:
        return np.asarray(self._data, dtype=dtype)


class GrowColumns:
    """``M(1:k, n) = [...]`` column appends for patch-index matrices (k x n)."""

    __slots__ = ("_cols", "k")

    def __init__(self, k: int, init=None):
        self.k = k
        self._cols: list[tuple[int, ...]] = []
        if init is not None:
            arr = np.asarray(init, dtype=int)
            if arr.size:
                self._cols = [tuple(col) for col in arr.reshape(k, -1).T.tolist()]

    def set(self, n: int, values) -> None:
        if n < 1:
            raise IndexError("MATLAB indices must be positive integers")
        values = tuple(int(v) for v in values)
        while len(self._cols) < n:
            self._cols.append((0,) * self.k)
        self._cols[n - 1] = values

    def to_numpy(self) -> np.ndarray:
        if not self._cols:
            return np.zeros((self.k, 0), dtype=int)
        return np.asarray(self._cols, dtype=int).T.copy()


# --------------------------------------------------------------------------
# Output formatting
# --------------------------------------------------------------------------
def num2str(x) -> str:
    """MATLAB ``num2str`` for a real scalar or row vector of reals."""
    arr = np.atleast_1d(np.asarray(x, dtype=float))
    if np.all(np.isfinite(arr)) and np.all(arr == np.round(arr)):
        return "  ".join(f"{int(v):d}" for v in arr)
    finite = np.abs(arr[np.isfinite(arr)])
    xmax = finite.max() if finite.size else 0.0
    ndgt = math.floor(math.log10(xmax)) if xmax > 0 else 0
    ndgt = max(ndgt + 5, 5)
    ndgt = min(ndgt, 16)
    return "  ".join(f"{v:.{ndgt}g}" for v in arr)
