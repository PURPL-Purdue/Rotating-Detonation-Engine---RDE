"""Small MATLAB-compatibility helpers shared by the HADES data-analysis ports.

Random numbers
--------------
The MATLAB scripts draw from MATLAB's global stream with ``rand``, ``randn`` and
``randi``.  Every Python port takes an ``rng`` argument instead:

* ``None``                      -> a fresh, unseeded ``numpy.random.Generator``
  (MATLAB scripts were run without a seed, so the data are random every run).
* ``int`` seed                  -> :func:`matlab_rng` (``numpy.random.RandomState``).
  Its uniform stream is bit-identical to MATLAB's ``rng(seed, 'twister');
  rand(...)`` (both are MT19937 with ``init_genrand`` and 53-bit
  ``genrand_res53``; MATLAB maps seed 0 to the reference seed 5489), so any
  script that only uses ``rand`` reproduces MATLAB exactly.  ``randn`` and
  ``randi`` cannot match MATLAB (MATLAB uses a Ziggurat normal generator and a
  different integer mapping).
* a ``Generator`` / ``RandomState`` (or any object exposing ``random``,
  ``standard_normal`` and ``integers``/``randint``) -> used as-is.

Draws are always made in the same order and with the same counts as the
MATLAB code (MATLAB fills ``rand(n,1)``/``rand(1,n)`` sequentially), so a
replayed stream gives identical results in both languages.

Candidates for later consolidation into ``rdelib``: ``matlab_round`` and the
random-stream adapter.
"""

from __future__ import annotations

from typing import Any

import numpy as np

RNGLike = Any  # None | int | numpy.random.Generator | numpy.random.RandomState | duck-typed


def matlab_rng(seed: int) -> np.random.RandomState:
    """Stream whose ``random`` values equal MATLAB ``rng(seed,'twister'); rand``.

    MATLAB ``rng(0)`` (also ``rng default``) seeds MT19937 with 5489, so e.g.
    ``matlab_rng(0).random(3)`` -> 0.8147, 0.9058, 0.1270 and
    ``matlab_rng(1).random()`` -> 0.4170.
    """
    seed = int(seed)
    return np.random.RandomState(5489 if seed == 0 else seed)


def as_rng(rng: RNGLike = None):
    """Normalise the ``rng`` argument accepted by the fake-data generators."""
    if rng is None:
        return np.random.default_rng()
    if isinstance(rng, (int, np.integer)):
        return matlab_rng(rng)
    return rng


def rand(rng, n: int | None = None):
    """MATLAB ``rand`` (scalar when ``n`` is None, else ``n`` values in order)."""
    return rng.random(n)


def randn(rng, n: int | None = None):
    """MATLAB ``randn`` (scalar when ``n`` is None, else ``n`` values in order)."""
    return rng.standard_normal(n)


def randi(rng, lo: int, hi: int, n: int | None = None):
    """MATLAB ``randi([lo hi])``: integers uniformly distributed on ``lo..hi`` inclusive."""
    if hasattr(rng, "integers"):
        out = rng.integers(lo, hi, size=n, endpoint=True)
    else:  # numpy.random.RandomState
        out = rng.randint(lo, hi + 1, size=n)
    return int(out) if n is None else np.asarray(out, dtype=int)


def matlab_round(x):
    """MATLAB ``round``: halves are rounded away from zero (Python rounds to even)."""
    x = np.asarray(x, dtype=float)
    r = np.sign(x) * np.floor(np.abs(x) + 0.5)
    return r.item() if r.ndim == 0 else r
