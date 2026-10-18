"""MATLAB ``inputParser``-style name/value handling for the CEA wrappers."""

from __future__ import annotations

import warnings


def resolve_name(name, names, func: str) -> str:
    """Match ``name`` like MATLAB ``inputParser`` (CaseSensitive=false, PartialMatching=true).

    Exact (case-insensitive) match first, otherwise a unique case-insensitive
    prefix. Unknown or ambiguous names raise ``TypeError``.
    """
    if not isinstance(name, str) or not name:
        raise TypeError(f"{func}: parameter names must be non-empty strings, got {name!r}")
    low = name.lower()
    exact = [n for n in names if n.lower() == low]
    if exact:
        return exact[0]
    pref = [n for n in names if n.lower().startswith(low)]
    if len(pref) == 1:
        return pref[0]
    if not pref:
        raise TypeError(f"{func}: unknown parameter {name!r}")
    raise TypeError(f"{func}: ambiguous parameter {name!r} matches {sorted(pref)}")


def parse_name_value(args: tuple, kwargs: dict, defaults: dict, func: str) -> dict:
    """Merge MATLAB-style ``('Name', value, ...)`` pairs and keyword args.

    Names are resolved like MATLAB ``inputParser`` defaults: case-insensitive,
    exact match first, then unambiguous prefix ('fu' -> 'fuel'; 'P0' is
    exact even though 'P0Units' starts with it). Unknown or ambiguous names
    raise ``TypeError``. A parameter given twice takes the last value.
    """
    if len(args) % 2:
        raise TypeError(f"{func}: name/value arguments must come in pairs")
    opts = dict(defaults)
    pairs = list(zip(args[0::2], args[1::2])) + list(kwargs.items())
    for name, value in pairs:
        opts[resolve_name(name, list(defaults), func)] = value
    if opts.get("ceaExe") is not None:
        warnings.warn(
            f"{func}: 'ceaExe' is ignored; CEA now runs in-process through the "
            "cross-platform 'cea' Python package.",
            stacklevel=3,
        )
    return opts


def deck_float(x: float) -> float:
    """Value as written to a CEA deck with MATLAB ``fprintf('%f')`` (6 decimals)."""
    return float("%f" % x)
