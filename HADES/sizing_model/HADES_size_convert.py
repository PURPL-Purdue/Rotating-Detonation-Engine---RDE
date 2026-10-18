"""Port of ``HADES_size_convert.m`` (unit conversions).

The implementation lives in :func:`rdelib.units.convert` so that other
projects can share it; this module keeps the HADES name/API.
"""

from __future__ import annotations

from rdelib.units import convert as _convert

__all__ = ["HADES_size_convert"]


def HADES_size_convert(val, from_unit: str, to_unit: str):
    """Convert ``val`` from ``from_unit`` to ``to_unit``.

    Pressure (bar, atm, psia, mmHg), speed (m/s, ft/s), temperature
    (K, C, F, R and aliases), distance (m, ft, in) and mass flow
    (kgs, lbms). See :func:`rdelib.units.convert`.
    """
    return _convert(val, from_unit, to_unit)
