"""Port of ``HADES_size_max_flux_range.m``."""

from __future__ import annotations

__all__ = ["HADES_size_max_flux_range"]


def HADES_size_max_flux_range(total_mdot, Area_Annulus):
    """Mass flux = total_mdot / Area_Annulus.

    Note: ``HADES_size_main`` passes ``total_mdot`` in lbm/s and the area in
    mm^2, so the result is in lbm/(s mm^2) even though the MATLAB label (kept
    verbatim) says kg/(s*mm^2). See PORTING_NOTES.md.
    """
    max_flux = total_mdot / Area_Annulus
    print("Max Flux Range: %.5f kg/(s*mm^2)" % max_flux)
    return max_flux
