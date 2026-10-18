"""Port of ``HADES_size_geometry.m`` (annulus geometry / wave path length)."""

from __future__ import annotations

import math

__all__ = ["HADES_size_geometry"]


def HADES_size_geometry(outer_r, wall_thickness, annulus_gap):
    """Detonation wave path length around the annulus mid-radius.

    Parameters
    ----------
    outer_r : outer radius of the annulus [in]
    wall_thickness : wall thickness [in] (printed only)
    annulus_gap : annulus gap [mm]

    Returns
    -------
    det_wave_path_length : circumference at the annulus mid-radius [in]
    """
    annulus_gap_in = annulus_gap / 25.4  # mm -> in
    inner_r = outer_r - annulus_gap_in
    annulus_r = (outer_r + inner_r) / 2
    det_wave_path_length = 2 * math.pi * annulus_r

    print("\nGeometry:")
    print("Outer Radius: %.3f in." % outer_r)
    print("Wall Thickness: %.3f in." % wall_thickness)
    print("Annulus Gap: %.3f mm (%.3f in). " % (annulus_gap, annulus_gap_in))
    print("Inner Radius: %.3f in." % inner_r)
    print("Annulus Radius: %.3f in." % annulus_r)
    print("Detonation Wave Path Length: %.3f in." % det_wave_path_length)
    return det_wave_path_length
