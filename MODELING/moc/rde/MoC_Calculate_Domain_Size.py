"""Unwrapped 2-D domain of the RDE annulus.

Python translation of ``RDE_Toolbox_2/src/MoC_Calculate_Domain_Size.m``.
"""

import math


def MoC_Calculate_Domain_Size(ChamberDim):
    """Return ``(xMax, yMax)`` [mm]: mean-diameter circumference and chamber length."""
    # Convert all variables to SI units (meters)
    OC_dia = ChamberDim.outer_chamber_diameter * 1e-3
    IC_dia = ChamberDim.inner_chamber_diameter * 1e-3
    Chamb_L = ChamberDim.chamber_length * 1e-3  # noqa: F841  (computed but unused in the MATLAB)

    # Annulus circumference and the unwrapped domain
    annulusMeanDiameter = (OC_dia + IC_dia) / 2
    annulusCircumference = math.pi * annulusMeanDiameter
    xMax = annulusCircumference * 1000
    yMax = ChamberDim.chamber_length
    return xMax, yMax
