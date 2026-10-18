"""Post-injection velocity from momentum-flux / area-change relations.

Python translation of ``RDE_Toolbox_2/src/MoC_Injection_Velocity.m``
(Blake's "Derivation of Momentum Flux-Area Change Relations").

Porting fix (see PORTING_NOTES.md): the MATLAB file does not parse (the last
``fprintf`` argument is ``out.``).  The 5th printed column ("Mach #") now
prints ``M4 = V4/sqrt(gamma1*R*T4)``; the "infill angle phi" header column
has no value in the original and is still not printed.
"""

import math

import numpy as np

from .structs import InjectionVelocityOut


def MoC_Injection_Velocity(InjV, ChamberDim):
    """Solve for the chamber-side velocity V4 (both roots), pressure P4 and temperature T4.

    SI units except ``InjV.injector_area`` [mm^2] and ``ChamberDim`` diameters [mm].
    Returns an :class:`~MODELING.moc.rde.structs.InjectionVelocityOut`.
    """
    P3 = InjV.P3
    T3 = InjV.T3
    mdot = InjV.mdot
    Pa = InjV.Pa
    R = InjV.R
    Cp = InjV.Cp
    N_holes = InjV.N_holes
    gamma1 = InjV.gamma1

    # Converting to SI and using ChamberDim to calculate A4
    A3 = InjV.injector_area * N_holes * 1e-6  # injector orifice area, mm^2 -> m^2
    A4 = (math.pi / 4) * (ChamberDim.outer_chamber_diameter**2 - ChamberDim.inner_chamber_diameter**2) * 1e-6

    V3 = mdot * R * T3 / (P3 * A3)

    T03 = T3 + V3**2 / (2 * Cp)
    beta = mdot / A4
    alpha = (P3 - Pa) * A3 / A4 + Pa + beta * V3

    a = beta * (R / (2 * Cp) - 1)
    b = alpha
    c = -beta * R * T03

    D = b**2 - 4 * a * c
    if D < 0:
        raise ValueError("No real velocity solutions for these inputs.")

    V4 = np.array([(-b + math.sqrt(D)) / (2 * a), (-b - math.sqrt(D)) / (2 * a)])
    P4 = alpha - beta * V4
    T4 = T03 - V4**2 / (2 * Cp)
    out = InjectionVelocityOut(V4=V4, P4=P4, T4=T4)

    # NaN marks a root with nonpositive velocity, pressure, or temperature.
    invalid = (out.V4 <= 0) | (out.P4 <= 0) | (out.T4 <= 0)
    # MATLAB ``out.Vwave(invalid) = NaN`` on a non-existent field: a row vector
    # long enough for the last invalid index, zero elsewhere.
    if np.any(invalid):
        last = int(np.flatnonzero(invalid)[-1]) + 1
        out.Vwave = np.zeros(last)
        out.Vwave[invalid[:last]] = np.nan
    else:
        out.Vwave = np.zeros(0)

    with np.errstate(invalid="ignore"):
        out.M4 = np.where(out.T4 > 0, out.V4 / np.sqrt(np.where(out.T4 > 0, gamma1 * R * out.T4, 1.0)), np.nan)

    print("Root     V4 (m/s)       P4 (Pa)        T4 (K)      Mach #   infill angle phi")
    for k in range(1, 1 + 1):
        print(
            "%d     %11.3f   %12.3f   %11.3f   %13.3f"
            % (k, out.V4[k - 1], out.P4[k - 1], out.T4[k - 1], out.M4[k - 1])
        )
    return out
