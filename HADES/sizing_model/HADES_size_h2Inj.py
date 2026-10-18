"""Port of ``HADES_size_h2Inj.m`` (hydrogen injector area sizing)."""

from __future__ import annotations

import math

from rdelib.matlab_compat import num2str

__all__ = ["HADES_size_h2Inj"]


def HADES_size_h2Inj(fuel_mdot, chamber_pressure):
    """Hydrogen injector area [mm^2] (incompressible orifice equation).

    Parameters
    ----------
    fuel_mdot : H2 mass flow [lbm/s]
    chamber_pressure : chamber pressure [psia]; the manifold is set at
        ``chamber_pressure * (1 + stiffness/100)`` with stiffness = 100 %.
    """
    stiffness = 100
    fuel_mdot = fuel_mdot * 0.453592  # lbm/s -> kg/s
    dischargeCoef = 0.7
    density = 0.91  # kg/m^3 ("UNSURE OF WHERE THIS NUM IS FROM" in the MATLAB)
    manifold_pressure = chamber_pressure * (1 + (stiffness / 100))  # psia
    pressureDrop = (manifold_pressure - chamber_pressure) * 6895  # Pa

    print("Given: ")
    print("gamma (specific heat ratio) = 1.4 ")
    print("R (specific gas constant for hydrogen in J/(kg*K) = 4124")

    h2_area = (fuel_mdot / dischargeCoef) * math.sqrt(1 / (2 * density * pressureDrop))
    h2_area = h2_area * 1000000  # m^2 -> mm^2

    # MATLAB: fprintf("\nHydrogen Area = " + h2_area + " mm^2 \n") (string concatenation)
    print("\nHydrogen Area = " + num2str(h2_area) + " mm^2 ")

    orificeSize = 1  # diameter, mm
    orificeArea = math.pi * (orificeSize / 2) ** 2
    orificeCount = h2_area / orificeArea
    # MATLAB: fprintf("Orifice count: " + orificeCount) -- no trailing newline
    print("Orifice count: " + num2str(orificeCount), end="")
    return h2_area
