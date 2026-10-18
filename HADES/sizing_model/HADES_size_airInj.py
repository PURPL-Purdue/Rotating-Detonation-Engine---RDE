"""Port of ``HADES_size_airInj.m`` (air injector area sizing)."""

from __future__ import annotations

import math

__all__ = ["HADES_size_airInj"]


def HADES_size_airInj(ox_mdot, an_gap, an_rad_outter, temp):
    """Air injector area and annulus area.

    Parameters
    ----------
    ox_mdot : air mass flow [lbm/s]
    an_gap : annulus gap [mm]
    an_rad_outter : annulus outer radius [in]
    temp : air temperature [K]

    Returns
    -------
    (air_area, Area_Annulus) : both in mm^2
    """
    R = 8.314  # J/(mol K)
    Cd = 0.6  # discharge coefficient
    Delta_P = 1000000  # Pa
    ambient_pressure_air = 1000000  # Pa
    molarmass_air = 28.97  # g/mol

    # Total area air
    Air_mass_FR = ox_mdot * 0.453592  # lbm/s -> kg/s
    density = (ambient_pressure_air * molarmass_air) / (R * temp * 1000)  # kg/m^3
    air_area = (Air_mass_FR / (Cd * math.sqrt(2 * density * Delta_P))) * 1e6  # mm^2
    print("\nTotal Area Air: %.5f mm^2" % air_area)

    # Annulus area
    an_rad_outter_mm = an_rad_outter * 25.4  # in -> mm
    Area_Annulus = math.pi * (an_rad_outter_mm**2) - math.pi * (an_rad_outter_mm - an_gap) ** 2  # mm^2
    print("\nArea Annulus: %.5f mm^2" % Area_Annulus)

    # Air injector width
    air_inj_width = an_rad_outter_mm - math.sqrt(an_rad_outter_mm**2 - (air_area / math.pi))  # mm
    print("\nAir Injector Width: %.5f mm" % air_inj_width)
    return air_area, Area_Annulus
