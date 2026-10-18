"""Injector mass flow from plenum pressure (port of ``PURPLE_MassFlow_from_PlenumPressure.m``).

The plenum static pressure is treated as the stagnation pressure upstream of
choked injectors at room temperature; density from the ideal-gas law and the
choked mass-flow equation ``mdot = Cd*A*sqrt(gamma*rho*P*K)`` give the H2 and
air mass flow rates, plotted in lb/s.

BUG FIX (see PORTING_NOTES.md): the MATLAB script read row 2 of
``PURPL_Fake_Plenum_pressure()`` as air and row 3 as hydrogen, but the
generator returns ``[time; h2; air]``.  The port reads row 2 as H2 and row 3
as air.  The mass-flow-vs-time figure was also titled "Static Plenum Pressure
Over Time" (copy/paste); it is now "Mass Flow Rate Over Time".

Run as a script::

    uv run python -m HADES.data_analysis.purpl_mass_flow_from_plenum_pressure
"""

from __future__ import annotations

import math

import numpy as np

from .choked_flow import choked_mass_flow
from .purpl_fake_data import PURPL_Fake_Plenum_pressure, _figure

C_d_h2 = 0.71  # Cd for hydrogen injector
C_d_air = 0.62  # Cd for air injector

temp = 293.15  # room temperature (K)

R_air = 287.05  # J/kg-K
R_h2 = 4124.2  # J/kg-K

A_h2 = math.pi * (0.001 / 2) ** 2 * 32  # area of the hydrogen injector: 32 x 1 mm holes (m^2)
A_air = 304.77 / 1000 / 1000  # area of the air injector (m^2)

gamma_h2 = 1.405  # specific heat ratio, hydrogen
gamma_air = 1.4  # specific heat ratio, air

# Defined but unused in the MATLAB script (density comes from the ideal-gas law)
stag_density_h2 = 0.08988  # stagnation density of hydrogen at 293 K (kg/m^3)
stag_density_air = 1.225  # density of air (kg/m^3)

KG_TO_LB = 2.20462


def plenum_mass_flow(air_plenum_pressure, h2_plenum_pressure):
    """Choked air and H2 injector mass flow [kg/s] from plenum pressures [Pa].

    Returns ``(mass_flow_air, mass_flow_h2)``.
    """
    air_plenum_pressure = np.asarray(air_plenum_pressure, dtype=float)
    h2_plenum_pressure = np.asarray(h2_plenum_pressure, dtype=float)
    density_air = air_plenum_pressure / (R_air * temp)  # ideal gas law
    density_h2 = h2_plenum_pressure / (R_h2 * temp)
    mass_flow_air = choked_mass_flow(C_d_air, A_air, gamma_air, density_air, air_plenum_pressure)
    mass_flow_h2 = choked_mass_flow(C_d_h2, A_h2, gamma_h2, density_h2, h2_plenum_pressure)
    return mass_flow_air, mass_flow_h2


def PURPLE_MassFlow_from_PlenumPressure(fake_plenum_pressure=None, *, plot: bool = True, rng=None):
    """Compute (and plot) mass flow from a ``[time; h2; air]`` plenum-pressure matrix.

    Defaults to :func:`PURPL_Fake_Plenum_pressure` data like the MATLAB script.
    Returns ``(time, mass_flow_air, mass_flow_h2)`` in s and kg/s.
    """
    if fake_plenum_pressure is None:
        fake_plenum_pressure = PURPL_Fake_Plenum_pressure(rng, plot=plot)
    fake_plenum_pressure = np.asarray(fake_plenum_pressure, dtype=float)

    time = fake_plenum_pressure[0, :]
    h2_plenum_pressure = fake_plenum_pressure[1, :]  # fixed: MATLAB read this row as air
    air_plenum_pressure = fake_plenum_pressure[2, :]  # fixed: MATLAB read this row as H2

    mass_flow_air, mass_flow_h2 = plenum_mass_flow(air_plenum_pressure, h2_plenum_pressure)

    if plot:  # mass flow rate in lb/s
        fig, ax = _figure(2)
        ax.plot(time, mass_flow_h2 * KG_TO_LB)
        ax.plot(time, mass_flow_air * KG_TO_LB)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Mass Flow Rate (lbs/s)")
        ax.set_title("Mass Flow Rate Over Time")
        ax.legend(["Hydrogen Plenum", "Air Plenum"])
        ax.grid(True)

    return time, mass_flow_air, mass_flow_h2


def main() -> None:
    import matplotlib.pyplot as plt

    PURPLE_MassFlow_from_PlenumPressure()
    plt.show()


if __name__ == "__main__":
    main()
