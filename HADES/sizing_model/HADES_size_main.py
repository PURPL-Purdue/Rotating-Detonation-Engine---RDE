"""Port of ``HADES_size_main_v1.m`` - the HADES RDE sizing chain, end to end.

Run::

    uv run python -m HADES.sizing_model.HADES_size_main

``HADES_size_main(...)`` takes the MATLAB script's inputs as keyword
arguments (defaults = the values hard-coded in ``HADES_size_main_v1.m``),
prints the same report, draws the same figures (not shown unless ``main()``
is used) and returns every intermediate result in a :class:`MainResults`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .HADES_size_airInj import HADES_size_airInj
from .HADES_size_annulusgap_fillheight import HADES_size_annulusgap_fillheight
from .HADES_size_ceaDet import CeaDetResults, HADES_size_ceaDet
from .HADES_size_ceaRocket import CeaRocketResults, HADES_size_ceaRocket
from .HADES_size_chamberPressure import HADES_size_chamberPressure
from .HADES_size_convert import HADES_size_convert
from .HADES_size_geometry import HADES_size_geometry
from .HADES_size_h2Inj import HADES_size_h2Inj
from .HADES_size_HoopStressTemps import HADES_size_HoopStressTemps
from .HADES_size_max_flux_range import HADES_size_max_flux_range
from .HADES_size_P0_v1 import HADES_size_P0_v1
from .HADES_size_P_plenum_fanno import HADES_size_P_plenum_fanno
from .HADES_size_thrust import HADES_size_thrust

__all__ = ["HADES_size_main", "MainResults", "main"]


@dataclass
class MainResults:
    phi: float
    ox_mdot: float  # lbm/s
    fuel_mdot: float  # lbm/s
    total_mdot: float  # lbm/s
    cellWidth: float  # mm
    annulus_gap: float  # mm
    chosen_fillHeight: float  # m
    det_wave_path_length: float  # in
    ceaDet_results: CeaDetResults
    avg_chamber_p: float  # psia
    ceaRock_results: CeaRocketResults
    thrust: float  # N
    Failure_temps: float  # C (first output of HoopStressTemps, as in the MATLAB)
    P_0: float  # bar
    h2_area: float  # mm^2
    air_area: float  # mm^2
    area_annulus: float  # mm^2
    max_flux: float
    P_plenum_req: float  # Pa


def HADES_size_main(
    *,
    phi: float = 0.9,
    outer_radius: float = 2,  # in
    wall_thickness: float = 0.625,  # in
    wave_modes: int = 1,
    azi_annulus_length: float = 0.2893,  # m
    det_wave_num: int = 1,
    CJ_det_speed: float = 1931.5,  # m/s (hand input for the fill height)
    burned_gas_mach: float = 1075.6,  # burned-gas sound speed [m/s]
    burned_gas_p: float = 32.33,  # bar
    inj_crit_p: float = 2,  # bar
    unburned_axial_vel: float = 391,  # m/s
    surface_roughness_304: float = 3.2e-6,  # m
    factor_of_safety: float = 1.75,
    ambient_pressure_Mpa: float = 0.101325,
    ox_mdot: float = 2,  # lbm/s
    initial_temp: float = 283.15,  # K
    temp_units: str = "K",
    initial_pressure: float = 26,  # psia
    pressure_units: str = "psia",
    fuel_type: str = "H2",
    ox_type: str = "Air",
    plot: bool = True,
) -> MainResults:
    """Run the HADES sizing chain (``HADES_size_main_v1.m``)."""
    # (injector_dia = 1e-3 m and injector_length = 0.25 in are defined but
    #  unused in the MATLAB; the Fanno call below hard-codes its geometry.)
    fuel_mdot = ox_mdot * phi / 34.3  # lbm/s
    total_mdot = ox_mdot + fuel_mdot

    print("Inputs: ")
    print("Equivalence Ratio = %.1f" % phi)
    print("Initial Pressure: %.3f %s " % (initial_pressure, pressure_units))
    print("Initial Temp: %.2f %s \n" % (initial_temp, temp_units))

    print("Outputs: ")
    cellWidth, annulus_gap, chosen_fillHeight = HADES_size_annulusgap_fillheight(
        phi, azi_annulus_length, det_wave_num, CJ_det_speed, burned_gas_mach, burned_gas_p,
        inj_crit_p, unburned_axial_vel, plot=plot)

    det_wave_path_length = HADES_size_geometry(outer_radius, wall_thickness, annulus_gap)

    ceaDet_results = HADES_size_ceaDet(ox=ox_type, fuel=fuel_type, phi=phi, P0=initial_pressure,
                                       P0Units=pressure_units, T0=initial_temp, T0Units=temp_units)

    avg_chamber_p = HADES_size_chamberPressure(wave_modes, det_wave_path_length, ceaDet_results.cjVel,
                                               ceaDet_results.P_ratio, initial_pressure, plot=plot)

    ceaRock_results = HADES_size_ceaRocket(ox=ox_type, fuel=fuel_type, phi=phi, Pc=avg_chamber_p,
                                           PcUnits=pressure_units)

    thrust = HADES_size_thrust(ceaRock_results.isp, total_mdot)
    Failure_temps = HADES_size_HoopStressTemps(
        outer_radius - wall_thickness, outer_radius, ambient_pressure_Mpa,
        HADES_size_convert(initial_pressure * ceaDet_results.P_ratio, pressure_units, "bar"),
        factor_of_safety, plot=plot)[0]

    P_0 = HADES_size_P0_v1()

    h2_area = HADES_size_h2Inj(fuel_mdot, avg_chamber_p)

    air_area, area_annulus = HADES_size_airInj(ox_mdot, annulus_gap, outer_radius, initial_temp)

    max_flux = HADES_size_max_flux_range(total_mdot, area_annulus)

    P_plenum_req = HADES_size_P_plenum_fanno(
        10 * 100000, fuel_mdot / (2.205 * 32), (math.pi * 0.25 * 0.00118**2),
        HADES_size_convert(0.25, "in", "m") * 1e3, 1.18, initial_temp, 0.6, 1.4, 4124.2, 84e-7, 273,
        72, surface_roughness_304 * 1e3, 0.5)

    return MainResults(
        phi=phi, ox_mdot=ox_mdot, fuel_mdot=fuel_mdot, total_mdot=total_mdot, cellWidth=cellWidth,
        annulus_gap=annulus_gap, chosen_fillHeight=chosen_fillHeight,
        det_wave_path_length=det_wave_path_length, ceaDet_results=ceaDet_results,
        avg_chamber_p=avg_chamber_p, ceaRock_results=ceaRock_results, thrust=thrust,
        Failure_temps=Failure_temps, P_0=P_0, h2_area=h2_area, air_area=air_area,
        area_annulus=area_annulus, max_flux=max_flux, P_plenum_req=P_plenum_req,
    )


def main():
    import matplotlib.pyplot as plt

    HADES_size_main()
    plt.show()


if __name__ == "__main__":
    main()
