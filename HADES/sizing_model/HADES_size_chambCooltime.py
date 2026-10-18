"""Port of ``HADES_size_chambCooltime.m`` (chamber cooldown, two-node model, v1).

PURPL RDE chamber cooldown time with a two-node (air + wall) lumped model,
constant convection coefficients. Author of the MATLAB: Jon DeSimone.
The integration uses :func:`rdelib.ode.ode45`, a reimplementation of
MATLAB's ``ode45`` (same step control and Refine=4 output grid), because the
reported cooldown time is the first *output sample* below the target.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from rdelib.ode import ode45

__all__ = ["HADES_size_chambCooltime", "CooltimeResult", "twoNodeModel", "main"]


@dataclass
class CooltimeResult:
    t: np.ndarray  # [s]
    T_air: np.ndarray  # [K]
    T_wall: np.ndarray  # [K]
    TTarget_C: float
    tCool95_hours: float
    tCool95_minutes: float
    T95_wall_C: float
    T95_air_C: float
    Bi: float
    tau_air: float
    C_wall: float
    C_air: float
    q_conv_outer_initial: float
    q_rad_initial: float
    q_total_initial: float


def twoNodeModel(t, Y, TSurr, mAir, cpAir, mWall, cpWall, hAir_inner, hAir_outer,
                 AsWall_inner, AsWall_outer, epsilon, sigma):
    """dY/dt for Y = [T_air, T_wall] (K)."""
    T_air, T_wall = Y[0], Y[1]
    q_air_to_wall = hAir_inner * AsWall_inner * (T_air - T_wall)
    q_wall_conv = hAir_outer * AsWall_outer * (T_wall - TSurr)
    q_wall_rad = epsilon * sigma * AsWall_outer * (T_wall**4 - TSurr**4)
    dT_air_dt = -q_air_to_wall / (mAir * cpAir)
    dT_wall_dt = (q_air_to_wall - q_wall_conv - q_wall_rad) / (mWall * cpWall)
    return np.array([dT_air_dt, dT_wall_dt])


def HADES_size_chambCooltime(*, plot: bool = True) -> CooltimeResult:
    """Run the v1 cooldown script (all inputs hard-coded as in the MATLAB)."""
    # Solid wall geometry
    ro = 0.066675  # outer wall radius [m], 2.625"
    ri = 0.0508  # inner wall radius [m], 2"
    L = 0.093345  # axial length [m], 3.675"
    rhoWall = 7930  # SS316 density [kg/m^3]
    cpWall = 500  # SS316 specific heat [J/kg-K]
    kWall = 15  # SS316 conductivity [W/m-K]

    # Temperatures
    TWalli_C = 1188.12
    TSurr_C = 10
    TWalli = TWalli_C + 273.15
    TSurr = TSurr_C + 273.15

    # Air properties
    hAir_outer = 10
    hAir_inner = 10
    rhoAir = 0.3
    cpAir = 1150

    # Radiation
    epsilon = 0.6
    sigma = 5.67e-8

    VWall = math.pi * L * (ro**2 - ri**2)
    mWall = rhoWall * VWall
    AsWall_outer = 2 * math.pi * ro * L
    AsWall_inner = 2 * math.pi * ri * L
    Lc = VWall / AsWall_outer

    VAir = math.pi * ri**2 * L
    mAir = rhoAir * VAir

    C_wall = mWall * cpWall
    C_air = mAir * cpAir
    print("\n=== THERMAL MASS ANALYSIS ===")
    print("Wall thermal mass: %.2f J/K" % C_wall)
    print("Air thermal mass: %.2f J/K" % C_air)
    print("Ratio (C_wall/C_air): %.1f" % (C_wall / C_air))
    print("Wall will cool %.1fx slower than air\n" % (C_wall / C_air))

    Bi = hAir_outer * Lc / kWall
    if Bi < 0.1:
        print("Bi = %.4f < 0.1, lumped capacitance valid for wall" % Bi)
    else:
        print("Bi = %.4f > 0.1, lumped capacitance may be invalid" % Bi)

    tau_air = (mAir * cpAir) / (hAir_inner * AsWall_inner)
    print("Air equilibration time constant: %.2f seconds\n" % tau_air)

    def ode(t, Y):
        return twoNodeModel(t, Y, TSurr, mAir, cpAir, mWall, cpWall, hAir_inner, hAir_outer,
                            AsWall_inner, AsWall_outer, epsilon, sigma)

    tspan = [0, 3 * 3600]
    Y0 = [TWalli, TWalli]
    sol = ode45(ode, tspan, Y0)
    t, Y = sol.t, sol.y
    T_air = Y[:, 0]
    T_wall = Y[:, 1]
    t_hours = t / 3600
    T_air_celsius = T_air - 273.15
    T_wall_celsius = T_wall - 273.15

    TTarget = TSurr + 0.05 * (TWalli - TSurr)
    TTarget_C = TTarget - 273.15

    idx = np.flatnonzero(T_wall <= TTarget)
    if idx.size:
        i = idx[0]
        tCool95 = t[i]
        tCool95_hours = tCool95 / 3600
        tCool95_minutes = tCool95 / 60
        T95_wall_celsius = T_wall_celsius[i]
        T95_air_celsius = T_air_celsius[i]
        print("\n=== COOLING TIME TO 95% OF ORIGINAL TEMPERATURE ===")
        print("Time to cool wall to %.2f°C: %.2f hours (%.2f minutes)" % (TTarget_C, tCool95_hours, tCool95_minutes))
        print("Wall temperature at this time: %.2f°C" % T95_wall_celsius)
        print("Air temperature at this time: %.2f°C\n" % T95_air_celsius)
    else:
        print("Wall did not reach target temperature within 3 hours.")
        print("Final wall temperature after 3 hours: %.2f°C" % T_wall_celsius[-1])
        print("Final air temperature after 3 hours: %.2f°C\n" % T_air_celsius[-1])
        tCool95_hours = tCool95_minutes = T95_wall_celsius = T95_air_celsius = math.nan

    q_conv_outer_initial = hAir_outer * AsWall_outer * (TWalli - TSurr)
    q_rad_initial = epsilon * sigma * AsWall_outer * (TWalli**4 - TSurr**4)
    q_total_initial = q_conv_outer_initial + q_rad_initial

    print("=== HEAT TRANSFER ANALYSIS AT INITIAL TEMPERATURE ===")
    print("Initial wall temperature: %.2f°C" % TWalli_C)
    print("Initial air temperature: %.2f°C" % TWalli_C)
    print("Ambient temperature: %.2f°C" % TSurr_C)
    print("External convective heat loss: %.2f W (%.1f%%)" % (q_conv_outer_initial, 100 * q_conv_outer_initial / q_total_initial))
    print("Radiative heat loss: %.2f W (%.1f%%)" % (q_rad_initial, 100 * q_rad_initial / q_total_initial))
    print("Total initial heat loss: %.2f W\n" % q_total_initial)

    if plot:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(10, 7))
        plt.plot(t_hours, T_wall_celsius, "b-", linewidth=2.5, label="Wall Temperature")
        plt.plot([0, t_hours.max()], [TTarget_C, TTarget_C], "r--", linewidth=1.5, label="95% Cooldown Target")
        if not math.isnan(tCool95_hours):
            plt.plot(tCool95_hours, T95_wall_celsius, "ro", markersize=12, markerfacecolor="r", label="Target Reached")
            plt.text(tCool95_hours * 1.05, T95_wall_celsius,
                     "  %.2f hrs (%.1f min)\n  %.1f°C" % (tCool95_hours, tCool95_minutes, T95_wall_celsius),
                     fontsize=11, color="r", fontweight="bold")
        plt.plot(0, TWalli_C, "go", markersize=12, markerfacecolor="g", label="Initial Temp")
        plt.text(0.05, TWalli_C, "  Initial: %.1f°C" % TWalli_C, fontsize=11, color="g", fontweight="bold")
        plt.plot(t_hours.max(), T_wall_celsius[-1], "mo", markersize=10, markerfacecolor="m", label="Final Temp")
        plt.text(t_hours.max() * 0.95, T_wall_celsius[-1], "%.1f°C  " % T_wall_celsius[-1], fontsize=10,
                 color="m", ha="right", fontweight="bold")
        plt.grid(True)
        plt.xlabel("Time (hours)", fontsize=13, fontweight="bold")
        plt.ylabel("Wall Temperature (°C)", fontsize=13, fontweight="bold")
        plt.title("RDE Chamber Cooldown: Two-Node Model", fontsize=15, fontweight="bold")
        plt.legend(loc="upper right", fontsize=11)
        plt.ylim(0, max(TWalli_C * 1.05, 1300))
        plt.xlim(0, 3)

    print("=== SIMULATION COMPLETE ===")
    print("Final wall temperature after 3 hours: %.2f°C" % T_wall_celsius[-1])
    print("Final air temperature after 3 hours: %.2f°C" % T_air_celsius[-1])
    print("Wall temperature drop: %.2f°C" % (TWalli_C - T_wall_celsius[-1]))
    print("Air temperature drop: %.2f°C" % (TWalli_C - T_air_celsius[-1]))

    return CooltimeResult(
        t=t, T_air=T_air, T_wall=T_wall, TTarget_C=TTarget_C, tCool95_hours=tCool95_hours,
        tCool95_minutes=tCool95_minutes, T95_wall_C=T95_wall_celsius, T95_air_C=T95_air_celsius,
        Bi=Bi, tau_air=tau_air, C_wall=C_wall, C_air=C_air,
        q_conv_outer_initial=q_conv_outer_initial, q_rad_initial=q_rad_initial,
        q_total_initial=q_total_initial,
    )


def main():
    import matplotlib.pyplot as plt

    HADES_size_chambCooltime()
    plt.show()


if __name__ == "__main__":
    main()
