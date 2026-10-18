"""Port of ``HADES_size_chamberPressure.m`` (average chamber pressure)."""

from __future__ import annotations

import numpy as np

from rdelib.matlab_compat import colon

from .HADES_size_convert import HADES_size_convert

__all__ = ["HADES_size_chamberPressure"]


def HADES_size_chamberPressure(wave_mode, detwave_pathlength, cj_vel_ms, p_p1, p_0, *, plot: bool = True):
    """Time-averaged chamber pressure [psia] for an exponential blow-down
    from the CJ pressure back to ``p_0`` over one wave cycle.

    Parameters
    ----------
    wave_mode : number of detonation waves [-]
    detwave_pathlength : wave path length [in]
    cj_vel_ms : CJ velocity [m/s]
    p_p1 : CJ pressure ratio P/P1 [-]
    p_0 : initial (fill) pressure [psia]
    """
    cj_vel = HADES_size_convert(cj_vel_ms, "m/s", "ft/s")
    cycle_time = (detwave_pathlength / (cj_vel * 12)) * 10**3 / wave_mode  # ms
    exp_decay_const = np.log(p_p1) / (cycle_time / 1000)

    time = colon(0, 0.001, cycle_time)  # ms, MATLAB 0:0.001:cycle_time
    chamber_pressure = p_p1 * p_0 * np.exp(-1 * exp_decay_const * time / 1000)  # psia
    avg_chamber_p = float(np.mean(chamber_pressure))

    print("\nAverage Chamber Pressure: %.2f psia (%.3f bar). " % (avg_chamber_p, avg_chamber_p / 14.504))

    if plot:
        import matplotlib.pyplot as plt

        plt.figure()
        plt.plot(time, chamber_pressure)
        plt.title("Chamber Pressure vs Time")
        plt.xlabel("Time (ms)")
        plt.ylabel("Chmaber Pressure (psia)")
        plt.grid(True)
    return avg_chamber_p
