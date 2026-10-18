"""Port of ``HADES_size_sonicnozzlearea_air.m`` (air sonic-nozzle choke diameter sweep).

Canonical version of four near-identical MATLAB scripts:
``HADES_size_sonicnozzlearea.m`` (differs only in an unused feed area),
``HADES_size_sonicnozzlearea_h2.m`` (byte-identical copy of the air script,
still air properties) and ``untitled.m`` (copy with merge-conflict markers).

BUG FIX (see PORTING_NOTES.md): the MATLAB used ``R = 8.134`` J/(mol K), a
typo for the universal gas constant 8.314. The default is now 8.314; pass
``R=8.134`` to reproduce the original numbers.
"""

from __future__ import annotations

import math

import numpy as np

from rdelib.matlab_compat import colon

__all__ = ["HADES_size_sonicnozzlearea_air", "main"]


def HADES_size_sonicnozzlearea_air(*, R: float = 8.314, plot: bool = True):
    """Choke diameter [mm] vs feed pressure (20..60 bar) for 2 lbm/s of air.

    Returns ``(p_feed_bar, dia_mm, A_choke_m2)``.
    """
    air_mdot_lb = 2  # lbm/s
    air_mdot = air_mdot_lb / 2.205  # kg/s

    p_feed_bar = colon(20, 1, 60)  # bar
    p_feed = p_feed_bar * 100000  # Pa

    c_d = 0.9
    c_f = 0.6989  # critical flow factor
    M = 0.02897  # kg/mol
    T_0 = 283.15  # K

    d_feed = 0.61 * 25.4  # mm (feed line ID, only used in commented-out prints)
    A_feed = math.pi * (d_feed / 2) ** 2 / 10**6  # m^2  # noqa: F841

    A_choke_better = air_mdot * np.sqrt(R / M * T_0) / (c_d * c_f * p_feed)
    dia_mm = 2 * np.sqrt(A_choke_better * 10**6 / math.pi)

    if plot:
        import matplotlib.pyplot as plt

        plt.figure()
        plt.plot(p_feed_bar, dia_mm)
        plt.xlabel("Feed Pressure (bar)")
        plt.ylabel("Choke Diameter (mm)")
        plt.axhline(11, color="r", linestyle="-")
        plt.text(p_feed_bar[-1], 11, "Minimum Diameter", color="r", ha="right", va="bottom")
        plt.title("Choke Diameter vs Feed Pressure")
        plt.grid(True)
    return p_feed_bar, dia_mm, A_choke_better


def main():
    import matplotlib.pyplot as plt

    HADES_size_sonicnozzlearea_air()
    plt.show()


if __name__ == "__main__":
    main()
