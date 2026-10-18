"""Port of ``HADES_size_HoopStressTemps.m`` (failure temperature from hoop stress)."""

from __future__ import annotations

import numpy as np

from rdelib.matlab_compat import fzero

__all__ = ["HADES_size_HoopStressTemps", "TEMP_F", "YIELD_KSI"]

# Temperature data [F] and yield strength [ksi] (before the factor of safety)
TEMP_F = np.array([1500, 1440, 1330, 1240, 1180, 1140, 1100, 1020, 980, 922, 835, 719, 661,
                   589, 502, 459, 401, 357, 285, 242, 213, 184, 155, 126, 68], dtype=float)
YIELD_KSI = np.array([25.8, 29.5, 37.4, 42.7, 45.9, 48.1, 49.9, 52.6, 53.7, 54.8, 55.7, 56.3,
                      56.3, 56.6, 57, 57.4, 58.3, 59.2, 62.1, 64, 65.5, 67.5, 69.8, 72.5, 79],
                     dtype=float)


def HADES_size_HoopStressTemps(r_i_in, r_o_in, p_o_MPa, P_cj_bar, fos, *, plot: bool = True):
    """Wall failure temperature where the allowed pressure equals the CJ pressure.

    Parameters
    ----------
    r_i_in, r_o_in : inner / outer wall radius [in]
    p_o_MPa : outside pressure [MPa]
    P_cj_bar : CJ pressure [bar]
    fos : yield-strength factor of safety [-]

    Returns
    -------
    (Fail_temp_C, Fail_temp_F)
    """
    yeild_points_ksi = YIELD_KSI / fos
    temp_points_K = ((TEMP_F - 32) * 5 / 9) + 273.15
    yeild_points_MPa = yeild_points_ksi * 6.89476
    temp_points_C = temp_points_K - 273.15

    r_i_mm = r_i_in * 25.4
    r_o_mm = r_o_in * 25.4
    r_mm = r_i_mm
    u_mm = (r_o_mm**2) - (r_i_mm**2)
    a_mm = (r_i_mm**2 / u_mm) + (r_i_mm**2 * r_o_mm**2 / (r_mm**2 * u_mm))
    b_mm = (r_o_mm**2 / u_mm) + (r_i_mm**2 * r_o_mm**2 / (r_mm**2 * u_mm))

    p_i_MPa = ((p_o_MPa * b_mm) + yeild_points_MPa) / a_mm
    p_i_psi = p_i_MPa * 145.03773773
    P_cj_psi = P_cj_bar * 14.5038

    p = np.polyfit(temp_points_C, p_i_psi, 3)
    temp_points_ideal_C = np.linspace(temp_points_C.min(), temp_points_C.max(), 100)
    p_i_ideal_psi = np.polyval(p, temp_points_C)
    p_i_curve_psi = np.polyval(p, temp_points_ideal_C)
    p_i_psi_bar = np.mean(p_i_psi)

    SSE = np.sum((p_i_psi - p_i_ideal_psi) ** 2)
    SST = np.sum((p_i_psi - p_i_psi_bar) ** 2)
    R_squared = 1 - (SSE / SST)
    print("The R squared value is %.3f" % R_squared)

    intersection_func = lambda T: np.polyval(p, T) - P_cj_psi  # noqa: E731
    T_initial_guess_C = np.mean(temp_points_C)
    T_intersection_C = float(fzero(intersection_func, T_initial_guess_C))

    print("Failure Temperature occurs at %.2f C" % T_intersection_C)
    print("Failure Temperature occurs at %.2f F" % ((T_intersection_C * 1.8) + 32))

    Fail_temp_C = T_intersection_C
    Fail_temp_F = (T_intersection_C * 1.8) + 32

    if plot:
        import matplotlib.pyplot as plt

        plt.figure()
        plt.plot((temp_points_C * 1.8) + 32, p_i_psi, "r*")
        plt.plot((temp_points_ideal_C * 1.8) + 32, p_i_curve_psi, linewidth=1.6, color=(0.5, 0, 0.8))
        plt.axhline(P_cj_psi, linestyle="-", color="k")
        plt.text(plt.xlim()[0], P_cj_psi, "CJ Presure", va="bottom")
        plt.title("Maximum Pressure of SS316 vs Temperature")
        plt.xlabel("Temperature [F]")
        plt.ylabel("Allowed Chamber Pressure [psi]")
        plt.grid(True)
    return Fail_temp_C, Fail_temp_F
