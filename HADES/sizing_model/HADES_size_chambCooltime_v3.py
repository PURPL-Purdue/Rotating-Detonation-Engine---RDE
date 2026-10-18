"""Port of ``HADES_size_chambCooltime_v3.m`` (chamber cooldown strategies, v3).

Four scenarios of the two-node (air + wall) lumped model with a
temperature-dependent internal convection coefficient: baseline, N2 forced
cooling after 30 min, rough outer surface, rough + N2. Author of the
MATLAB: Jon DeSimone. Supersedes ``HADES_size_chambCooltime_v2.m`` (same
physics; v3 renamed variables and changed plot annotations).

Integration uses :func:`rdelib.ode.ode45` (MATLAB ``ode45`` algorithm) since
the reported times are first output samples below the targets.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from rdelib.ode import ode45

__all__ = [
    "HADES_size_chambCooltime_v3",
    "CooldownScenario",
    "findCoolingTime",
    "ternary",
    "twoNodeModel_tempDependent",
    "twoNodeModel_N2",
    "main",
]


@dataclass
class CooldownScenario:
    name: str
    t: np.ndarray  # [s]
    T_air: np.ndarray  # [K]
    T_wall: np.ndarray  # [K]
    t_95: float  # [s] first output time with T_wall <= 95 % target (NaN if never)
    T_95: float  # [K]
    t_90: float  # [s]
    T_90: float  # [K]


def findCoolingTime(t, T, T_target):
    """First sample with ``T <= T_target``: (idx, t_cool, T_cool) or NaNs.

    ``idx`` is 1-based like the MATLAB (NaN when not found).
    """
    idx = np.flatnonzero(np.asarray(T) <= T_target)
    if idx.size:
        i = int(idx[0])
        return i + 1, float(t[i]), float(T[i])
    return math.nan, math.nan, math.nan


def ternary(condition, true_val, false_val):
    return true_val if condition else false_val


def _h_inner(t_air, t_surr, t_wall_i, h_min, h_max):
    temp_ratio = max(0.0, (t_air - t_surr) / (t_wall_i - t_surr))
    return h_min + (h_max - h_min) * temp_ratio**0.5


def twoNodeModel_tempDependent(t, Y, t_surr, t_wall_i, m_air, cp_air, m_wall, cp_wall,
                               h_air_inner_min, h_air_inner_max, h_air_outer, as_wall_inner,
                               as_wall_outer, epsilon, sigma):
    """Baseline / rough-surface model: dY/dt for Y = [T_air, T_wall] (K)."""
    t_air, t_wall = Y[0], Y[1]
    h_air_inner = _h_inner(t_air, t_surr, t_wall_i, h_air_inner_min, h_air_inner_max)
    q_air_to_wall = h_air_inner * as_wall_inner * (t_air - t_wall)
    q_wall_conv = h_air_outer * as_wall_outer * (t_wall - t_surr)
    q_wall_rad = epsilon * sigma * as_wall_outer * (t_wall**4 - t_surr**4)
    return np.array([-q_air_to_wall / (m_air * cp_air),
                     (q_air_to_wall - q_wall_conv - q_wall_rad) / (m_wall * cp_wall)])


def twoNodeModel_N2(t, Y, t_surr, t_wall_i, m_air, cp_air, m_wall, cp_wall, h_air_inner_min,
                    h_air_inner_max, h_air_outer_base, h_air_n2, t_n2_start, as_wall_inner,
                    as_wall_outer, epsilon, sigma):
    """N2-cooling model: h_outer switches to ``h_air_n2`` at ``t >= t_n2_start``."""
    t_air, t_wall = Y[0], Y[1]
    h_air_inner = _h_inner(t_air, t_surr, t_wall_i, h_air_inner_min, h_air_inner_max)
    h_air_outer = h_air_n2 if t >= t_n2_start else h_air_outer_base
    q_air_to_wall = h_air_inner * as_wall_inner * (t_air - t_wall)
    q_wall_conv = h_air_outer * as_wall_outer * (t_wall - t_surr)
    q_wall_rad = epsilon * sigma * as_wall_outer * (t_wall**4 - t_surr**4)
    return np.array([-q_air_to_wall / (m_air * cp_air),
                     (q_air_to_wall - q_wall_conv - q_wall_rad) / (m_wall * cp_wall)])


def HADES_size_chambCooltime_v3(*, plot: bool = True) -> dict[str, CooldownScenario]:
    """Run the v3 cooldown script; returns the four scenarios keyed by name."""
    r_outer = 0.066675
    r_inner = 0.0508
    length = 0.093345
    rho_wall = 7930
    cp_wall = 500
    k_wall = 15

    t_wall_i_c = 1188.12
    t_surr_c = 10
    t_wall_i = t_wall_i_c + 273.15
    t_surr = t_surr_c + 273.15

    h_air_outer_base = 10
    h_air_inner_max = 25
    h_air_inner_min = 10
    rho_air = 0.3
    cp_air = 1150

    h_air_outer_rough = h_air_outer_base * 1.8
    h_air_n2 = 60
    t_n2_start = 0.5 * 3600

    epsilon = 0.6
    sigma = 5.67e-8

    v_wall = math.pi * length * (r_outer**2 - r_inner**2)
    m_wall = rho_wall * v_wall
    as_wall_outer = 2 * math.pi * r_outer * length
    as_wall_inner = 2 * math.pi * r_inner * length
    lc = v_wall / as_wall_outer

    v_air = math.pi * r_inner**2 * length
    m_air = rho_air * v_air

    bi = h_air_outer_base * lc / k_wall
    print("Bi = %.4f (lumped capacitance %s)" % (bi, ternary(bi < 0.1, "valid", "questionable")))

    tspan = [0, 3 * 3600]
    Y0 = [t_wall_i, t_wall_i]
    common = (t_surr, t_wall_i, m_air, cp_air, m_wall, cp_wall, h_air_inner_min, h_air_inner_max)
    geo = (as_wall_inner, as_wall_outer, epsilon, sigma)

    print("\n=== RUNNING SIMULATIONS ===\n")
    print("Scenario 1: Baseline")
    s_base = ode45(lambda t, Y: twoNodeModel_tempDependent(t, Y, *common, h_air_outer_base, *geo), tspan, Y0)
    print("Scenario 2: N2 cooling after %.1f minutes" % (t_n2_start / 60))
    s_n2 = ode45(lambda t, Y: twoNodeModel_N2(t, Y, *common, h_air_outer_base, h_air_n2, t_n2_start, *geo), tspan, Y0)
    print("Scenario 3: Rough surface (h_outer = %.1f W/m²-K)" % h_air_outer_rough)
    s_rough = ode45(lambda t, Y: twoNodeModel_tempDependent(t, Y, *common, h_air_outer_rough, *geo), tspan, Y0)
    print("Scenario 4: Rough surface + N2 cooling")
    s_combo = ode45(lambda t, Y: twoNodeModel_N2(t, Y, *common, h_air_outer_rough, h_air_n2, t_n2_start, *geo), tspan, Y0)

    t_target_95 = t_surr + 0.05 * (t_wall_i - t_surr)
    t_target_95_c = t_target_95 - 273.15
    t_target_90 = t_surr + 0.10 * (t_wall_i - t_surr)
    t_target_90_c = t_target_90 - 273.15

    out: dict[str, CooldownScenario] = {}
    for key, s in (("base", s_base), ("N2", s_n2), ("rough", s_rough), ("combo", s_combo)):
        _, t95, T95 = findCoolingTime(s.t, s.y[:, 1], t_target_95)
        _, t90, T90 = findCoolingTime(s.t, s.y[:, 1], t_target_90)
        out[key] = CooldownScenario(key, s.t, s.y[:, 0], s.y[:, 1], t95, T95, t90, T90)
    b, n, r, c = out["base"], out["N2"], out["rough"], out["combo"]

    print("\n=== COOLING TIME COMPARISON ===\n")
    print("To 95%% Cooling (%.1f°C):" % t_target_95_c)
    print("  Baseline:      %.2f hrs (%.1f min)" % (b.t_95 / 3600, b.t_95 / 60))
    print("  N2 Cooling:    %.2f hrs (%.1f min)" % (n.t_95 / 3600, n.t_95 / 60))
    print("  Rough Surface: %.2f hrs (%.1f min)" % (r.t_95 / 3600, r.t_95 / 60))
    print("  Rough + N2:    %.2f hrs (%.1f min)\n" % (c.t_95 / 3600, c.t_95 / 60))
    print("To 90%% Cooling (%.1f°C):" % t_target_90_c)
    print("  Baseline:      %.2f hrs (%.1f min)" % (b.t_90 / 3600, b.t_90 / 60))
    print("  N2 Cooling:    %.2f hrs (%.1f min)" % (n.t_90 / 3600, n.t_90 / 60))
    print("  Rough Surface: %.2f hrs (%.1f min)" % (r.t_90 / 3600, r.t_90 / 60))
    print("  Rough + N2:    %.2f hrs (%.1f min)\n" % (c.t_90 / 3600, c.t_90 / 60))
    print("Time Savings vs Baseline (to 90%):")
    for label, s in (("N2 Cooling:   ", n), ("Rough Surface:", r), ("Rough + N2:   ", c)):
        print("  %s %.1f min (%.1f%% faster)" % (label, (b.t_90 - s.t_90) / 60, 100 * (b.t_90 - s.t_90) / b.t_90))
    print("")

    if plot:
        _plot(out, t_target_95_c, t_target_90_c, t_n2_start, t_wall_i_c)
    return out


def _plot(out, t95c, t90c, t_n2_start, t_wall_i_c):
    import matplotlib.pyplot as plt

    green = (0, 0.6, 0)
    styles = {"base": ("b", "Baseline"), "N2": ("r", "N₂ Cooling"),
              "rough": (green, "Rough Surface"), "combo": ("m", "Rough + N₂")}
    off95 = {"base": 100, "N2": 80, "rough": 60, "combo": 40}
    off90 = {"base": 80, "N2": 60, "rough": 40, "combo": 20}
    plt.figure(figsize=(12, 7))
    for k, s in out.items():
        col, lab = styles[k]
        plt.plot(s.t / 3600, s.T_wall - 273.15, color=col, linewidth=2.5, label=lab)
    plt.plot([0, 3], [t95c, t95c], "--", linewidth=1, color=(0.5, 0.5, 0.5), label="95% Target")
    plt.plot([0, 3], [t90c, t90c], ":", linewidth=1.5, color=(0.3, 0.3, 0.3), label="90% Target")
    for k, s in out.items():
        col, _ = styles[k]
        if not math.isnan(s.t_95):
            plt.plot(s.t_95 / 3600, s.T_95 - 273.15, "o", markersize=8, color=col, markerfacecolor=col)
            plt.text(s.t_95 / 3600, s.T_95 - 273.15 + off95[k], "%.0f min" % (s.t_95 / 60), fontsize=8,
                     color=col, ha="center", bbox=dict(facecolor="white", edgecolor=col, pad=1))
        if not math.isnan(s.t_90):
            plt.plot(s.t_90 / 3600, s.T_90 - 273.15, "s", markersize=11, color=col, markerfacecolor=col,
                     markeredgewidth=1.5)
            plt.text(s.t_90 / 3600, s.T_90 - 273.15 + off90[k], "%.0f min" % (s.t_90 / 60), fontsize=9,
                     color=col, fontweight="bold", ha="center",
                     bbox=dict(facecolor="white", edgecolor=col, pad=2))
    plt.axvline(t_n2_start / 3600, linestyle="--", linewidth=1, alpha=0.4, color=(0.7, 0.3, 0.3))
    plt.text(t_n2_start / 3600, 1150, "N₂ Start", fontsize=9, color=(0.7, 0.3, 0.3), ha="center", rotation=90)
    plt.plot(0, t_wall_i_c, "ko", markersize=8, markerfacecolor=(0.3, 0.3, 0.3))
    plt.grid(True)
    plt.xlabel("Time (hours)", fontsize=13, fontweight="bold")
    plt.ylabel("Wall Temperature (°C)", fontsize=13, fontweight="bold")
    plt.title("RDE Chamber Cooldown Strategies", fontsize=15, fontweight="bold")
    plt.legend(loc="upper right", fontsize=10)
    plt.ylim(0, 1250)
    plt.xlim(0, 3)


def main():
    import matplotlib.pyplot as plt

    HADES_size_chambCooltime_v3()
    plt.show()


if __name__ == "__main__":
    main()
