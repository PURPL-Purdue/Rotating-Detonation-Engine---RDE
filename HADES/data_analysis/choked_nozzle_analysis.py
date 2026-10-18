"""Choked-flow nozzle sweeps for multiple equivalence ratios.

Port of two MATLAB scripts that share every constant and the same equation
``mdot = Cd*A*sqrt(gamma*rho0*P0*(2/(gamma+1))^((gamma+1)/(gamma-1)))``
solved for the stagnation pressure ``P0`` needed to choke the design flow:

* ``cD_analysis.m``   -> :func:`cD_analysis`   : P0 vs discharge coefficient (T0 = 283 K)
* ``temp_analysis.m`` -> :func:`temp_analysis` : P0 vs stagnation temperature (Cd = 0.7)

Each function returns the computed curves and draws the same two-panel figure
(air on the left, one H2 curve per equivalence ratio on the right).

Run as a script::

    uv run python -m HADES.data_analysis.choked_nozzle_analysis          # both figures
    uv run python -m HADES.data_analysis.choked_nozzle_analysis cd       # only cD_analysis
    uv run python -m HADES.data_analysis.choked_nozzle_analysis temp     # only temp_analysis
"""

from __future__ import annotations

import sys
from dataclasses import dataclass

import numpy as np

from .choked_flow import mass_flow_factor

# ---------------- Given constants (shared by both MATLAB scripts) ----------------
mdot_air = 2 * 0.453592  # 2 lbm/s converted to kg/s
FAR_stoich = 0.0292  # stoichiometric fuel-to-air ratio for H2 (kg H2 / kg air)
phi_values = np.linspace(0.8, 1.2, 5)  # MATLAB 0.8:0.1:1.2
mdot_h2_vec = phi_values * FAR_stoich * mdot_air  # mdot_h2 = phi * FAR_stoich * mdot_air [kg/s]

# Throat areas: in^2 -> m^2 (1 in^2 = 6.4516e-4 m^2)
A_air = 0.472 * 6.4516e-4  # air throat area [m^2]
A_h2 = 0.03875 * 6.4516e-4  # H2 throat area [m^2]

# Gas properties
gamma_air = 1.4
R_air = 287.0  # J/(kg K)
gamma_h2 = 1.405
R_h2 = 4157.0  # J/(kg K)

T0_CD_SWEEP = 283.0  # stagnation temperature for cD_analysis [K]
CD_TEMP_SWEEP = 0.7  # fixed discharge coefficient for temp_analysis

# MATLAB ``lines`` colormap (default axes colour order)
MATLAB_LINES = np.array(
    [
        [0.0000, 0.4470, 0.7410],
        [0.8500, 0.3250, 0.0980],
        [0.9290, 0.6940, 0.1250],
        [0.4940, 0.1840, 0.5560],
        [0.4660, 0.6740, 0.1880],
        [0.3010, 0.7450, 0.9330],
        [0.6350, 0.0780, 0.1840],
    ]
)


@dataclass
class NozzleSweep:
    """Sweep variable ``x`` (Cd or T0) and the stagnation pressures in bar."""

    x: np.ndarray
    MFF_air: np.ndarray
    MFF_h2: np.ndarray
    P0_air_bar: np.ndarray  # shape (len(x),)
    P0_h2_bar: np.ndarray  # shape (len(phi_values), len(x))
    phi_values: np.ndarray
    mdot_h2_vec: np.ndarray


def _plot(sweep: NozzleSweep, xlabel, air_title, h2_title, air_sub, h2_sub, suptitle):
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.plot(sweep.x, sweep.P0_air_bar, "b-", linewidth=2.5)
    for ax, title, sub in ((ax1, air_title, air_sub), (ax2, h2_title, h2_sub)):
        ax.set_title(title, fontsize=14, fontweight="bold", pad=20)
        ax.text(0.5, 1.01, sub, transform=ax.transAxes, ha="center", va="bottom", fontsize=10)

    for i, (phi, mdot) in enumerate(zip(sweep.phi_values, sweep.mdot_h2_vec)):
        ax2.plot(
            sweep.x,
            sweep.P0_h2_bar[i],
            color=MATLAB_LINES[i % len(MATLAB_LINES)],
            linewidth=2.5,
            label=r"$\phi$ = %.1f  ($\dot{m}$ = %.4f kg/s)" % (phi, mdot),
        )
    ax2.legend(loc="best", fontsize=10)

    for ax in (ax1, ax2):
        ax.set_xlabel(xlabel, fontsize=13)
        ax.set_ylabel("Stagnation Pressure $P_0$ [bar]", fontsize=13)
        ax.grid(True, which="major")
        ax.minorticks_on()
        ax.grid(True, which="minor", alpha=0.3)
        ax.tick_params(labelsize=11)

    fig.suptitle(suptitle, fontsize=15, fontweight="bold")
    fig.tight_layout()
    return fig


def cD_analysis(*, plot: bool = True) -> NozzleSweep:
    """P0 required to choke the design flows vs discharge coefficient (0.5-0.9), T0 = 283 K."""
    T0 = T0_CD_SWEEP
    Cd = np.linspace(0.5, 0.9, 200)

    MFF_air = mass_flow_factor(gamma_air, R_air, T0)
    MFF_h2 = mass_flow_factor(gamma_h2, R_h2, T0)

    P0_air_bar = mdot_air / (Cd * A_air * MFF_air) * 1e-5
    P0_h2_bar = np.array([m / (Cd * A_h2 * MFF_h2) * 1e-5 for m in mdot_h2_vec])

    sweep = NozzleSweep(Cd, MFF_air, MFF_h2, P0_air_bar, P0_h2_bar, phi_values, mdot_h2_vec)
    if plot:
        _plot(
            sweep,
            "Discharge Coefficient $C_D$",
            "Air: $C_D$ vs. $P_0$",
            "H$_2$: $C_D$ vs. $P_0$",
            r"$\dot{m}$ = %.3f kg/s,  $T_0$ = %d K,  A = %.6f m$^2$" % (mdot_air, T0, A_air),
            r"$T_0$ = %d K,  A = %.8f m$^2$,  FAR$_{stoich}$ = %.4f" % (T0, A_h2, FAR_stoich),
            "Choked Nozzle Flow: $C_D$ vs. Stagnation Pressure",
        )
    return sweep


def temp_analysis(*, plot: bool = True) -> NozzleSweep:
    """P0 required to choke the design flows vs stagnation temperature (273-300 K), Cd = 0.7."""
    Cd = CD_TEMP_SWEEP
    T0 = np.linspace(273, 300, 200)

    MFF_air = mass_flow_factor(gamma_air, R_air, T0)  # vector over T0
    MFF_h2 = mass_flow_factor(gamma_h2, R_h2, T0)

    P0_air_bar = mdot_air / (Cd * A_air * MFF_air) * 1e-5
    P0_h2_bar = np.array([m / (Cd * A_h2 * MFF_h2) * 1e-5 for m in mdot_h2_vec])

    sweep = NozzleSweep(T0, MFF_air, MFF_h2, P0_air_bar, P0_h2_bar, phi_values, mdot_h2_vec)
    if plot:
        _plot(
            sweep,
            "Stagnation Temperature $T_0$ [K]",
            "Air: $T_0$ vs. $P_0$",
            "H$_2$: $T_0$ vs. $P_0$",
            r"$\dot{m}$ = %.3f kg/s,  $C_D$ = %.1f,  A = %.6f m$^2$" % (mdot_air, Cd, A_air),
            r"$C_D$ = %.1f,  A = %.8f m$^2$,  FAR$_{stoich}$ = %.4f" % (Cd, A_h2, FAR_stoich),
            "Choked Nozzle Flow: Stagnation Pressure vs. Temperature",
        )
    return sweep


def main(argv: list[str] | None = None) -> None:
    import matplotlib.pyplot as plt

    argv = sys.argv[1:] if argv is None else argv
    which = argv[0].lower() if argv else "both"
    if which in ("cd", "both"):
        cD_analysis()
    if which in ("temp", "both"):
        temp_analysis()
    plt.show()


if __name__ == "__main__":
    main()
