"""Sonic-orifice discharge coefficient from a supply-tank blowdown test.

Consolidates the two HADES feed-sizing scripts that read
``orifice_blowdown_input.csv``:

* ``blowdown_orifice_cd.py`` - run-average tank mass flow divided by the
  instantaneous ideal (Cd = 1) choked mass flow at every sample, giving
  ``Cds`` (Cd vs. time) and their arithmetic mean ``average_Cd``.
* ``blowdown_orifice_cd_chat.py`` - integrates the ideal choked flow
  (trapezoidal rule) and fits ONE constant ``Cd`` to the measured cumulative
  tank mass loss by least squares through the origin; also reports
  ``Cd_integrated = total measured mass loss / total ideal mass loss``.

Both scripts shared the same preprocessing and ideal-flow model, which is done
once here:

* CSV pressures are gauge [psi]; ``P_abs = (psi + 14.6959) * 6894.75729`` Pa.
* Tank mass = CoolProp real-gas density at the supply P, T times the tank volume.
* Ideal choked flow per sample uses the SUPPLY P, T as stagnation conditions,
  ``gamma = cp/cv`` from CoolProp at those conditions and ``R = R_u / M``.

Invalid samples (non-finite or non-positive time span, supply pressure or
temperature, or a state CoolProp cannot evaluate) raise ``ValueError`` naming the
offending rows, as the original per-sample CoolProp calls did.  The script
prints and plots only; it writes no files.

Input CSV columns (see ``orifice_blowdown_input_template.csv``)::

    time_ms, supply_pressure_psi, upstream_pressure_psi, downstream_pressure_psi,
    supply_temperature_K, upstream_temperature_K, downstream_temperature_K

Run from the repository root::

    uv run python -m HADES.feed_sizing.blowdown_orifice_cd
    uv run python -m HADES.feed_sizing.blowdown_orifice_cd --input test.csv --tank-volume-L 50 --orifice-diameter-in 0.24 --no-show
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import CoolProp.CoolProp as cp
import numpy as np
import pandas as pd
from rdelib.units import M_PER_IN as IN_TO_M  # 0.0254 m/in

MODULE_DIR = Path(__file__).resolve().parent

#: Default blowdown test data (resolved relative to this module, not the CWD).
BLOWDOWN_INPUT_CSV = MODULE_DIR / "orifice_blowdown_input.csv"

# Conversion factors (as in the original scripts; this psi->Pa value is specific
# to the blowdown scripts, see PORTING_NOTES.md issue 3)
PSI_TO_PA = 6894.75729
L_TO_M3 = 0.001
P_ATM_PSI = 14.6959

REQUIRED_COLUMNS = (
    "time_ms",
    "supply_pressure_psi",
    "upstream_pressure_psi",
    "downstream_pressure_psi",
    "supply_temperature_K",
    "upstream_temperature_K",
    "downstream_temperature_K",
)

# Default test article (original script inputs)
DEFAULT_GAS = "Air"
DEFAULT_TANK_VOLUME_L = 50.0
DEFAULT_ORIFICE_DIAMETER_IN = 0.24


def read_blowdown_csv(input_file: str | Path = BLOWDOWN_INPUT_CSV) -> pd.DataFrame:
    """Read blowdown test data and check the required columns are present."""
    data = pd.read_csv(input_file, float_precision="round_trip")
    missing = [c for c in REQUIRED_COLUMNS if c not in data.columns]
    if missing:
        raise ValueError(f"Blowdown CSV is missing required columns: {missing}")
    return data


def _row_list(bad: np.ndarray) -> str:
    rows = np.flatnonzero(bad)
    shown = ", ".join(str(int(i)) for i in rows[:20])
    more = f" ... ({len(rows)} rows)" if len(rows) > 20 else ""
    return f"data row index(es) [{shown}{more}] (0-based; CSV line = index + 2)"


def _require(ok: np.ndarray, message: str) -> None:
    """Raise ``ValueError`` listing the rows where ``ok`` is False."""
    bad = ~np.asarray(ok, dtype=bool)
    if bad.any():
        raise ValueError(f"{message}; invalid at {_row_list(bad)}.")


def _numeric_column(data: pd.DataFrame, column: str) -> np.ndarray:
    return pd.to_numeric(data[column], errors="coerce").to_numpy(dtype=float)


def _scalar_propssi_or_nan(output: str, p: float, t: float, gas: str) -> float:
    try:
        return cp.PropsSI(output, "P", float(p), "T", float(t), gas)
    except ValueError:
        return np.nan


def _coolprop(output: str, pressure_pa: np.ndarray, temperature_k: np.ndarray, gas: str) -> np.ndarray:
    """Vectorized ``PropsSI(output, 'P', p, 'T', T, gas)`` that rejects non-finite results.

    Scalar ``PropsSI`` (used by the original scripts) raises on an invalid state;
    the vectorized call returns ``inf`` for that sample instead, so check it.  If
    NO sample can be evaluated the vectorized call raises ("No outputs were able
    to be calculated"); then fall back to per-sample calls to find the bad rows.
    """
    try:
        values = np.asarray(cp.PropsSI(output, "P", pressure_pa, "T", temperature_k, gas), dtype=float)
    except ValueError:
        values = np.array(
            [_scalar_propssi_or_nan(output, p, t, gas) for p, t in zip(pressure_pa, temperature_k)],
            dtype=float,
        )
    bad = ~np.isfinite(values) | (values <= 0)
    if bad.any():
        i = int(np.flatnonzero(bad)[0])
        p, t = float(pressure_pa[i]), float(temperature_k[i])
        try:
            cp.PropsSI(output, "P", p, "T", t, gas)
            reason = f"returned {values[i]!r}"
        except ValueError as exc:
            reason = str(exc)
        raise ValueError(
            f"CoolProp could not evaluate {output} for {gas} at {_row_list(bad)}; "
            f"first bad row {i}: P = {p:.6g} Pa, T = {t:.6g} K: {reason}"
        )
    return values


@dataclass
class BlowdownCdResult:
    """Outputs of :func:`blowdown_orifice_cd` (original script variable names, SI units)."""

    gas: str
    tank_volume: float  # [m^3]
    orifice_diameter: float  # [in]
    orifice_area: float  # [m^2]
    R: float  # specific gas constant from CoolProp [J/(kg*K)]
    time_ms: np.ndarray
    time_s: np.ndarray
    supply_pressure: np.ndarray  # absolute [Pa]
    upstream_pressure: np.ndarray  # absolute [Pa]
    downstream_pressure: np.ndarray  # absolute [Pa]
    supply_temperature_K: np.ndarray
    tank_masses: np.ndarray  # supply-tank gas mass [kg]
    mass_start: float  # [kg]
    mass_end: float  # [kg]
    m_dot: float  # run-average blowdown mass flow [kg/s]
    gammas: np.ndarray  # cp/cv at supply conditions
    m_dot_ideals: np.ndarray  # ideal choked mass flow, Cd = 1 [kg/s]
    Cds: np.ndarray  # m_dot / m_dot_ideals (blowdown_orifice_cd.py)
    average_Cd: float  # mean of Cds (blowdown_orifice_cd.py)
    mass_actual_cumulative: np.ndarray  # measured cumulative mass loss [kg]
    mass_ideal_cumulative: np.ndarray  # integrated ideal mass loss, Cd = 1 [kg]
    Cd: float  # least-squares calibrated Cd (blowdown_orifice_cd_chat.py)
    Cd_integrated: float  # total measured / total ideal mass loss
    mass_actual_total: float  # [kg]
    mass_ideal_total: float  # [kg]
    mass_predicted_cumulative: np.ndarray  # Cd * mass_ideal_cumulative [kg]
    m_dot_predicted: np.ndarray  # Cd * m_dot_ideals [kg/s]

    def summary(self) -> str:
        """Terminal output of ``blowdown_orifice_cd_chat.py`` plus the average Cd
        that ``blowdown_orifice_cd.py`` showed in its plot legend."""
        return (
            f"Calibrated Cd: {self.Cd:.4f}\n"
            f"Integrated Cd: {self.Cd_integrated:.4f}\n"
            f"Actual Mass Loss: {self.mass_actual_total:.6f} kg\n"
            f"Ideal Mass Loss (Cd = 1): {self.mass_ideal_total:.6f} kg\n"
            f"Average Cd (run-average mdot / ideal mdot): {self.average_Cd:.4f}\n"
        )


def blowdown_orifice_cd(
    data: pd.DataFrame | str | Path = BLOWDOWN_INPUT_CSV,
    gas: str = DEFAULT_GAS,
    tank_volume_L: float = DEFAULT_TANK_VOLUME_L,
    orifice_diameter_in: float = DEFAULT_ORIFICE_DIAMETER_IN,
) -> BlowdownCdResult:
    """Orifice discharge coefficient from a rigid-tank blowdown.

    Parameters
    ----------
    data : blowdown table or CSV path (gauge pressures in psi, temperatures in K,
        time in ms; see :data:`REQUIRED_COLUMNS`).
    gas : CoolProp fluid name (e.g. ``"Air"``, ``"Oxygen"``, ``"Methane"``).
    tank_volume_L : rigid supply-tank internal volume [L].
    orifice_diameter_in : orifice diameter [in].
    """
    if not isinstance(data, pd.DataFrame):
        data = read_blowdown_csv(data)
    else:
        missing = [c for c in REQUIRED_COLUMNS if c not in data.columns]
        if missing:
            raise ValueError(f"Blowdown data is missing required columns: {missing}")

    if not (np.isfinite(tank_volume_L) and tank_volume_L > 0):
        raise ValueError(f"tank_volume_L must be finite and > 0, got {tank_volume_L!r}.")
    if not (np.isfinite(orifice_diameter_in) and orifice_diameter_in > 0):
        raise ValueError(f"orifice_diameter_in must be finite and > 0, got {orifice_diameter_in!r}.")

    # ---- data pre-processing ------------------------------------------------
    time_ms = _numeric_column(data, "time_ms")
    supply_pressure_psi = _numeric_column(data, "supply_pressure_psi")
    upstream_pressure_psi = _numeric_column(data, "upstream_pressure_psi")
    downstream_pressure_psi = _numeric_column(data, "downstream_pressure_psi")
    supply_temperature_K = _numeric_column(data, "supply_temperature_K")

    supply_pressure = (supply_pressure_psi + P_ATM_PSI) * PSI_TO_PA
    upstream_pressure = (upstream_pressure_psi + P_ATM_PSI) * PSI_TO_PA
    downstream_pressure = (downstream_pressure_psi + P_ATM_PSI) * PSI_TO_PA

    # The original scripts evaluated CoolProp one sample at a time, so an invalid
    # sample raised.  Vectorized PropsSI returns inf instead, so every input the
    # Cd calculation uses is checked here and the offending rows are reported.
    if len(time_ms) < 2:
        raise ValueError("Blowdown data needs at least two samples.")
    _require(np.isfinite(time_ms), "time_ms must be finite")
    if not time_ms[-1] > time_ms[0]:
        raise ValueError("time_ms must end later than it starts (run duration must be > 0).")
    _require(np.isfinite(supply_pressure) & (supply_pressure > 0),
             "supply_pressure_psi must be finite and above -14.6959 psig (absolute pressure > 0)")
    _require(np.isfinite(supply_temperature_K) & (supply_temperature_K > 0),
             "supply_temperature_K must be finite and > 0")

    time_s = time_ms / 1000

    # ---- gas properties -----------------------------------------------------
    R_u = cp.PropsSI("GAS_CONSTANT", gas)  # J/(mol*K)
    M = cp.PropsSI("MOLAR_MASS", gas)  # kg/mol
    R = R_u / M  # J/(kg*K)

    tank_volume = tank_volume_L * L_TO_M3
    orifice_radius = (orifice_diameter_in * IN_TO_M) / 2
    orifice_area = np.pi * orifice_radius**2

    # ---- tank mass and run-average mass flow -------------------------------
    density = _coolprop("Dmass", supply_pressure, supply_temperature_K, gas)
    tank_masses = density * tank_volume
    mass_start = float(tank_masses[0])
    mass_end = float(tank_masses[-1])
    m_dot = (mass_start - mass_end) / ((time_ms[-1] - time_ms[0]) / 1000)

    # ---- ideal choked mass flow (Cd = 1) at supply conditions --------------
    c_p = _coolprop("Cpmass", supply_pressure, supply_temperature_K, gas)
    c_v = _coolprop("Cvmass", supply_pressure, supply_temperature_K, gas)
    gammas = c_p / c_v
    _require(gammas > 1, "CoolProp gamma = cp/cv must be > 1")
    # = choked_flow.choked_mass_flow(1, A, P0, T0, gamma, R); written out to keep
    # the original operation order.
    m_dot_ideals = (
        orifice_area
        * supply_pressure
        / np.sqrt(supply_temperature_K)
        * np.sqrt(gammas / R)
        * (2 / (gammas + 1)) ** ((gammas + 1) / (2 * (gammas - 1)))
    )

    # ---- method 1: run-average mdot / instantaneous ideal mdot -------------
    Cds = m_dot / m_dot_ideals
    average_Cd = sum(Cds.tolist()) / len(Cds)

    # ---- method 2: least-squares fit of cumulative mass loss ---------------
    mass_actual_cumulative = tank_masses[0] - tank_masses

    mass_ideal_cumulative = np.zeros(len(time_s))
    for i in range(1, len(time_s)):
        dt = time_s[i] - time_s[i - 1]
        mass_ideal_cumulative[i] = (
            mass_ideal_cumulative[i - 1] + 0.5 * (m_dot_ideals[i - 1] + m_dot_ideals[i]) * dt
        )

    x = mass_ideal_cumulative[1:]
    y = mass_actual_cumulative[1:]
    Cd = float(np.sum(x * y) / np.sum(x**2))

    mass_actual_total = float(mass_actual_cumulative[-1])
    mass_ideal_total = float(mass_ideal_cumulative[-1])
    Cd_integrated = mass_actual_total / mass_ideal_total

    return BlowdownCdResult(
        gas=gas,
        tank_volume=tank_volume,
        orifice_diameter=orifice_diameter_in,
        orifice_area=float(orifice_area),
        R=R,
        time_ms=time_ms,
        time_s=time_s,
        supply_pressure=supply_pressure,
        upstream_pressure=upstream_pressure,
        downstream_pressure=downstream_pressure,
        supply_temperature_K=supply_temperature_K,
        tank_masses=tank_masses,
        mass_start=mass_start,
        mass_end=mass_end,
        m_dot=float(m_dot),
        gammas=gammas,
        m_dot_ideals=m_dot_ideals,
        Cds=Cds,
        average_Cd=float(average_Cd),
        mass_actual_cumulative=mass_actual_cumulative,
        mass_ideal_cumulative=mass_ideal_cumulative,
        Cd=Cd,
        Cd_integrated=float(Cd_integrated),
        mass_actual_total=mass_actual_total,
        mass_ideal_total=mass_ideal_total,
        mass_predicted_cumulative=Cd * mass_ideal_cumulative,
        m_dot_predicted=Cd * m_dot_ideals,
    )


# ===========================================================================
# Plots
# ===========================================================================


def plot_cd_vs_time(result: BlowdownCdResult):
    """Discharge coefficient vs. time with its average (blowdown_orifice_cd.py)."""
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(10, 6))
    plt.plot(result.time_ms, result.Cds, linewidth=2, label="$C_d$")
    plt.axhline(
        result.average_Cd,
        linestyle="--",
        linewidth=2,
        label=f"Average $C_d$ = {result.average_Cd:.3f}",
    )
    plt.xlabel("Time (ms)", fontsize=16)
    plt.ylabel("Discharge Coefficient, $C_d$", fontsize=16)
    plt.title("Discharge Coefficient vs. Time", fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=13)
    plt.tight_layout()
    return fig


def plot_mass_calibration(result: BlowdownCdResult):
    """Measured vs. calibrated-model cumulative mass loss (blowdown_orifice_cd_chat.py)."""
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(10, 6))
    plt.plot(result.time_ms, result.mass_actual_cumulative, linewidth=2, label="Measured Mass Loss")
    plt.plot(
        result.time_ms,
        result.mass_predicted_cumulative,
        linestyle="--",
        linewidth=2,
        label=f"Model, $C_d$ = {result.Cd:.3f}",
    )
    plt.xlabel("Time (ms)", fontsize=16)
    plt.ylabel("Cumulative Mass Loss (kg)", fontsize=16)
    plt.title("Orifice Blowdown Calibration", fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=13)
    plt.tight_layout()
    return fig


def main(argv: list[str] | None = None) -> BlowdownCdResult:
    parser = argparse.ArgumentParser(description="Orifice Cd from a supply-tank blowdown test.")
    parser.add_argument("--input", type=Path, default=BLOWDOWN_INPUT_CSV, help="blowdown CSV")
    parser.add_argument("--gas", default=DEFAULT_GAS, help="CoolProp fluid name (default Air)")
    parser.add_argument("--tank-volume-L", type=float, default=DEFAULT_TANK_VOLUME_L)
    parser.add_argument("--orifice-diameter-in", type=float, default=DEFAULT_ORIFICE_DIAMETER_IN)
    parser.add_argument("--no-show", action="store_true", help="do not open plot windows")
    args = parser.parse_args(argv)

    try:
        result = blowdown_orifice_cd(args.input, args.gas, args.tank_volume_L, args.orifice_diameter_in)
    except ValueError as exc:
        parser.exit(1, f"{parser.prog}: error: {exc}\n")
    print(result.summary(), end="")

    plot_cd_vs_time(result)
    plot_mass_calibration(result)
    if not args.no_show:
        import matplotlib.pyplot as plt

        plt.show()
    return result


if __name__ == "__main__":
    main()
