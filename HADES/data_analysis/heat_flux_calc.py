"""Two-thermocouple surface temperature and heat flux estimate (port of ``HeatFluxCalc.m``).

Two thermocouples sit at depths ``x1 < x2`` (metres) below the heated surface.
Both signals are smoothed with a Savitzky-Golay filter and differentiated, then
combined to estimate the surface temperature ``T0`` and the surface heat flux
``q0`` (positive into the wall).

There is NO default model. The caller must choose ``equations``:

* ``"cubic"``: derived from the 1D heat equation ``dT/dt = alpha * d2T/dx2``,
  with ``dT/dt`` linear in depth (a cubic temperature profile). Exact for any
  temperature field that is cubic in depth, e.g. uniform heating
  ``T = 300 - 20x + 50x^2 + 100*alpha*t`` (q0 = +7800 W/m^2).
* ``"legacy"``: Eq. 3 / Eq. 4 exactly as written in ``HeatFluxCalc.m``.
  KNOWN INVALID: the transient terms ``alpha*(...)*T_dot`` have the wrong units
  (K m^3/s^2 instead of K/m) and the wrong sign. On the uniform-heating solution
  above they give q0 = -8287.5 W/m^2. Kept only to reproduce the original
  numbers. Using it emits :class:`LegacyEquationWarning`.

Sensor depths are explicit inputs in metres (``depths_m``). The MATLAB file gave
``x1 = 0.25``, ``x2 = 0.575`` with no unit. That pair is kept as
:data:`MATLAB_SENSOR_DEPTHS_M` (read as metres, as the MATLAB formulas did),
but the unit is unconfirmed and is probably mm. Passing these values emits
:class:`SensorDepthUnitsWarning`.

Other porting changes (numbers in ``PORTING_NOTES.md``):

* BUG FIX: MATLAB computed the "temperature derivative" as ``T_s / dt``. The
  port uses the Savitzky-Golay first derivative with the same order and window
  (``savgol_filter(..., deriv=1, delta=dt)``).
* The hard-coded ``/Users/Abby/Desktop/RDE/...`` data path now points to the
  copy of ``sample_thermocouple_data.csv`` next to this module.

Run as a script (both choices are required)::

    uv run python -m HADES.data_analysis.heat_flux_calc --equations cubic --depths-m 0.25e-3 0.575e-3 [data.csv]
    uv run python -m HADES.data_analysis.heat_flux_calc --equations legacy --matlab-depths
"""

from __future__ import annotations

import argparse
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

# ---------------- CONSTANTS ----------------
#: Sensor depths from HeatFluxCalc.m (x1 = 0.25, x2 = 0.575), read as metres as
#: the MATLAB formulas did. UNITS UNCONFIRMED (probably mm); using them warns.
MATLAB_SENSOR_DEPTHS_M = (0.25, 0.575)
K = 390.0  # thermal conductivity (W/m-K)
ALPHA = 1.11e-4  # thermal diffusivity (m^2/s)

SGOLAY_ORDER = 2
SGOLAY_WINDOW = 21

EQUATIONS = ("legacy", "cubic")

DATA_FILE = Path(__file__).resolve().parent / "sample_thermocouple_data.csv"


class LegacyEquationWarning(UserWarning):
    """The original (dimensionally inconsistent) Eq. 3/4 from HeatFluxCalc.m were used."""


class SensorDepthUnitsWarning(UserWarning):
    """The unconfirmed MATLAB sensor depths (0.25 / 0.575 "m") were used."""


@dataclass
class HeatFluxResult:
    t: np.ndarray
    T1: np.ndarray
    T2: np.ndarray
    dt: float
    T1_s: np.ndarray
    T2_s: np.ndarray
    T1_dot: np.ndarray
    T2_dot: np.ndarray
    T0: np.ndarray
    q0: np.ndarray
    equations: str
    depths_m: tuple[float, float]


def readmatrix(path) -> np.ndarray:
    """Minimal ``readmatrix`` for a numeric CSV file with an optional header line.

    Values are parsed exactly (``float_precision="round_trip"``), as MATLAB does.
    """
    with open(path, encoding="utf-8") as fh:
        first = fh.readline()
    try:
        [float(v) for v in first.strip().split(",")]
        skip = 0
    except ValueError:  # text header (e.g. "time_s,T1_K,T2_K")
        skip = 1
    df = pd.read_csv(path, header=None, skiprows=skip, float_precision="round_trip")
    return df.to_numpy(dtype=float)


def sgolayfilt(x, order: int, framelen: int) -> np.ndarray:
    """MATLAB ``sgolayfilt(x, order, framelen)`` for a vector.

    scipy's ``mode='interp'`` handles the edges the way MATLAB does: the first
    and last ``(framelen-1)/2`` outputs come from the polynomial fitted to the
    first and last frame.
    """
    return savgol_filter(np.asarray(x, dtype=float), framelen, order, mode="interp")


def sgolay_derivative(x, order: int, framelen: int, dt: float) -> np.ndarray:
    """First time derivative of the same Savitzky-Golay fit (``d/dt`` of the local polynomial)."""
    return savgol_filter(np.asarray(x, dtype=float), framelen, order, deriv=1, delta=dt, mode="interp")


# ---------------------------------------------------------------- legacy (invalid)
def _warn_legacy():
    warnings.warn(
        "Using the original HeatFluxCalc.m Eq. 3/4, which are dimensionally inconsistent "
        "(alpha*(...)*dT/dt terms) and fail exact heat-equation solutions; results are "
        "known to be physically invalid. Use equations='cubic'.",
        LegacyEquationWarning,
        stacklevel=3,
    )


def surface_temperature_legacy(T1_s, T2_s, T1_dot, T2_dot, x1_m, x2_m, alpha=ALPHA):
    """Eq. 3 exactly as written in ``HeatFluxCalc.m``. KNOWN INVALID (see module docstring)."""
    _warn_legacy()
    x1, x2 = x1_m, x2_m
    return (
        (x2 / (x2 - x1)) * T1_s
        - (x1 / (x2 - x1)) * T2_s
        + alpha * ((x2 * (x2 - 2 * x1)) / (6 * (x2 - x1))) * T1_dot
        + alpha * ((x1 * (2 * x2 - x1)) / (6 * (x2 - x1))) * T2_dot
    )


def surface_heat_flux_legacy(T1_s, T2_s, T1_dot, T2_dot, x1_m, x2_m, k=K, alpha=ALPHA):
    """Eq. 4 exactly as written in ``HeatFluxCalc.m`` (W/m^2). KNOWN INVALID (see module docstring)."""
    _warn_legacy()
    x1, x2 = x1_m, x2_m
    return -k * (
        (T2_s - T1_s) / (x2 - x1)
        + alpha * ((x2**2 + x1 * x2 - 2 * x1**2) / (6 * (x2 - x1))) * T1_dot
        + alpha * ((2 * x2**2 - x1 * x2 - x1**2) / (6 * (x2 - x1))) * T2_dot
    )


# ---------------------------------------------------------------- cubic (derived)
def surface_temperature_cubic(T1_s, T2_s, T1_dot, T2_dot, x1_m, x2_m, alpha=ALPHA):
    """Surface temperature from 1D conduction with ``dT/dt`` linear in depth.

    ``T0 = (x2*T1 - x1*T2)/d + x1*x2*((2*x2 - x1)*T1_dot + (x2 - 2*x1)*T2_dot) / (6*alpha*d)``,
    where ``d = x2 - x1``. Exact for any temperature field that is cubic in depth.
    """
    x1, x2 = x1_m, x2_m
    d = x2 - x1
    return (x2 * T1_s - x1 * T2_s) / d + x1 * x2 * ((2 * x2 - x1) * T1_dot + (x2 - 2 * x1) * T2_dot) / (6 * alpha * d)


def surface_heat_flux_cubic(T1_s, T2_s, T1_dot, T2_dot, x1_m, x2_m, k=K, alpha=ALPHA):
    """Surface heat flux [W/m^2] into the wall, from the same cubic-profile model.

    ``q0 = -k*(T2 - T1)/d + k*((2*x2^2 + 2*x1*x2 - x1^2)*T1_dot + (x2^2 - 2*x1*x2 - 2*x1^2)*T2_dot) / (6*alpha*d)``.
    """
    x1, x2 = x1_m, x2_m
    d = x2 - x1
    return -k * (T2_s - T1_s) / d + k * (
        (2 * x2**2 + 2 * x1 * x2 - x1**2) * T1_dot + (x2**2 - 2 * x1 * x2 - 2 * x1**2) * T2_dot
    ) / (6 * alpha * d)


def _check_depths(depths_m) -> tuple[float, float]:
    x1, x2 = (float(v) for v in depths_m)
    if not 0 <= x1 < x2:
        raise ValueError("depths_m must be (x1, x2) in metres with 0 <= x1 < x2")
    if (x1, x2) == MATLAB_SENSOR_DEPTHS_M:
        warnings.warn(
            "Sensor depths 0.25 / 0.575 come from HeatFluxCalc.m without units and are used as "
            "metres; the units are unconfirmed (probably mm, i.e. depths_m=(0.25e-3, 0.575e-3)).",
            SensorDepthUnitsWarning,
            stacklevel=3,
        )
    return x1, x2


def HeatFluxCalc(
    data_file: str | Path = DATA_FILE,
    *,
    equations: str,
    depths_m: tuple[float, float],
    k: float = K,
    alpha: float = ALPHA,
    sgolay_order: int = SGOLAY_ORDER,
    sgolay_window: int = SGOLAY_WINDOW,
    plot: bool = True,
) -> HeatFluxResult:
    """Read ``[t, T1, T2]`` from ``data_file``, compute T0 and q0, and optionally plot.

    ``equations`` ("cubic" or "legacy") and ``depths_m`` (x1, x2 in metres) are
    required keywords. See the module docstring.
    """
    if equations not in EQUATIONS:
        raise ValueError(f"equations must be one of {EQUATIONS}, got {equations!r}")
    x1, x2 = _check_depths(depths_m)

    # ---------------- DATA READING ----------------
    data = readmatrix(data_file)
    t = data[:, 0]
    T1 = data[:, 1]
    T2 = data[:, 2]
    dt = float(np.mean(np.diff(t)))

    # ---------------- SAVITZKY-GOLAY FILTER ----------------
    T1_s = sgolayfilt(T1, sgolay_order, sgolay_window)
    T2_s = sgolayfilt(T2, sgolay_order, sgolay_window)

    # ---------------- TEMPERATURE DERIVATIVES (fixed: was T_s / dt) ----------------
    T1_dot = sgolay_derivative(T1, sgolay_order, sgolay_window, dt)
    T2_dot = sgolay_derivative(T2, sgolay_order, sgolay_window, dt)

    # ---------------- SURFACE TEMPERATURE AND HEAT FLUX ----------------
    if equations == "cubic":
        T0 = surface_temperature_cubic(T1_s, T2_s, T1_dot, T2_dot, x1, x2, alpha)
        q0 = surface_heat_flux_cubic(T1_s, T2_s, T1_dot, T2_dot, x1, x2, k, alpha)
    else:
        T0 = surface_temperature_legacy(T1_s, T2_s, T1_dot, T2_dot, x1, x2, alpha)
        q0 = surface_heat_flux_legacy(T1_s, T2_s, T1_dot, T2_dot, x1, x2, k, alpha)

    if plot:
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots()
        ax.plot(t, T1, "b", t, T2, "r")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Temperature (K)")
        ax.legend(["T1", "T2"])
        ax.set_title("Measured Thermocouple Temperatures")

        fig, ax = plt.subplots()
        ax.plot(t, T0, "c")
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Surface Temperature (K)")
        ax.set_title(f"Calculated Surface Temperature ({equations})")

        fig, ax = plt.subplots()
        ax.plot(t, q0)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Heat Flux (W/m^2)")
        ax.set_title(f"Calculated Surface Heat Flux ({equations})")

    return HeatFluxResult(t, T1, T2, dt, T1_s, T2_s, T1_dot, T2_dot, T0, q0, equations, (x1, x2))


def _parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="python -m HADES.data_analysis.heat_flux_calc", description=__doc__.split("\n")[0])
    p.add_argument("data_file", nargs="?", default=str(DATA_FILE), help="CSV with columns t, T1, T2")
    p.add_argument("--equations", required=True, choices=EQUATIONS, help="'cubic' (derived) or 'legacy' (original, invalid)")
    depth = p.add_mutually_exclusive_group(required=True)
    depth.add_argument("--depths-m", nargs=2, type=float, metavar=("X1_M", "X2_M"), help="sensor depths in metres")
    depth.add_argument("--matlab-depths", action="store_true", help="use the unconfirmed MATLAB depths 0.25/0.575 (as metres)")
    return p


def main(argv: list[str] | None = None) -> None:
    import matplotlib.pyplot as plt

    args = _parser().parse_args(argv)
    depths = MATLAB_SENSOR_DEPTHS_M if args.matlab_depths else tuple(args.depths_m)
    HeatFluxCalc(args.data_file, equations=args.equations, depths_m=depths)
    plt.show()


if __name__ == "__main__":
    main()
