"""Unit conversions shared across the RDE projects.

``convert`` is the Python port of ``HADES_size_convert.m``; it keeps the
same lookup tables, the same aliases and the same order of checks so that
results are identical to the MATLAB version (with the one documented fix
for ``'atm'``, see ``HADES/sizing_model/PORTING_NOTES.md``).

The module-level constants are the conversion factors that were hard-coded
inline throughout the MATLAB sizing scripts. They are kept with their
original (sometimes rounded) values so that ported scripts reproduce the
MATLAB numbers exactly; do not "improve" them without updating the tests.
"""

from __future__ import annotations

import math

import numpy as np

__all__ = [
    "convert",
    "LBM_TO_KG",
    "LBM_PER_KG",
    "PSI_TO_PA",
    "PSI_TO_PA_ROUNDED",
    "PA_PER_BAR",
    "PSI_PER_BAR",
    "MPA_TO_PSI",
    "KSI_TO_MPA",
    "N_TO_LBF",
    "MM_PER_IN",
    "M_PER_IN",
    "FT_PER_M",
    "ATM_TO_PA",
    "MMHG_TO_PA",
    "CEA2_PSIA_TO_BAR",
    "R_UNIVERSAL",
    "STEFAN_BOLTZMANN",
    "convert_pressure",
    "convert_temperature",
    "psia_to_bar",
    "bar_to_psia",
    "c_to_k",
    "k_to_c",
    "c_to_f",
    "f_to_c",
]

# --- Factors used inline in the MATLAB sources (kept verbatim) -------------
LBM_TO_KG = 0.453592  # kg per lbm (HADES_size_thrust/h2Inj/airInj, feed sizing)
LBM_PER_KG = 2.205  # lbm per kg (HADES_size_convert 'lbms', main, sonic nozzle)
PSI_TO_PA = 6894.76  # Pa per psi (CEA wrappers, feed sizing V1-V3; V5 and blowdown keep their own)
PSI_TO_PA_ROUNDED = 6895.0  # Pa per psi (HADES_size_h2Inj)
PA_PER_BAR = 1.0e5  # (DEIMOS)
PSI_PER_BAR = 14.5038  # HADES_size_convert / HoopStressTemps
MPA_TO_PSI = 145.03773773  # HoopStressTemps
KSI_TO_MPA = 6.89476  # HoopStressTemps
N_TO_LBF = 0.2248  # HADES_size_thrust
MM_PER_IN = 25.4
M_PER_IN = 0.0254  # (feed sizing, DEIMOS CdCalc)
FT_PER_M = 3.28084
ATM_TO_PA = 101325.0
MMHG_TO_PA = 133.322
#: psia -> bar factor applied internally by NASA CEA2 when it reads
#: ``p,psia=`` in an input deck (1/14.696006 atm per psi, 1.01325 bar/atm).
CEA2_PSIA_TO_BAR = 1.01325 / 14.696006
R_UNIVERSAL = 8.314  # J/(mol K), value used in HADES_size_airInj
STEFAN_BOLTZMANN = 5.67e-8  # W/(m^2 K^4), value used in the cooldown scripts


# --- HADES_size_convert port -----------------------------------------------
# Pressure: units per bar (bar is the pivot).  NOTE: the MATLAB table had
# 'atm' = 1.01325 (that is bar-per-atm, i.e. inverted).  It is fixed here to
# atm-per-bar; see PORTING_NOTES.md.
_PRESSURE_PER_BAR = {"bar": 1.0, "atm": 1.0 / 1.01325, "psia": 14.5038, "mmhg": 750.062}
_SPEED_PER_MS = {"ms": 1.0, "fts": 3.28084}
_DIST_PER_M = {"m": 1.0, "ft": 3.28084, "in": 39.37}
_MDOT_PER_KGS = {"kgs": 1.0, "lbms": 2.205}

_TEMP_K = ("k", "kelvin")
_TEMP_C = ("c", "degc", "celsius")
_TEMP_F = ("f", "degf", "fahrenheit")
_TEMP_R = ("r", "rankine")


def _isnan(x) -> bool:
    return bool(np.any(np.isnan(x)))


def convert(val, from_unit: str, to_unit: str):
    """Convert pressure, speed, temperature, distance or mass-flow units.

    Port of ``HADES_size_convert(val, from_unit, to_unit)``. Unit names are
    case-insensitive.

    Supported units
    ---------------
    pressure  : ``bar``, ``atm``, ``psia``, ``mmHg``
    speed     : ``m/s`` (``ms``), ``ft/s`` (``fts``)
    temperature: ``K``/``kelvin``, ``C``/``degC``/``celsius``,
                 ``F``/``degF``/``fahrenheit``, ``R``/``rankine``
    distance  : ``m``, ``ft``, ``in`` (1 m = 39.37 in, as in the MATLAB)
    mass flow : ``kgs``, ``lbms`` (1 kg = 2.205 lbm, as in the MATLAB)

    ``val`` may be a scalar or a numpy array (element-wise).

    Raises
    ------
    ValueError
        If the unit pair is not supported (MATLAB ``error``).
    """
    f = from_unit.lower()
    t = to_unit.lower()

    # Pressure (bar pivot)
    if f in _PRESSURE_PER_BAR and t in _PRESSURE_PER_BAR:
        return val / _PRESSURE_PER_BAR[f] * _PRESSURE_PER_BAR[t]

    # Speed (aliases first, exactly as in the MATLAB)
    if f == "m/s":
        f = "ms"
    if t == "m/s":
        t = "ms"
    if f == "ft/s":
        f = "fts"
    if t == "ft/s":
        t = "fts"
    if f in _SPEED_PER_MS and t in _SPEED_PER_MS:
        return val / _SPEED_PER_MS[f] * _SPEED_PER_MS[t]

    # Temperature (via kelvin)
    if f in _TEMP_K:
        temp_k = val
    elif f in _TEMP_C:
        temp_k = val + 273.15
    elif f in _TEMP_F:
        temp_k = (val - 32) * 5 / 9 + 273.15
    elif f in _TEMP_R:
        temp_k = val * 5 / 9
    else:
        temp_k = math.nan
    if not _isnan(temp_k):
        if t in _TEMP_K:
            return temp_k
        if t in _TEMP_C:
            return temp_k - 273.15
        if t in _TEMP_F:
            return (temp_k - 273.15) * 9 / 5 + 32
        if t in _TEMP_R:
            return temp_k * 9 / 5
        raise ValueError("Target unit not recognized.")

    # Distance
    if f in _DIST_PER_M and t in _DIST_PER_M:
        return val / _DIST_PER_M[f] * _DIST_PER_M[t]

    # Mass flow rate
    if f in _MDOT_PER_KGS and t in _MDOT_PER_KGS:
        return val / _MDOT_PER_KGS[f] * _MDOT_PER_KGS[t]

    raise ValueError(f"Unit conversion from {f} to {t} is not supported.")


# --- Local helpers of the CEA wrappers (HADES_size_ceaDet/ceaRocket) --------
_PA_PER_PUNIT = {"psia": PSI_TO_PA, "atm": ATM_TO_PA, "bar": PA_PER_BAR, "mmhg": MMHG_TO_PA}


def convert_pressure(p_in, unit_in: str, unit_out: str):
    """Pressure conversion used inside the CEA wrappers (``convertPressure``).

    Units: ``psia`` (6894.76 Pa), ``atm`` (101325 Pa), ``bar`` (1e5 Pa),
    ``mmHg`` (133.322 Pa). Case-insensitive.
    """
    ui, uo = unit_in.lower(), unit_out.lower()
    if ui not in _PA_PER_PUNIT or uo not in _PA_PER_PUNIT:
        raise ValueError("Unsupported pressure unit")
    return p_in * _PA_PER_PUNIT[ui] / _PA_PER_PUNIT[uo]


def convert_temperature(T_in, unit_in: str, unit_out: str):
    """Temperature conversion used inside the CEA wrappers (``convertTemperature``).

    Units: ``K``, ``C``, ``F``, ``R``. Case-insensitive.
    """
    ui, uo = unit_in.upper(), unit_out.upper()
    if ui == "K":
        T_K = T_in
    elif ui == "C":
        T_K = T_in + 273.15
    elif ui == "F":
        T_K = (T_in - 32) * 5 / 9 + 273.15
    elif ui == "R":
        T_K = T_in * 5 / 9
    else:
        raise ValueError("Unsupported temperature unit")
    if uo == "K":
        return T_K
    if uo == "C":
        return T_K - 273.15
    if uo == "F":
        return (T_K - 273.15) * 9 / 5 + 32
    if uo == "R":
        return T_K * 9 / 5
    raise ValueError("Unsupported temperature unit")


def psia_to_bar(p_psia):
    """psia -> bar using the HADES 14.5038 psi/bar factor."""
    return p_psia / PSI_PER_BAR


def bar_to_psia(p_bar):
    """bar -> psia using the HADES 14.5038 psi/bar factor."""
    return p_bar * PSI_PER_BAR


def c_to_k(T_C):
    """Celsius -> kelvin."""
    return T_C + 273.15


def k_to_c(T_K):
    """Kelvin -> Celsius."""
    return T_K - 273.15


def c_to_f(T_C):
    """Celsius -> Fahrenheit."""
    return T_C * 9 / 5 + 32


def f_to_c(T_F):
    """Fahrenheit -> Celsius."""
    return (T_F - 32) * 5 / 9
