"""HADES feed-system sonic nozzle / sonic orifice sizing.

Python port of ``HADES/FEED SIZING/SonicNozzleSizingV1.m`` ... ``V5.m`` and the
earlier ``SonicNozzleSizing.py`` prototype (programmers: Deepesh Balwani,
Noah Ha).  Each original is a different sizing workflow, so each one is a
function here (CLI subcommand in the first column); they share the ideal-gas
choked-flow relations in :mod:`HADES.feed_sizing.choked_flow`.

=========  ==============================================  ==========================
Original   Workflow                                        Python function
=========  ==============================================  ==========================
feedline   prototype ``SonicNozzleSizing.py``: throat      :func:`size_sonic_nozzle_with_feed_line`
           from feed-line Mach number (area-Mach bug
           corrected, see PORTING_NOTES.md)
V1         single-point sonic nozzle for a fixed           :func:`size_single_sonic_nozzle`
           downstream pressure (air)
V2         one fixed sonic orifice per gas, fixed          :func:`size_fixed_sonic_orifice`
           downstream pressure per gas (CSV in/out)
V3         two-stage: size injector area from chamber      :func:`size_two_stage_choked`
           pressure, then the upstream sonic orifice
V4         manifold pressure from a known injector area    :func:`manifold_pressure_from_injector`
V5         known injector area -> manifold pressure ->     :func:`size_injector_sonic_orifice`
           fixed sonic orifice -> supply pressure          (canonical, CLI default)
=========  ==============================================  ==========================

Inputs are CSV files with columns ``H2, Air, CH4, GOx`` holding required mass
flows in kg/s (blank cells allowed).  Result CSVs keep the MATLAB
``writetable`` layout: the input columns followed by the result columns,
logicals written as ``1``/``0``, missing values written as ``NaN`` and numbers
written with 15 significant digits.

The bundled default inputs are resolved relative to this module; generated
result CSVs are written to the CURRENT WORKING DIRECTORY unless ``--output``
is given (so the bundled reference results are never overwritten).

Run from the repository root::

    uv run python -m HADES.feed_sizing.sonic_nozzle_sizing            # V5 -> ./flow_requirements_manifold_input_combined_results.csv
    uv run python -m HADES.feed_sizing.sonic_nozzle_sizing v2         # any version
    uv run python -m HADES.feed_sizing.sonic_nozzle_sizing feedline   # prototype, prints only
    uv run python -m HADES.feed_sizing.sonic_nozzle_sizing v3 --input my.csv --output out.csv
"""

from __future__ import annotations

import argparse
import re
import warnings
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from HADES.feed_sizing.choked_flow import (
    area_to_diameter,
    choked_area,
    choked_flow_parameter,
    choked_upstream_pressure,
    critical_pressure_ratio,
)

MODULE_DIR = Path(__file__).resolve().parent

#: Default CSV inputs (resolved relative to this module, not the CWD).
FLOW_REQUIREMENTS_CSV = MODULE_DIR / "flow_requirements.csv"
MANIFOLD_INPUT_CSV = MODULE_DIR / "flow_requirements_manifold_input.csv"

# ---------------------------------------------------------------------------
# Conversions (each MATLAB version carried its own psi->Pa constant; they are
# kept separate so outputs match the originals).
# ---------------------------------------------------------------------------

IN_TO_M = 0.0254
LB_TO_KG = 0.453592
PSI_TO_PA = 6894.76  # V1, V2, V3
PSI_TO_PA_V5 = 6894.757293  # V5

GAS_NAMES = ("H2", "Air", "CH4", "GOx")


@dataclass(frozen=True)
class GasProperties:
    """Ideal-gas properties of one feed gas (identical in V2-V5)."""

    R: float  # specific gas constant [J/(kg*K)]
    gamma: float  # specific heat ratio [-]
    T0: float  # upstream stagnation temperature [K]


#: Approximate gas properties near room temperature (SonicNozzleSizingV2-V5).
GAS_PROPERTIES: dict[str, GasProperties] = {
    "H2": GasProperties(R=4124.0, gamma=1.405, T0=283.0),
    "Air": GasProperties(R=287.05, gamma=1.400, T0=283.0),
    "CH4": GasProperties(R=518.3, gamma=1.310, T0=283.0),
    "GOx": GasProperties(R=259.84, gamma=1.400, T0=283.0),
}

# --- V2 defaults -----------------------------------------------------------
V2_CD = {"H2": 0.85, "Air": 0.85, "CH4": 0.85, "GOx": 0.85}
V2_P_DOWN_BAR = {"H2": 20.0, "Air": 20.0, "CH4": 24.0, "GOx": 36.0}  # [bar abs]
V2_P_MARGIN_PSI = {"H2": 20.0, "Air": 20.0, "CH4": 20.0, "GOx": 20.0}  # [psi]

# --- V3 defaults -----------------------------------------------------------
V3_CD_INJ = {"H2": 0.85, "Air": 0.85, "CH4": 0.85, "GOx": 0.85}
V3_CD_ORIFICE = {"H2": 0.85, "Air": 0.85, "CH4": 0.85, "GOx": 0.85}

# --- V4 / V5 defaults ------------------------------------------------------
#: Total physical injector flow area per gas [mm^2].
INJECTOR_AREA_MM2 = {"H2": 25.13, "Air": 304.77, "CH4": 56.55, "GOx": 86.59}
#: Effective injector discharge coefficients.
INJECTOR_CD = {"H2": 0.71, "Air": 0.62, "CH4": 0.60, "GOx": 0.60}
#: Sonic-orifice discharge coefficients (V5).
SONIC_ORIFICE_CD = {"H2": 0.85, "Air": 0.85, "CH4": 0.85, "GOx": 0.85}


# ===========================================================================
# CSV input / output helpers
# ===========================================================================


def read_flow_requirements(input_file: str | Path = FLOW_REQUIREMENTS_CSV) -> pd.DataFrame:
    """Read a flow-requirements CSV (MATLAB ``readtable(..., 'preserve')``).

    Raises ``ValueError`` if a required ``H2/Air/CH4/GOx`` column is missing
    (matched case-insensitively, like ``strcmpi``).
    """
    data = pd.read_csv(input_file, float_precision="round_trip")
    _gas_column_indices(data)
    return data


def _gas_column_indices(data: pd.DataFrame) -> dict[str, int]:
    names = [str(c) for c in data.columns]
    indices = {}
    for gas in GAS_NAMES:
        matches = [i for i, name in enumerate(names) if name.lower() == gas.lower()]
        if not matches:
            raise ValueError(f'CSV must contain a column named "{gas}".')
        indices[gas] = matches[0]
    return indices


def _mass_flows(data: pd.DataFrame) -> dict[str, np.ndarray]:
    """Extract the H2/Air/CH4/GOx mass-flow columns as float arrays [kg/s]."""
    indices = _gas_column_indices(data)
    return {gas: data.iloc[:, idx].to_numpy(dtype=float) for gas, idx in indices.items()}


def _as_table(data: pd.DataFrame | str | Path) -> pd.DataFrame:
    if isinstance(data, pd.DataFrame):
        _gas_column_indices(data)
        return data.copy()
    return read_flow_requirements(data)


def write_results_csv(table: pd.DataFrame, output_file: str | Path) -> Path:
    """Write a results table like MATLAB ``writetable`` (15 significant digits,
    ``NaN`` for missing values, logicals as ``1``/``0``)."""
    out = table.copy()
    for column in out.columns:
        if out[column].dtype == bool:
            out[column] = out[column].astype(int)
    output_file = Path(output_file)
    out.to_csv(output_file, index=False, float_format="%.15g", na_rep="NaN", lineterminator="\n")
    return output_file


def _merged(defaults: Mapping[str, float], override: Mapping[str, float] | None) -> dict[str, float]:
    merged = dict(defaults)
    if override:
        unknown = set(override) - set(GAS_NAMES)
        if unknown:
            raise ValueError(f"Unknown gas name(s): {sorted(unknown)}; expected {GAS_NAMES}.")
        merged.update(override)
    return merged


def _gas_props(gas_properties: Mapping[str, GasProperties] | None) -> dict[str, GasProperties]:
    return {**GAS_PROPERTIES, **(gas_properties or {})}


_NONFINITE = re.compile(r"(?<![A-Za-z])(nan|inf)(?![A-Za-z])")


def _fprintf(fmt: str, *args) -> str:
    """``sprintf`` with MATLAB's spelling of non-finite numbers (NaN, Inf)."""
    return _NONFINITE.sub(lambda m: "NaN" if m.group(1) == "nan" else "Inf", fmt % args)


def _valid(m_dot: np.ndarray) -> np.ndarray:
    """Positive, non-NaN mass-flow points (MATLAB ``(m_dot > 0) & ~isnan(m_dot)``)."""
    with np.errstate(invalid="ignore"):
        return (m_dot > 0) & ~np.isnan(m_dot)


# ===========================================================================
# Prototype - sonic nozzle sized from feed-line geometry (SonicNozzleSizing.py)
# ===========================================================================


@dataclass
class FeedLineNozzleResult:
    """Outputs of :func:`size_sonic_nozzle_with_feed_line` (SonicNozzleSizing.py variables, SI)."""

    m_dot: float  # required mass flow [kg/s]
    P_m: float  # manifold / feed-line static pressure downstream of the nozzle [Pa]
    R: float  # [J/(kg*K)]
    gamma: float  # [-]
    T: float  # feed-line static temperature [K]
    mass_flow_margin: float  # factor applied to m_dot in the feed-line Mach number [-]
    crit_ratio: float  # critical pressure ratio [-]
    P_t_u_req: float  # legacy threshold P_m / crit_ratio [Pa] (NOT a valid choking requirement)
    wall_thickness: float  # feed-line wall thickness [in]
    feed_diameter: float  # feed-line inner diameter [m]
    A_2: float  # feed-line flow area [m^2]
    M_2: float  # feed-line Mach number [-]
    A_factor: float  # M_2 * A_2/A_t = [(2/(g+1)) (1 + (g-1)/2 M_2^2)]^((g+1)/(2(g-1)))
    A_t: float  # throat area [m^2]
    D_t: float  # throat diameter [m]
    T_t: float  # total temperature [K]
    P_t_u: float  # upstream total pressure for choked m_dot through A_t [Pa]
    #: Legacy diagnostic ``P_t_u > P_m / crit_ratio`` from the prototype's
    #: "Flow is choked!" test.  NOT a valid choking criterion: it compares the
    #: feed-line (downstream) pressure with the THROAT critical pressure and so
    #: ignores the isentropic pressure recovery between throat and feed line that
    #: this same model assumes; see PORTING_NOTES.md.
    legacy_pressure_check_passed: bool

    LEGACY_CHECK_NOTE = (
        "NOTE: legacy diagnostic from SonicNozzleSizing.py, not a valid choking criterion "
        "(it ignores pressure recovery downstream of the throat); "
        "see HADES/feed_sizing/PORTING_NOTES.md."
    )

    def summary(self) -> str:
        """Terminal output of SonicNozzleSizing.py, with the old choke verdict relabelled
        as a legacy diagnostic."""
        outcome = "passed" if self.legacy_pressure_check_passed else "not passed"
        return (
            "\n"
            f"M_2 = {self.M_2}\n"
            f"A_2 = {self.A_2} m^2\n"
            f"A_t = {self.A_t} m^2\n"
            f"Throat Diameter: {self.D_t * 1000} mm\n"
            f"Upstream Pressure: {self.P_t_u * 10 ** (-5)} bar\n"
            f"Legacy threshold P_m / critical ratio: {self.P_t_u_req * 10 ** (-5)} bar\n"
            f"Legacy pressure-ratio check (Upstream Pressure > threshold): {outcome}\n"
            f"{self.LEGACY_CHECK_NOTE}\n"
            "\n"
        )


def size_sonic_nozzle_with_feed_line(
    m_dot: float = 3 * LB_TO_KG,
    P_m: float = 2 * 10**6,
    R: float = 287.0,
    gamma: float = 1.4,
    T: float = 283.0,
    mass_flow_margin: float = 0.85,
    wall_thickness_in: float = 0.095,
    feed_OD_in: float = 16 / 16,
) -> FeedLineNozzleResult:
    """Sonic nozzle sized from the downstream feed-line geometry (prototype ``SonicNozzleSizing.py``).

    Unlike V1-V5 this workflow makes the throat depend on the feed line:

    1. Feed-line inner diameter ``OD - 2*wall`` gives the flow area ``A_2``.
    2. Feed-line Mach number from continuity at the static manifold pressure
       ``P_m`` and static temperature ``T``:
       ``M_2 = m_dot*mass_flow_margin / (P_m*A_2) * sqrt(R*T/gamma)``.
    3. Throat = isentropic sonic area for that Mach number,
       ``A_t = A_2 / (A/A*)(M_2)`` with
       ``A/A* = (1/M) [(2/(g+1)) (1 + (g-1)/2 M^2)]^((g+1)/(2(g-1)))``.
    4. Total temperature ``T_t = T (1 + (g-1)/2 M_2^2)`` and the upstream total
       pressure that passes ``m_dot`` choked through ``A_t`` (no Cd, no margin).
    5. Legacy diagnostic ``legacy_pressure_check_passed = P_t_u > P_m / crit_ratio``
       (printed by the prototype as "Flow is choked!" / "Flow does not choke :(").

    Step 5 is kept only as a historical diagnostic and is NOT a choking
    determination: steps 3-4 already put the throat exactly at M = 1 with
    isentropic recovery to ``P_m`` downstream, so comparing ``P_m`` with the
    throat critical pressure ignores that recovery.  E.g. with
    ``mass_flow_margin=1`` the model is self-consistently sonic (``A_t = A*``,
    ``P_t_u`` = isentropic total pressure of the feed-line state = 23.613659 bar)
    yet the legacy check fails (threshold 37.858583 bar).  No replacement choking
    criterion is provided.

    PORTING FIX: the prototype had ``2/(gamma+1)`` outside the power in step 3
    (``A_factor = 2/(g+1) * (1 + ...)**(...)``).  The corrected form here is the one
    the authors used in the MATLAB translation (commit f8c4cd3).  On the default
    inputs D_t goes 13.852781 -> 16.623337 mm and P_t_u 38.232314 -> 26.550218 bar;
    see PORTING_NOTES.md.

    The other quirks of the prototype are kept as written: ``mass_flow_margin``
    (0.85) only scales the flow used for ``M_2`` while ``P_t_u`` uses the unscaled
    ``m_dot`` (so by default ``P_t_u`` is 1/0.85 times the isentropic total pressure
    of the feed-line state), and ``R = 287``.  A ``RuntimeWarning`` is issued if
    ``M_2 >= 1`` (the subsonic feed-line relation no longer applies).
    """
    crit_ratio = critical_pressure_ratio(gamma)
    P_t_u_req = P_m / crit_ratio

    feed_diameter = (feed_OD_in - wall_thickness_in * 2) * IN_TO_M
    A_2 = np.pi * feed_diameter**2 / 4

    M_2 = (m_dot * mass_flow_margin) / (P_m * A_2) * np.sqrt(R * T / gamma)
    if M_2 >= 1:
        warnings.warn(
            f"Feed-line Mach number M_2 = {M_2:.4f} >= 1; the subsonic feed-line model is not valid.",
            RuntimeWarning,
            stacklevel=2,
        )

    A_factor = (2 / (gamma + 1) * (1 + (gamma - 1) / 2 * M_2**2)) ** (0.5 * (gamma + 1) / (gamma - 1))
    A_t = A_2 / A_factor * M_2
    D_t = 2 * np.sqrt(A_t / np.pi)

    T_t = T * (1 + (gamma - 1) / 2 * M_2**2)
    P_t_u = (
        m_dot / A_t * np.sqrt(T_t * R / gamma) * ((gamma + 1) / 2) ** ((gamma + 1) / (2 * (gamma - 1)))
    )

    return FeedLineNozzleResult(
        m_dot=float(m_dot),
        P_m=float(P_m),
        R=float(R),
        gamma=float(gamma),
        T=float(T),
        mass_flow_margin=float(mass_flow_margin),
        crit_ratio=float(crit_ratio),
        P_t_u_req=float(P_t_u_req),
        wall_thickness=float(wall_thickness_in),
        feed_diameter=float(feed_diameter),
        A_2=float(A_2),
        M_2=float(M_2),
        A_factor=float(A_factor),
        A_t=float(A_t),
        D_t=float(D_t),
        T_t=float(T_t),
        P_t_u=float(P_t_u),
        legacy_pressure_check_passed=bool(P_t_u > P_t_u_req),
    )


# ===========================================================================
# V1 - single-point sonic nozzle sizing
# ===========================================================================


@dataclass
class SingleNozzleResult:
    """Outputs of :func:`size_single_sonic_nozzle` (SonicNozzleSizingV1.m variables, SI)."""

    m_dot_req: float  # required mass flow [kg/s]
    flow_margin: float  # design mass-flow margin [-]
    m_dot: float  # design mass flow [kg/s]
    Cd: float  # discharge coefficient [-]
    P_down: float  # downstream/back pressure [Pa]
    R: float  # [J/(kg*K)]
    gamma: float  # [-]
    T_tu: float  # upstream total temperature [K]
    crit_ratio: float  # critical pressure ratio [-]
    P_tu_min: float  # minimum upstream total pressure for choking [Pa]
    choked_param: float  # choked mass-flow parameter
    A_t: float  # throat area [m^2]
    D_t_m: float  # throat diameter [m]
    T_throat: float  # throat static temperature at M=1 [K]
    P_throat: float  # throat static pressure at M=1 [Pa]
    wall_thickness: float  # feed-line wall thickness [m]
    feed_OD: float  # feed-line outer diameter [m]
    feed_ID: float  # feed-line inner diameter [m]

    def summary(self) -> str:
        """Command-window printout of SonicNozzleSizingV1.m."""
        lines = [
            "\n",
            "=============================================\n",
            "       SONIC NOZZLE SIZING RESULTS\n",
            "=============================================\n",
            "\n--- DESIGN REQUIREMENTS ---\n",
            _fprintf("Required Mass Flow          : %.3f lb/s\n", self.m_dot_req / LB_TO_KG),
            _fprintf("Design Mass Flow            : %.3f lb/s\n", self.m_dot / LB_TO_KG),
            _fprintf("Design Mass Flow            : %.3f kg/s\n", self.m_dot),
            _fprintf("Downstream Pressure         : %.2f psi\n", self.P_down / PSI_TO_PA),
            _fprintf("Downstream Pressure         : %.2f bar\n", self.P_down / 1e5),
            _fprintf("Upstream Total Temperature  : %.2f K\n", self.T_tu),
            "\n--- CHOKING CONDITION ---\n",
            _fprintf("Critical Pressure Ratio     : %.4f\n", self.crit_ratio),
            _fprintf("Minimum Upstream Total P    : %.2f psi\n", self.P_tu_min / PSI_TO_PA),
            _fprintf("Minimum Upstream Total P    : %.2f bar\n", self.P_tu_min / 1e5),
            "\n--- THROAT ---\n",
            _fprintf("Throat Diameter             : %.2f thou\n", (self.D_t_m / IN_TO_M) * 1000),
            _fprintf("Throat Diameter             : %.4f in\n", self.D_t_m / IN_TO_M),
            _fprintf("Throat Diameter             : %.3f mm\n", self.D_t_m * 1000),
            "\n--- THROAT STATIC CONDITIONS ---\n",
            _fprintf("Throat Pressure (M=1)       : %.2f psi\n", self.P_throat / PSI_TO_PA),
            _fprintf("Throat Pressure (M=1)       : %.2f bar\n", self.P_throat / 1e5),
            _fprintf("Throat Temperature (M=1)    : %.2f K\n", self.T_throat),
            "\n--- FEED LINE ---\n",
            _fprintf("Feed OD                     : %.3f in\n", self.feed_OD / IN_TO_M),
            _fprintf("Feed ID                     : %.3f in\n", self.feed_ID / IN_TO_M),
            "=============================================\n",
        ]
        return "".join(lines)


def size_single_sonic_nozzle(
    m_dot_req: float = 0.9072,
    flow_margin: float = 1.15,
    Cd: float = 1.0,
    P_down: float = 20e5 + 20 * PSI_TO_PA,
    R: float = 287.05,
    gamma: float = 1.4,
    T_tu: float = 283.0,
    wall_thickness_in: float = 0.049,
    feed_OD_in: float = 16 / 16,
) -> SingleNozzleResult:
    """Size a sonic nozzle from a downstream pressure and mass-flow requirement (V1).

    The nozzle is sized at the limiting choking condition: the minimum upstream
    total pressure is ``P_down / crit_ratio`` and the throat passes
    ``m_dot_req * flow_margin`` there.

    Parameters
    ----------
    m_dot_req : required mass flow [kg/s] (default 0.9072 kg/s = 2 lb/s).
    flow_margin : positive design margin on the mass flow (1.15 = +15 %).
    Cd : discharge coefficient.
    P_down : downstream/back pressure [Pa].  The MATLAB default is
        ``20e5 + 20*psi_to_Pa`` (20 bar plus a 20 psi margin = 21.38 bar,
        although the MATLAB comment says "= 20 bar").
    R, gamma, T_tu : gas constant [J/(kg*K)], specific heat ratio, upstream
        total temperature [K].
    wall_thickness_in, feed_OD_in : feed-line tube wall and outer diameter [in].
    """
    m_dot = m_dot_req * flow_margin

    crit_ratio = critical_pressure_ratio(gamma)
    P_tu_min = P_down / crit_ratio
    choked_param = choked_flow_parameter(gamma, R)

    A_t = (m_dot * np.sqrt(T_tu)) / (Cd * P_tu_min * choked_param)
    D_t_m = 2 * np.sqrt(A_t / np.pi)

    T_throat = T_tu * (2 / (gamma + 1))
    P_throat = P_tu_min * crit_ratio

    wall_thickness = wall_thickness_in * IN_TO_M
    feed_OD = feed_OD_in * IN_TO_M
    feed_ID = feed_OD - 2 * wall_thickness

    return SingleNozzleResult(
        m_dot_req=m_dot_req,
        flow_margin=flow_margin,
        m_dot=m_dot,
        Cd=Cd,
        P_down=P_down,
        R=R,
        gamma=gamma,
        T_tu=T_tu,
        crit_ratio=float(crit_ratio),
        P_tu_min=float(P_tu_min),
        choked_param=float(choked_param),
        A_t=float(A_t),
        D_t_m=float(D_t_m),
        T_throat=float(T_throat),
        P_throat=float(P_throat),
        wall_thickness=wall_thickness,
        feed_OD=feed_OD,
        feed_ID=feed_ID,
    )


# ===========================================================================
# V2 - fixed sonic orifice per gas against a fixed downstream pressure
# ===========================================================================


@dataclass
class FixedSonicOrificeResult:
    """Outputs of :func:`size_fixed_sonic_orifice` (SonicNozzleSizingV2.m).

    Per-gas arrays are ordered like :data:`GAS_NAMES`; per-point arrays are
    ``(num_rows, 4)``.  Pressures in Pa unless the name says otherwise.
    """

    table: pd.DataFrame
    gas_names: tuple[str, ...]
    pressure_margin: float
    Cd_all: np.ndarray
    P_down_all: np.ndarray
    crit_ratio_all: np.ndarray
    P0_critical_all: np.ndarray
    P0_design_min_all: np.ndarray
    orifice_area: np.ndarray  # [m^2]
    orifice_diameter_m: np.ndarray  # [m]
    P0_bar_all: np.ndarray
    P0_psi_all: np.ndarray
    choked_all: np.ndarray

    def summary(self, output_file: str | Path | None = None) -> str:
        """Command-window printout of SonicNozzleSizingV2.m."""
        rule = "===============================================================================\n"
        lines = [
            "\n",
            rule,
            "                   FIXED SONIC ORIFICE SIZING SUMMARY\n",
            rule,
            "\n",
            _fprintf("Pressure Margin Factor : %.3f\n", self.pressure_margin),
            "\n",
            "Gas       Cd      Pdown [bar]     D [mm]       D [in]      "
            "P0crit [bar]   P0,min [bar]\n",
            "-------------------------------------------------------------------------------\n",
        ]
        for g, gas in enumerate(self.gas_names):
            lines.append(
                _fprintf(
                    "%-5s   %6.3f     %9.2f     %9.3f     %9.4f     %11.2f     %11.2f\n",
                    gas,
                    self.Cd_all[g],
                    self.P_down_all[g] / 1e5,
                    self.orifice_diameter_m[g] * 1000,
                    self.orifice_diameter_m[g] / IN_TO_M,
                    self.P0_critical_all[g] / 1e5,
                    self.P0_design_min_all[g] / 1e5,
                )
            )
        lines += [rule, "\n"]
        if output_file is not None:
            lines.append(f"Results written to:\n{output_file}\n\n")
        return "".join(lines)


def size_fixed_sonic_orifice(
    data: pd.DataFrame | str | Path = FLOW_REQUIREMENTS_CSV,
    Cd: Mapping[str, float] | None = None,
    pressure_margin: float = 1.15,
    P_down_bar: Mapping[str, float] | None = None,
    P_margin_psi: Mapping[str, float] | None = None,
    gas_properties: Mapping[str, GasProperties] | None = None,
) -> FixedSonicOrificeResult:
    """Size ONE fixed sonic orifice per gas from a flow-requirements CSV (V2).

    For each gas the orifice is sized at the minimum positive mass flow with an
    upstream stagnation pressure ``pressure_margin`` times the choking limit
    ``P_down / crit_ratio``; the upstream pressure needed for every other flow
    point through that fixed orifice is then computed and checked for choking.

    Parameters
    ----------
    data : flow-requirements table or CSV path (columns H2, Air, CH4, GOx [kg/s]).
    Cd : sonic-orifice discharge coefficient per gas (default 0.85 for all).
    pressure_margin : factor above the theoretical choking pressure at minimum flow.
    P_down_bar : target downstream pressure per gas [bar abs]
        (default H2 20, Air 20, CH4 24, GOx 36).
    P_margin_psi : additional downstream pressure margin per gas [psi] (default 20).
    gas_properties : per-gas :class:`GasProperties` overrides.
    """
    table = _as_table(data)
    m_dot_all = _mass_flows(table)
    Cd = _merged(V2_CD, Cd)
    P_down_bar = _merged(V2_P_DOWN_BAR, P_down_bar)
    P_margin_psi = _merged(V2_P_MARGIN_PSI, P_margin_psi)
    props = _gas_props(gas_properties)

    num_rows = len(table)
    n = len(GAS_NAMES)
    P0_bar_all = np.full((num_rows, n), np.nan)
    P0_psi_all = np.full((num_rows, n), np.nan)
    choked_all = np.zeros((num_rows, n), dtype=bool)
    orifice_area = np.full(n, np.nan)
    orifice_diameter_m = np.full(n, np.nan)
    P0_critical_all = np.full(n, np.nan)
    P0_design_min_all = np.full(n, np.nan)
    crit_ratio_all = np.full(n, np.nan)

    Cd_all = np.array([Cd[gas] for gas in GAS_NAMES], dtype=float)
    P_down_all = np.array(
        [P_down_bar[gas] * 1e5 + P_margin_psi[gas] * PSI_TO_PA for gas in GAS_NAMES], dtype=float
    )

    for g, gas in enumerate(GAS_NAMES):
        m_dot = m_dot_all[gas]
        R, gamma, T0 = props[gas].R, props[gas].gamma, props[gas].T0
        Cd_g = Cd_all[g]
        P_down = P_down_all[g]

        valid = _valid(m_dot)
        if not valid.any():
            warnings.warn(f"No positive flow rates found for {gas}.", RuntimeWarning, stacklevel=2)
            continue

        m_dot_min = np.min(m_dot[valid])

        crit_ratio = critical_pressure_ratio(gamma)
        crit_ratio_all[g] = crit_ratio

        P0_critical = P_down / crit_ratio
        P0_critical_all[g] = P0_critical

        P0_design_min = pressure_margin * P0_critical
        P0_design_min_all[g] = P0_design_min

        choked_param = choked_flow_parameter(gamma, R)

        A_t = choked_area(m_dot_min, Cd_g, P0_design_min, T0, choked_param)
        orifice_area[g] = A_t
        orifice_diameter_m[g] = area_to_diameter(A_t)

        P0_required = np.full(m_dot.shape, np.nan)
        P0_required[valid] = choked_upstream_pressure(m_dot[valid], Cd_g, A_t, T0, choked_param)

        is_choked = np.zeros(m_dot.shape, dtype=bool)
        is_choked[valid] = (P_down / P0_required[valid]) <= crit_ratio

        P0_bar_all[:, g] = P0_required / 1e5
        P0_psi_all[:, g] = P0_required / PSI_TO_PA
        choked_all[:, g] = is_choked

    for g, gas in enumerate(GAS_NAMES):
        table[f"{gas}_P0_bar"] = P0_bar_all[:, g]
        table[f"{gas}_P0_psi"] = P0_psi_all[:, g]
        table[f"{gas}_Choked"] = choked_all[:, g]
    for g, gas in enumerate(GAS_NAMES):
        table[f"{gas}_Pdown_bar"] = np.full(num_rows, P_down_all[g] / 1e5)
    for g, gas in enumerate(GAS_NAMES):
        table[f"{gas}_Orifice_mm"] = np.full(num_rows, orifice_diameter_m[g] * 1000)
        table[f"{gas}_Orifice_in"] = np.full(num_rows, orifice_diameter_m[g] / IN_TO_M)

    return FixedSonicOrificeResult(
        table=table,
        gas_names=GAS_NAMES,
        pressure_margin=pressure_margin,
        Cd_all=Cd_all,
        P_down_all=P_down_all,
        crit_ratio_all=crit_ratio_all,
        P0_critical_all=P0_critical_all,
        P0_design_min_all=P0_design_min_all,
        orifice_area=orifice_area,
        orifice_diameter_m=orifice_diameter_m,
        P0_bar_all=P0_bar_all,
        P0_psi_all=P0_psi_all,
        choked_all=choked_all,
    )


# ===========================================================================
# V3 - two-stage choked sizing (injector area, then sonic orifice)
# ===========================================================================


@dataclass
class TwoStageChokedResult:
    """Outputs of :func:`size_two_stage_choked` (SonicNozzleSizingV3.m).

    Per-gas arrays are ordered like :data:`GAS_NAMES`; per-point arrays are
    ``(num_rows, 4)``.  Pressures in Pa, areas in m^2, diameters in m unless
    the name says otherwise.
    """

    table: pd.DataFrame
    gas_names: tuple[str, ...]
    P_chamber: float
    injector_pressure_margin: float
    orifice_pressure_margin: float
    crit_ratio_all: np.ndarray
    Pman_critical_all: np.ndarray
    Pman_design_min_all: np.ndarray
    injector_area: np.ndarray
    injector_diameter_m: np.ndarray
    P0_orifice_critical_all: np.ndarray
    P0_orifice_design_min_all: np.ndarray
    orifice_area: np.ndarray
    orifice_diameter_m: np.ndarray
    P_manifold_bar_all: np.ndarray
    P_manifold_psi_all: np.ndarray
    P0_bar_all: np.ndarray
    P0_psi_all: np.ndarray
    injector_choked_all: np.ndarray
    orifice_choked_all: np.ndarray

    def summary(self, output_file: str | Path | None = None) -> str:
        """Command-window printout of SonicNozzleSizingV3.m."""
        rule = "=" * 94 + "\n"
        lines = [
            "\n",
            rule,
            "                        TWO-STAGE CHOKED FLOW SIZING SUMMARY\n",
            rule,
            "\n",
            _fprintf("Chamber Pressure              : %.2f bar absolute\n", self.P_chamber / 1e5),
            _fprintf("Injector Pressure Margin      : %.3f\n", self.injector_pressure_margin),
            _fprintf("Orifice Pressure Margin       : %.3f\n", self.orifice_pressure_margin),
            "\n",
            "Gas      Pman,min    Inj Area     Inj Eq D      Orifice D     P0,min\n",
            "         [bar]       [mm^2]       [mm]          [mm]          [bar]\n",
            "-" * 94 + "\n",
        ]
        for g, gas in enumerate(self.gas_names):
            lines.append(
                _fprintf(
                    "%-5s   %9.2f   %10.3f   %10.3f   %11.3f   %10.2f\n",
                    gas,
                    self.Pman_design_min_all[g] / 1e5,
                    self.injector_area[g] * 1e6,
                    self.injector_diameter_m[g] * 1000,
                    self.orifice_diameter_m[g] * 1000,
                    self.P0_orifice_design_min_all[g] / 1e5,
                )
            )
        lines += [rule, "\n"]
        if output_file is not None:
            lines.append(f"Results written to:\n{output_file}\n\n")
        return "".join(lines)


def size_two_stage_choked(
    data: pd.DataFrame | str | Path = FLOW_REQUIREMENTS_CSV,
    P_chamber_bar: float = 10.0,
    Cd_inj: Mapping[str, float] | None = None,
    Cd_orifice: Mapping[str, float] | None = None,
    injector_pressure_margin: float = 1.15,
    orifice_pressure_margin: float = 1.15,
    gas_properties: Mapping[str, GasProperties] | None = None,
) -> TwoStageChokedResult:
    """Two-stage choked-flow sizing from a flow-requirements CSV (V3).

    Stage 1 sizes one fixed injector effective area per gas so that, at the
    minimum positive flow, the manifold sits ``injector_pressure_margin`` above
    the pressure needed to choke into the chamber.  Stage 2 sizes one fixed
    upstream sonic orifice that chokes into that manifold pressure with
    ``orifice_pressure_margin``.  Manifold and supply pressures and both choking
    checks are then evaluated at every flow point.

    Parameters
    ----------
    data : flow-requirements table or CSV path (columns H2, Air, CH4, GOx [kg/s]).
    P_chamber_bar : chamber pressure [bar abs], same for every gas.
    Cd_inj, Cd_orifice : injector / sonic-orifice discharge coefficients per gas
        (default 0.85 for all).
    injector_pressure_margin, orifice_pressure_margin : factors above the
        theoretical choking pressures at minimum flow.
    gas_properties : per-gas :class:`GasProperties` overrides.
    """
    table = _as_table(data)
    m_dot_all = _mass_flows(table)
    Cd_inj = _merged(V3_CD_INJ, Cd_inj)
    Cd_orifice = _merged(V3_CD_ORIFICE, Cd_orifice)
    props = _gas_props(gas_properties)

    P_chamber = P_chamber_bar * 1e5

    num_rows = len(table)
    n = len(GAS_NAMES)
    P_manifold_bar_all = np.full((num_rows, n), np.nan)
    P_manifold_psi_all = np.full((num_rows, n), np.nan)
    P0_bar_all = np.full((num_rows, n), np.nan)
    P0_psi_all = np.full((num_rows, n), np.nan)
    injector_choked_all = np.zeros((num_rows, n), dtype=bool)
    orifice_choked_all = np.zeros((num_rows, n), dtype=bool)
    injector_area = np.full(n, np.nan)
    injector_diameter_m = np.full(n, np.nan)
    orifice_area = np.full(n, np.nan)
    orifice_diameter_m = np.full(n, np.nan)
    Pman_critical_all = np.full(n, np.nan)
    Pman_design_min_all = np.full(n, np.nan)
    P0_orifice_critical_all = np.full(n, np.nan)
    P0_orifice_design_min_all = np.full(n, np.nan)
    crit_ratio_all = np.full(n, np.nan)

    for g, gas in enumerate(GAS_NAMES):
        m_dot = m_dot_all[gas]
        R, gamma, T0 = props[gas].R, props[gas].gamma, props[gas].T0
        Cd_inj_g = Cd_inj[gas]
        Cd_orifice_g = Cd_orifice[gas]

        valid = _valid(m_dot)
        if not valid.any():
            warnings.warn(f"No positive flow rates found for {gas}.", RuntimeWarning, stacklevel=2)
            continue

        m_dot_min = np.min(m_dot[valid])

        crit_ratio = critical_pressure_ratio(gamma)
        crit_ratio_all[g] = crit_ratio
        choked_param = choked_flow_parameter(gamma, R)

        # Stage 1: chamber -> injector -> manifold
        Pman_critical = P_chamber / crit_ratio
        Pman_critical_all[g] = Pman_critical
        Pman_design_min = injector_pressure_margin * Pman_critical
        Pman_design_min_all[g] = Pman_design_min

        A_inj = choked_area(m_dot_min, Cd_inj_g, Pman_design_min, T0, choked_param)
        injector_area[g] = A_inj
        injector_diameter_m[g] = area_to_diameter(A_inj)

        P_manifold_required = np.full(m_dot.shape, np.nan)
        P_manifold_required[valid] = choked_upstream_pressure(
            m_dot[valid], Cd_inj_g, A_inj, T0, choked_param
        )

        injector_choked = np.zeros(m_dot.shape, dtype=bool)
        injector_choked[valid] = (P_chamber / P_manifold_required[valid]) <= crit_ratio

        # Stage 2: manifold -> sonic orifice -> supply
        P_down_orifice_min = Pman_design_min
        P0_orifice_critical = P_down_orifice_min / crit_ratio
        P0_orifice_critical_all[g] = P0_orifice_critical
        P0_orifice_design_min = orifice_pressure_margin * P0_orifice_critical
        P0_orifice_design_min_all[g] = P0_orifice_design_min

        A_orifice = choked_area(m_dot_min, Cd_orifice_g, P0_orifice_design_min, T0, choked_param)
        orifice_area[g] = A_orifice
        orifice_diameter_m[g] = area_to_diameter(A_orifice)

        P0_required = np.full(m_dot.shape, np.nan)
        P0_required[valid] = choked_upstream_pressure(
            m_dot[valid], Cd_orifice_g, A_orifice, T0, choked_param
        )

        # Orifice back pressure is the (flow-dependent) manifold pressure.
        orifice_choked = np.zeros(m_dot.shape, dtype=bool)
        orifice_choked[valid] = (P_manifold_required[valid] / P0_required[valid]) <= crit_ratio

        P_manifold_bar_all[:, g] = P_manifold_required / 1e5
        P_manifold_psi_all[:, g] = P_manifold_required / PSI_TO_PA
        P0_bar_all[:, g] = P0_required / 1e5
        P0_psi_all[:, g] = P0_required / PSI_TO_PA
        injector_choked_all[:, g] = injector_choked
        orifice_choked_all[:, g] = orifice_choked

    for g, gas in enumerate(GAS_NAMES):
        table[f"{gas}_Manifold_bar"] = P_manifold_bar_all[:, g]
        table[f"{gas}_Manifold_psi"] = P_manifold_psi_all[:, g]
        table[f"{gas}_Injector_Choked"] = injector_choked_all[:, g]
    for g, gas in enumerate(GAS_NAMES):
        table[f"{gas}_P0_bar"] = P0_bar_all[:, g]
        table[f"{gas}_P0_psi"] = P0_psi_all[:, g]
        table[f"{gas}_Orifice_Choked"] = orifice_choked_all[:, g]
    for g, gas in enumerate(GAS_NAMES):
        table[f"{gas}_InjectorArea_mm2"] = np.full(num_rows, injector_area[g] * 1e6)
    for g, gas in enumerate(GAS_NAMES):
        table[f"{gas}_InjectorEqDiameter_mm"] = np.full(num_rows, injector_diameter_m[g] * 1000)
    for g, gas in enumerate(GAS_NAMES):
        table[f"{gas}_Orifice_mm"] = np.full(num_rows, orifice_diameter_m[g] * 1000)
    for g, gas in enumerate(GAS_NAMES):
        table[f"{gas}_Orifice_in"] = np.full(num_rows, orifice_diameter_m[g] / IN_TO_M)

    return TwoStageChokedResult(
        table=table,
        gas_names=GAS_NAMES,
        P_chamber=P_chamber,
        injector_pressure_margin=injector_pressure_margin,
        orifice_pressure_margin=orifice_pressure_margin,
        crit_ratio_all=crit_ratio_all,
        Pman_critical_all=Pman_critical_all,
        Pman_design_min_all=Pman_design_min_all,
        injector_area=injector_area,
        injector_diameter_m=injector_diameter_m,
        P0_orifice_critical_all=P0_orifice_critical_all,
        P0_orifice_design_min_all=P0_orifice_design_min_all,
        orifice_area=orifice_area,
        orifice_diameter_m=orifice_diameter_m,
        P_manifold_bar_all=P_manifold_bar_all,
        P_manifold_psi_all=P_manifold_psi_all,
        P0_bar_all=P0_bar_all,
        P0_psi_all=P0_psi_all,
        injector_choked_all=injector_choked_all,
        orifice_choked_all=orifice_choked_all,
    )


# ===========================================================================
# V4 - manifold pressure from known injector area
# ===========================================================================


@dataclass
class ManifoldPressureResult:
    """Outputs of :func:`manifold_pressure_from_injector` (SonicNozzleSizingV4.m)."""

    table: pd.DataFrame
    gas_names: tuple[str, ...]
    P_chamber_bar: float
    A_mm2_all: np.ndarray  # injector area [mm^2]
    Cd_all: np.ndarray
    K_all: np.ndarray
    crit_all: np.ndarray
    Pman_min_choke_bar_all: np.ndarray  # minimum manifold pressure for choking [bar]
    Pman_bar_all: np.ndarray  # (num_rows, 4) [bar abs]
    choked_all: np.ndarray  # (num_rows, 4)

    def summary(self, output_file: str | Path | None = None) -> str:
        """Command-window printout of SonicNozzleSizingV4.m."""
        rule = "============================================================\n"
        lines = [
            "\n",
            rule,
            "              RDE MANIFOLD PRESSURE MODEL\n",
            rule,
            "\n",
            _fprintf("Chamber Pressure: %.2f bar absolute\n", self.P_chamber_bar),
            "\n",
            "Minimum manifold pressures for ideal choking:\n",
        ]
        labels = {"H2": "H2  ", "Air": "Air ", "CH4": "CH4 ", "GOx": "GOx "}
        for g, gas in enumerate(self.gas_names):
            lines.append(_fprintf(f"{labels[gas]}: %.3f bar\n", self.Pman_min_choke_bar_all[g]))
        lines += ["\n", "Injector areas:\n"]
        for g, gas in enumerate(self.gas_names):
            lines.append(_fprintf(f"{labels[gas]}: %.2f mm^2\n", self.A_mm2_all[g]))
        lines += ["\n", "Effective injector Cd values:\n"]
        for g, gas in enumerate(self.gas_names):
            lines.append(_fprintf(f"{labels[gas]}: %.3f\n", self.Cd_all[g]))
        lines.append("\n")
        if output_file is not None:
            lines.append(f"Results written to:\n{output_file}\n")
        lines.append(rule)
        return "".join(lines)


def manifold_pressure_from_injector(
    data: pd.DataFrame | str | Path = MANIFOLD_INPUT_CSV,
    P_chamber_bar: float = 10.0,
    A_mm2: Mapping[str, float] | None = None,
    Cd: Mapping[str, float] | None = None,
    gas_properties: Mapping[str, GasProperties] | None = None,
) -> ManifoldPressureResult:
    """Required RDE manifold pressure from injector area and mass flow (V4).

    Assumes each injector array is choked:
    ``Pman = m_dot * sqrt(T0) / (Cd * A * K)``.  ``P_chamber_bar`` is used only
    to check whether the injector actually chokes at that manifold pressure.
    Unlike V2/V3/V5 there is no positive-flow filter: every row is evaluated and
    only NaN rows are forced to "not choked" (as in the MATLAB).

    Parameters
    ----------
    data : flow-requirements table or CSV path (columns H2, Air, CH4, GOx [kg/s]).
    P_chamber_bar : mean chamber pressure [bar abs].
    A_mm2 : total physical injector flow area per gas [mm^2].
    Cd : effective injector discharge coefficient per gas.
    gas_properties : per-gas :class:`GasProperties` overrides.
    """
    table = _as_table(data)
    m_dot_all = _mass_flows(table)
    A_mm2 = _merged(INJECTOR_AREA_MM2, A_mm2)
    Cd = _merged(INJECTOR_CD, Cd)
    props = _gas_props(gas_properties)

    num_rows = len(table)
    n = len(GAS_NAMES)
    A_mm2_all = np.array([A_mm2[gas] for gas in GAS_NAMES], dtype=float)
    Cd_all = np.array([Cd[gas] for gas in GAS_NAMES], dtype=float)
    K_all = np.full(n, np.nan)
    crit_all = np.full(n, np.nan)
    Pman_min_choke_bar_all = np.full(n, np.nan)
    Pman_bar_all = np.full((num_rows, n), np.nan)
    choked_all = np.zeros((num_rows, n), dtype=bool)

    for g, gas in enumerate(GAS_NAMES):
        m_dot = m_dot_all[gas]
        R, gamma, T0 = props[gas].R, props[gas].gamma, props[gas].T0
        A = A_mm2_all[g] * 1e-6

        K = choked_flow_parameter(gamma, R)
        K_all[g] = K
        Pman = choked_upstream_pressure(m_dot, Cd_all[g], A, T0, K)
        Pman_bar = Pman / 1e5

        crit = critical_pressure_ratio(gamma)
        crit_all[g] = crit
        Pman_min_choke_bar_all[g] = P_chamber_bar / crit

        with np.errstate(divide="ignore", invalid="ignore"):
            choked = (P_chamber_bar / Pman_bar) <= crit
        choked[np.isnan(m_dot)] = False

        Pman_bar_all[:, g] = Pman_bar
        choked_all[:, g] = choked

    for g, gas in enumerate(GAS_NAMES):
        table[f"{gas}_ManifoldPressure_bar"] = Pman_bar_all[:, g]
        table[f"{gas}_Choked"] = choked_all[:, g]

    return ManifoldPressureResult(
        table=table,
        gas_names=GAS_NAMES,
        P_chamber_bar=P_chamber_bar,
        A_mm2_all=A_mm2_all,
        Cd_all=Cd_all,
        K_all=K_all,
        crit_all=crit_all,
        Pman_min_choke_bar_all=Pman_min_choke_bar_all,
        Pman_bar_all=Pman_bar_all,
        choked_all=choked_all,
    )


# ===========================================================================
# V5 - injector + fixed sonic orifice (canonical)
# ===========================================================================


@dataclass
class InjectorSonicOrificeResult:
    """Outputs of :func:`size_injector_sonic_orifice` (SonicNozzleSizingV5.m).

    Per-gas arrays are ordered like :data:`GAS_NAMES`; per-point arrays are
    ``(num_rows, 4)``.  Pressures in Pa unless the name says otherwise.
    """

    table: pd.DataFrame
    gas_names: tuple[str, ...]
    P_chamber_bar: float
    pressure_margin: float
    Cd_inj_all: np.ndarray
    Cd_so_all: np.ndarray
    valid_all: np.ndarray  # (num_rows, 4) positive, non-NaN flow points
    critical_ratio_all: np.ndarray
    min_flow_all: np.ndarray  # [kg/s]
    Pman_minflow_all: np.ndarray
    P0_critical_minflow_all: np.ndarray
    P0_design_minflow_all: np.ndarray
    sonic_area_all: np.ndarray  # [m^2]
    sonic_diameter_all: np.ndarray  # [m]
    Pman_bar_all: np.ndarray
    injector_choked_all: np.ndarray
    P0_supply_bar_all: np.ndarray
    P0_supply_psia_all: np.ndarray
    P0_supply_psig_all: np.ndarray
    sonic_choked_all: np.ndarray
    model_valid_all: np.ndarray

    def summary(self, output_file: str | Path | None = None) -> str:
        """Command-window printout of SonicNozzleSizingV5.m."""
        rule = "=" * 88 + "\n"
        lines = [
            "\n",
            rule,
            "             RDE INJECTOR + FIXED SONIC ORIFICE SIZING SUMMARY\n",
            rule,
            "\n",
            _fprintf("Chamber Pressure             : %.3f bar absolute\n", self.P_chamber_bar),
            _fprintf("Sonic Orifice Pressure Margin: %.3f\n", self.pressure_margin),
            "\n",
            "Gas      Cd_inj   Cd_SO    Min Flow     Pman@Min    "
            "SO D [mm]   SO D [in]   P0@Min [bar]\n",
            "                           [kg/s]        [bar abs]    "
            "                        [bar abs]\n",
            "-" * 88 + "\n",
        ]
        for g, gas in enumerate(self.gas_names):
            lines.append(
                _fprintf(
                    "%-5s    %6.3f   %6.3f    %8.4f     %8.3f     %8.3f    %8.4f     %8.3f\n",
                    gas,
                    self.Cd_inj_all[g],
                    self.Cd_so_all[g],
                    self.min_flow_all[g],
                    self.Pman_minflow_all[g] / 1e5,
                    self.sonic_diameter_all[g] * 1000,
                    self.sonic_diameter_all[g] / IN_TO_M,
                    self.P0_design_minflow_all[g] / 1e5,
                )
            )
        lines += [rule, "\n"]
        for g, gas in enumerate(self.gas_names):
            valid = self.valid_all[:, g]
            n_valid = int(valid.sum())
            lines += [
                f"{gas}:\n",
                f"  Injector choked at {int(self.injector_choked_all[valid, g].sum())} "
                f"of {n_valid} valid points\n",
                f"  Sonic orifice choked at {int(self.sonic_choked_all[valid, g].sum())} "
                f"of {n_valid} valid points\n",
                f"  Full simplified model valid at {int(self.model_valid_all[valid, g].sum())} "
                f"of {n_valid} valid points\n",
                "\n",
            ]
        if output_file is not None:
            lines.append(f"Results written to:\n{output_file}\n\n")
        return "".join(lines)


def size_injector_sonic_orifice(
    data: pd.DataFrame | str | Path = MANIFOLD_INPUT_CSV,
    P_chamber_bar: float = 10.0,
    pressure_margin: float = 1.15,
    A_inj_mm2: Mapping[str, float] | None = None,
    Cd_inj: Mapping[str, float] | None = None,
    Cd_so: Mapping[str, float] | None = None,
    gas_properties: Mapping[str, GasProperties] | None = None,
    P_atm: float = 101325.0,
) -> InjectorSonicOrificeResult:
    """Combined RDE injector + fixed sonic orifice model (V5, canonical).

    1. Required manifold pressure at every flow point from the (known) injector
       area and Cd, assuming a choked injector; warn where the injector is not
       actually choked against ``P_chamber_bar``.
    2. At the minimum positive flow, size ONE fixed sonic orifice per gas whose
       upstream pressure is ``pressure_margin`` above the pressure needed to
       choke into that manifold pressure.
    3. Required supply stagnation pressure through that orifice at every flow
       point, its choking check, and overall model validity (both choked).

    Parameters
    ----------
    data : flow-requirements table or CSV path (columns H2, Air, CH4, GOx [kg/s]).
    P_chamber_bar : mean chamber pressure [bar abs].
    pressure_margin : sonic-orifice pressure margin at minimum flow.
    A_inj_mm2 : total physical injector flow area per gas [mm^2].
    Cd_inj : effective injector discharge coefficient per gas.
    Cd_so : sonic-orifice discharge coefficient per gas (default 0.85).
    gas_properties : per-gas :class:`GasProperties` overrides.
    P_atm : atmospheric pressure [Pa], only for absolute -> gauge conversion.
    """
    table = _as_table(data)
    m_dot_all = _mass_flows(table)
    A_inj_mm2 = _merged(INJECTOR_AREA_MM2, A_inj_mm2)
    Cd_inj = _merged(INJECTOR_CD, Cd_inj)
    Cd_so = _merged(SONIC_ORIFICE_CD, Cd_so)
    props = _gas_props(gas_properties)

    num_rows = len(table)
    n = len(GAS_NAMES)
    A_inj_all = np.array([A_inj_mm2[gas] * 1e-6 for gas in GAS_NAMES], dtype=float)
    Cd_inj_all = np.array([Cd_inj[gas] for gas in GAS_NAMES], dtype=float)
    Cd_so_all = np.array([Cd_so[gas] for gas in GAS_NAMES], dtype=float)

    valid_all = np.zeros((num_rows, n), dtype=bool)
    Pman_bar_all = np.full((num_rows, n), np.nan)
    injector_choked_all = np.zeros((num_rows, n), dtype=bool)
    P0_supply_bar_all = np.full((num_rows, n), np.nan)
    P0_supply_psia_all = np.full((num_rows, n), np.nan)
    P0_supply_psig_all = np.full((num_rows, n), np.nan)
    sonic_choked_all = np.zeros((num_rows, n), dtype=bool)
    model_valid_all = np.zeros((num_rows, n), dtype=bool)
    sonic_area_all = np.full(n, np.nan)
    sonic_diameter_all = np.full(n, np.nan)
    critical_ratio_all = np.full(n, np.nan)
    min_flow_all = np.full(n, np.nan)
    Pman_minflow_all = np.full(n, np.nan)
    P0_critical_minflow_all = np.full(n, np.nan)
    P0_design_minflow_all = np.full(n, np.nan)

    for g, gas in enumerate(GAS_NAMES):
        m_dot = m_dot_all[gas]
        A_inj = A_inj_all[g]
        Cd_inj_g = Cd_inj_all[g]
        Cd_so_g = Cd_so_all[g]
        R, gamma, T0 = props[gas].R, props[gas].gamma, props[gas].T0

        valid = _valid(m_dot)
        valid_all[:, g] = valid
        if not valid.any():
            warnings.warn(
                f"No valid positive mass flow rates for {gas}.", RuntimeWarning, stacklevel=2
            )
            continue

        K = choked_flow_parameter(gamma, R)
        critical_ratio = critical_pressure_ratio(gamma)
        critical_ratio_all[g] = critical_ratio

        # Step 1: required manifold pressure (choked injector)
        Pman = np.full(m_dot.shape, np.nan)
        Pman[valid] = choked_upstream_pressure(m_dot[valid], Cd_inj_g, A_inj, T0, K)
        Pman_bar = Pman / 1e5

        injector_choked = np.zeros(m_dot.shape, dtype=bool)
        injector_choked[valid] = (P_chamber_bar / Pman_bar[valid]) <= critical_ratio

        bad_injector = valid & ~injector_choked
        if bad_injector.any():
            warnings.warn(
                f"{gas} injector is NOT choked at {int(bad_injector.sum())} test point(s). "
                "The choked-injector manifold pressure equation is not valid at those points.",
                RuntimeWarning,
                stacklevel=2,
            )

        # Step 2: minimum mass-flow point
        valid_indices = np.flatnonzero(valid)
        local_min_index = int(np.argmin(m_dot[valid]))
        m_dot_min = m_dot[valid][local_min_index]
        min_index = valid_indices[local_min_index]
        Pman_min = Pman[min_index]
        min_flow_all[g] = m_dot_min
        Pman_minflow_all[g] = Pman_min

        # Steps 3-4: sonic-orifice choking limit + margin
        P0_critical_min = Pman_min / critical_ratio
        P0_critical_minflow_all[g] = P0_critical_min
        P0_design_min = pressure_margin * P0_critical_min
        P0_design_minflow_all[g] = P0_design_min

        # Step 5: one fixed sonic orifice
        A_so = choked_area(m_dot_min, Cd_so_g, P0_design_min, T0, K)
        sonic_area_all[g] = A_so
        sonic_diameter_all[g] = area_to_diameter(A_so)

        # Step 6: supply pressure through the fixed orifice
        P0_supply = np.full(m_dot.shape, np.nan)
        P0_supply[valid] = choked_upstream_pressure(m_dot[valid], Cd_so_g, A_so, T0, K)

        # Step 7: sonic-orifice choking (back pressure = manifold pressure)
        sonic_choked = np.zeros(m_dot.shape, dtype=bool)
        sonic_choked[valid] = (Pman[valid] / P0_supply[valid]) <= critical_ratio

        # Step 8: overall model validity
        model_valid = injector_choked & sonic_choked

        Pman_bar_all[:, g] = Pman_bar
        injector_choked_all[:, g] = injector_choked
        P0_supply_bar_all[:, g] = P0_supply / 1e5
        P0_supply_psia_all[:, g] = P0_supply / PSI_TO_PA_V5
        P0_supply_psig_all[:, g] = (P0_supply - P_atm) / PSI_TO_PA_V5
        sonic_choked_all[:, g] = sonic_choked
        model_valid_all[:, g] = model_valid

    for g, gas in enumerate(GAS_NAMES):
        table[f"{gas}_ManifoldPressure_bar"] = Pman_bar_all[:, g]
        table[f"{gas}_InjectorChoked"] = injector_choked_all[:, g]
        table[f"{gas}_SupplyP0_bar_abs"] = P0_supply_bar_all[:, g]
        table[f"{gas}_SupplyP0_psia"] = P0_supply_psia_all[:, g]
        table[f"{gas}_SupplyP0_psig"] = P0_supply_psig_all[:, g]
        table[f"{gas}_SonicOrificeChoked"] = sonic_choked_all[:, g]
        table[f"{gas}_ModelValid"] = model_valid_all[:, g]
        table[f"{gas}_SonicOrificeArea_mm2"] = np.full(num_rows, sonic_area_all[g] * 1e6)
        table[f"{gas}_SonicOrificeDiameter_mm"] = np.full(num_rows, sonic_diameter_all[g] * 1000)
        table[f"{gas}_SonicOrificeDiameter_in"] = np.full(num_rows, sonic_diameter_all[g] / IN_TO_M)

    return InjectorSonicOrificeResult(
        table=table,
        gas_names=GAS_NAMES,
        P_chamber_bar=P_chamber_bar,
        pressure_margin=pressure_margin,
        Cd_inj_all=Cd_inj_all,
        Cd_so_all=Cd_so_all,
        valid_all=valid_all,
        critical_ratio_all=critical_ratio_all,
        min_flow_all=min_flow_all,
        Pman_minflow_all=Pman_minflow_all,
        P0_critical_minflow_all=P0_critical_minflow_all,
        P0_design_minflow_all=P0_design_minflow_all,
        sonic_area_all=sonic_area_all,
        sonic_diameter_all=sonic_diameter_all,
        Pman_bar_all=Pman_bar_all,
        injector_choked_all=injector_choked_all,
        P0_supply_bar_all=P0_supply_bar_all,
        P0_supply_psia_all=P0_supply_psia_all,
        P0_supply_psig_all=P0_supply_psig_all,
        sonic_choked_all=sonic_choked_all,
        model_valid_all=model_valid_all,
    )


# ===========================================================================
# Script entry point
# ===========================================================================

#: model key -> (function, default input CSV, output-file rule)
MODELS = {
    "feedline": (size_sonic_nozzle_with_feed_line, None, None),
    "v1": (size_single_sonic_nozzle, None, None),
    "v2": (size_fixed_sonic_orifice, FLOW_REQUIREMENTS_CSV, "_sonic_orifice_results"),
    "v3": (size_two_stage_choked, FLOW_REQUIREMENTS_CSV, "_two_stage_choked_results"),
    "v4": (manifold_pressure_from_injector, MANIFOLD_INPUT_CSV, "manifold_pressure_results.csv"),
    "v5": (size_injector_sonic_orifice, MANIFOLD_INPUT_CSV, "_combined_results"),
}


def default_output_path(model: str, input_file: str | Path) -> Path:
    """Default results CSV for ``model``: ``<input name><suffix>.csv`` (V4: the fixed
    name ``manifold_pressure_results.csv``) in the CURRENT WORKING DIRECTORY.

    The bundled inputs live next to this module, but generated results go to the
    CWD so a run never overwrites the bundled reference CSVs or needs a writable
    package folder.  (MATLAB, run on its default bare-filename input, also wrote
    into ``pwd``.)
    """
    rule = MODELS[model][2]
    if rule is None:
        raise ValueError(f"Model {model!r} does not write a CSV.")
    name = rule if rule.endswith(".csv") else Path(input_file).stem + rule + ".csv"
    return Path.cwd() / name


def run(
    model: str = "v5",
    input_file: str | Path | None = None,
    output_file: str | Path | None = None,
    write: bool = True,
    verbose: bool = True,
):
    """Run one sizing workflow like the original script: compute, write CSV, print summary.

    ``input_file`` defaults to the bundled CSV next to this module; ``output_file``
    defaults to :func:`default_output_path` (current working directory).
    ``feedline`` and ``v1`` take no CSV and write nothing.  Returns the result
    dataclass.
    """
    model = model.lower()
    if model not in MODELS:
        raise ValueError(f"Unknown model {model!r}; choose from {sorted(MODELS)}.")
    func, default_input, _ = MODELS[model]

    if default_input is None:
        result = func()
        if verbose:
            print(result.summary(), end="")
        return result

    input_file = Path(input_file) if input_file is not None else default_input
    result = func(input_file)
    written = None
    if write:
        written = Path(output_file) if output_file is not None else default_output_path(model, input_file)
        write_results_csv(result.table, written)
    if verbose:
        print(result.summary(written), end="")
    return result


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="HADES feed-system sonic nozzle / orifice sizing "
        "(SonicNozzleSizing.py prototype = feedline, SonicNozzleSizingV1-V5.m = v1-v5)."
    )
    parser.add_argument(
        "model",
        nargs="?",
        default="v5",
        choices=sorted(MODELS),
        help="sizing workflow; default v5",
    )
    parser.add_argument(
        "--input", type=Path, help="flow-requirements CSV (H2, Air, CH4, GOx [kg/s]); default: bundled CSV"
    )
    parser.add_argument(
        "--output", type=Path, help="results CSV (default: ./<input name>_<suffix>.csv in the current directory)"
    )
    parser.add_argument("--no-write", action="store_true", help="do not write the results CSV")
    args = parser.parse_args(argv)
    run(args.model, args.input, args.output, write=not args.no_write)


if __name__ == "__main__":
    main()
