"""
Blowdown Sonic-Orifice Cd Calibration
=====================================

Reads a CSV containing:
    time_ms
    supply_pressure_psi
    upstream_pressure_psi
    downstream_pressure_psi
    supply_temperature_K
    upstream_temperature_K
    downstream_temperature_K

Calculates:
    - Orifice flow area
    - Upstream-to-downstream pressure drop
    - Supply and upstream real-gas density with CoolProp
    - Gas mass remaining in a rigid supply tank
    - Blowdown mass-flow rate from -dm_tank/dt
    - Real-gas compressible discharge coefficient, Cd
    - Ideal-gas compressible Cd for comparison
    - Incompressible-equivalent Cd for comparison
    - Choked / unchoked status

IMPORTANT:
    The real-gas Cd calculation treats the upstream pressure/temperature
    immediately before the orifice as stagnation conditions and computes
    the ideal isentropic nozzle mass flux with CoolProp.

Dependencies:
    pip install numpy pandas scipy CoolProp
"""

from pathlib import Path
import math

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from CoolProp.CoolProp import PropsSI


# ============================================================
# USER INPUTS
# ============================================================

from pathlib import Path

# Folder containing this Python script
SCRIPT_DIR = Path(__file__).resolve().parent

# Input/output files in the same FEED SIZING folder
INPUT_CSV = SCRIPT_DIR / "orifice_blowdown_input.csv"
OUTPUT_CSV = SCRIPT_DIR / "orifice_blowdown_results.csv"

# CoolProp fluid name.
# Common examples:
#   "Air"
#   "Oxygen"
#   "Methane"
#   "Hydrogen"
COOLPROP_FLUID = "Air"

# Orifice diameter
ORIFICE_DIAMETER_IN = 0.240

# Specific gas constant [J/(kg*K)]
#
# Used for the ideal-gas comparison calculation.
# Examples near room temperature:
#   Air      ~ 287.05
#   Oxygen   ~ 259.84
#   Methane  ~ 518.3
#   Hydrogen ~ 4124
GAS_CONSTANT_J_KG_K = 287.05

# Rigid supply-tank internal volume.
#
# This is REQUIRED because:
#
#     tank mass = real-gas density * tank volume
#
# Example below is 10 L. Replace with your actual tank volume.
TANK_VOLUME_L = 50.0

# Set True if ALL pressure columns in the CSV are psig.
# Set False if they are already psia.
PRESSURES_ARE_GAUGE = True

# Atmospheric pressure used for psig -> psia conversion.
ATM_PRESSURE_PSI = 14.6959

# Optional centered smoothing of calculated tank mass before
# differentiating it.
#
# 1 = no smoothing
# 3, 5, 7, ... = centered rolling-average window in rows
#
# Differentiating noisy pressure measurements can create very noisy mdot.
SMOOTHING_WINDOW_ROWS = 1


# ============================================================
# CONSTANTS / CONVERSIONS
# ============================================================

PSI_TO_PA = 6894.757293168
IN_TO_M = 0.0254
L_TO_M3 = 1.0e-3

ORIFICE_DIAMETER_M = ORIFICE_DIAMETER_IN * IN_TO_M
ORIFICE_AREA_M2 = math.pi * ORIFICE_DIAMETER_M**2 / 4.0
TANK_VOLUME_M3 = TANK_VOLUME_L * L_TO_M3


# ============================================================
# REQUIRED CSV COLUMNS
# ============================================================

REQUIRED_COLUMNS = [
    "time_ms",
    "supply_pressure_psi",
    "upstream_pressure_psi",
    "downstream_pressure_psi",
    "supply_temperature_K",
    "upstream_temperature_K",
    "downstream_temperature_K",
]


# ============================================================
# HELPER FUNCTIONS
# ============================================================

def pressure_psi_to_absolute_pa(p_psi):
    """
    Convert a pressure from the CSV to absolute Pa.
    """
    p_psi = np.asarray(p_psi, dtype=float)

    if PRESSURES_ARE_GAUGE:
        p_abs_psi = p_psi + ATM_PRESSURE_PSI
    else:
        p_abs_psi = p_psi

    return p_abs_psi * PSI_TO_PA


def coolprop_property(output, pressure_pa, temperature_k):
    """
    Calculate a CoolProp property from P,T for every row.

    Returns a NumPy array.

    Raises an informative error if CoolProp cannot evaluate a row.
    """
    values = np.empty(len(pressure_pa), dtype=float)

    for i, (p, t) in enumerate(zip(pressure_pa, temperature_k)):
        try:
            values[i] = PropsSI(
                output,
                "P",
                float(p),
                "T",
                float(t),
                COOLPROP_FLUID,
            )
        except Exception as exc:
            raise RuntimeError(
                f"CoolProp failed at row {i} for "
                f"{COOLPROP_FLUID}: P={p:.6g} Pa, T={t:.6g} K"
            ) from exc

    return values


def real_gas_sonic_state(p0_pa, t0_k):
    """
    Find the sonic state on an isentrope starting from upstream
    stagnation conditions P0,T0.

    At the sonic point:

        V = a

    and for an adiabatic isentropic nozzle:

        h0 = h* + a*^2 / 2

    Therefore solve:

        2*(h0 - h*) - a*^2 = 0

    Returns:
        p_star_pa
        rho_star_kg_m3
        a_star_m_s
        ideal_choked_mass_flux_kg_m2_s
    """
    h0 = PropsSI(
        "Hmass", "P", p0_pa, "T", t0_k, COOLPROP_FLUID
    )
    s0 = PropsSI(
        "Smass", "P", p0_pa, "T", t0_k, COOLPROP_FLUID
    )

    def sonic_function(p_pa):
        h = PropsSI(
            "Hmass", "P", p_pa, "Smass", s0, COOLPROP_FLUID
        )
        a = PropsSI(
            "speed_of_sound",
            "P",
            p_pa,
            "Smass",
            s0,
            COOLPROP_FLUID,
        )

        return 2.0 * (h0 - h) - a**2

    # Search for a sign change over a wide pressure range.
    #
    # The upper end is just below P0 because V=0 at exactly P0.
    p_min = max(100.0, 1.0e-4 * p0_pa)
    p_max = 0.999999 * p0_pa

    pressure_grid = np.geomspace(
        p_min,
        p_max,
        120,
    )

    f_grid = np.full_like(pressure_grid, np.nan, dtype=float)

    for i, p in enumerate(pressure_grid):
        try:
            f_grid[i] = sonic_function(float(p))
        except Exception:
            # A state may be unavailable for some fluids/conditions.
            # Leave NaN and continue searching.
            pass

    p_star = None

    for i in range(len(pressure_grid) - 1):
        f1 = f_grid[i]
        f2 = f_grid[i + 1]

        if not (np.isfinite(f1) and np.isfinite(f2)):
            continue

        if f1 == 0.0:
            p_star = pressure_grid[i]
            break

        if f1 * f2 < 0.0:
            p_star = brentq(
                sonic_function,
                pressure_grid[i],
                pressure_grid[i + 1],
                xtol=1.0e-6,
                rtol=1.0e-10,
                maxiter=100,
            )
            break

    if p_star is None:
        raise RuntimeError(
            "Could not locate a real-gas sonic state with CoolProp "
            f"for P0={p0_pa:.6g} Pa, T0={t0_k:.6g} K."
        )

    rho_star = PropsSI(
        "Dmass", "P", p_star, "Smass", s0, COOLPROP_FLUID
    )
    a_star = PropsSI(
        "speed_of_sound",
        "P",
        p_star,
        "Smass",
        s0,
        COOLPROP_FLUID,
    )

    mass_flux_star = rho_star * a_star

    return p_star, rho_star, a_star, mass_flux_star


def real_gas_ideal_mass_flux(p0_pa, t0_k, p2_pa):
    """
    Ideal isentropic mass flux through the orifice using real-gas
    properties from CoolProp.

    If downstream pressure is below the sonic critical pressure,
    the ideal flow is choked and the sonic mass flux is used.

    If downstream pressure is above the critical pressure,
    the isentropic exit state at P2 is used.

    Returns:
        ideal_mass_flux [kg/(m^2*s)]
        p_star [Pa]
        is_choked [bool]
    """
    if not (
        np.isfinite(p0_pa)
        and np.isfinite(t0_k)
        and np.isfinite(p2_pa)
    ):
        return np.nan, np.nan, False

    if p0_pa <= 0.0 or t0_k <= 0.0:
        return np.nan, np.nan, False

    if p2_pa >= p0_pa:
        return 0.0, np.nan, False

    h0 = PropsSI(
        "Hmass", "P", p0_pa, "T", t0_k, COOLPROP_FLUID
    )
    s0 = PropsSI(
        "Smass", "P", p0_pa, "T", t0_k, COOLPROP_FLUID
    )

    p_star, _, _, choked_mass_flux = real_gas_sonic_state(
        p0_pa,
        t0_k,
    )

    if p2_pa <= p_star:
        return choked_mass_flux, p_star, True

    # Subcritical / unchoked condition.
    h2s = PropsSI(
        "Hmass", "P", p2_pa, "Smass", s0, COOLPROP_FLUID
    )
    rho2s = PropsSI(
        "Dmass", "P", p2_pa, "Smass", s0, COOLPROP_FLUID
    )

    delta_h = h0 - h2s

    if delta_h <= 0.0:
        return 0.0, p_star, False

    velocity_ideal = math.sqrt(2.0 * delta_h)
    mass_flux = rho2s * velocity_ideal

    return mass_flux, p_star, False


def ideal_gas_mass_flux(p0_pa, t0_k, p2_pa, gamma):
    """
    Ideal-gas compressible orifice/nozzle mass flux.

    Uses the user-entered specific gas constant and local gamma.

    Returns:
        ideal_mass_flux [kg/(m^2*s)]
        critical_pressure_ratio [-]
        is_choked [bool]
    """
    if (
        p0_pa <= 0.0
        or t0_k <= 0.0
        or p2_pa >= p0_pa
        or gamma <= 1.0
    ):
        return 0.0, np.nan, False

    pressure_ratio = p2_pa / p0_pa

    critical_ratio = (
        2.0 / (gamma + 1.0)
    ) ** (
        gamma / (gamma - 1.0)
    )

    if pressure_ratio <= critical_ratio:
        mass_flux = (
            p0_pa
            / math.sqrt(t0_k)
            * math.sqrt(gamma / GAS_CONSTANT_J_KG_K)
            * (
                2.0 / (gamma + 1.0)
            ) ** (
                (gamma + 1.0)
                / (2.0 * (gamma - 1.0))
            )
        )

        return mass_flux, critical_ratio, True

    mass_flux = (
        p0_pa
        / math.sqrt(t0_k)
        * math.sqrt(
            (2.0 * gamma)
            / (
                GAS_CONSTANT_J_KG_K
                * (gamma - 1.0)
            )
            * (
                pressure_ratio ** (2.0 / gamma)
                - pressure_ratio
                ** ((gamma + 1.0) / gamma)
            )
        )
    )

    return mass_flux, critical_ratio, False


# ============================================================
# LOAD CSV
# ============================================================

input_path = Path(INPUT_CSV)

if not input_path.exists():
    raise FileNotFoundError(
        f"Input CSV not found: {input_path.resolve()}"
    )

df = pd.read_csv(input_path)

missing = [
    column
    for column in REQUIRED_COLUMNS
    if column not in df.columns
]

if missing:
    raise ValueError(
        "Input CSV is missing required columns:\n"
        + "\n".join(f"  - {column}" for column in missing)
    )

# Force required columns to numeric.
for column in REQUIRED_COLUMNS:
    df[column] = pd.to_numeric(
        df[column],
        errors="coerce",
    )

if df[REQUIRED_COLUMNS].isna().any().any():
    bad_rows = df.index[
        df[REQUIRED_COLUMNS].isna().any(axis=1)
    ].tolist()

    raise ValueError(
        "Required columns contain blank/non-numeric values "
        f"at rows: {bad_rows[:20]}"
    )


# ============================================================
# TIME
# ============================================================

time_s = df["time_ms"].to_numpy(dtype=float) / 1000.0

if len(time_s) < 2:
    raise ValueError(
        "At least two time points are required."
    )

if np.any(np.diff(time_s) <= 0.0):
    raise ValueError(
        "time_ms must be strictly increasing."
    )

df["time_s"] = time_s


# ============================================================
# PRESSURES
# ============================================================

p_supply_abs_pa = pressure_psi_to_absolute_pa(
    df["supply_pressure_psi"]
)

p_upstream_abs_pa = pressure_psi_to_absolute_pa(
    df["upstream_pressure_psi"]
)

p_downstream_abs_pa = pressure_psi_to_absolute_pa(
    df["downstream_pressure_psi"]
)

df["supply_pressure_abs_pa"] = p_supply_abs_pa
df["upstream_pressure_abs_pa"] = p_upstream_abs_pa
df["downstream_pressure_abs_pa"] = p_downstream_abs_pa

df["supply_pressure_psia"] = (
    p_supply_abs_pa / PSI_TO_PA
)
df["upstream_pressure_psia"] = (
    p_upstream_abs_pa / PSI_TO_PA
)
df["downstream_pressure_psia"] = (
    p_downstream_abs_pa / PSI_TO_PA
)


# ============================================================
# PRESSURE DROP
# ============================================================

delta_p_pa = (
    p_upstream_abs_pa
    - p_downstream_abs_pa
)

df["pressure_drop_psi"] = (
    delta_p_pa / PSI_TO_PA
)
df["pressure_drop_pa"] = delta_p_pa


# ============================================================
# REAL-GAS DENSITIES FROM COOLPROP
# ============================================================

t_supply_k = df[
    "supply_temperature_K"
].to_numpy(dtype=float)

t_upstream_k = df[
    "upstream_temperature_K"
].to_numpy(dtype=float)

rho_supply = coolprop_property(
    "Dmass",
    p_supply_abs_pa,
    t_supply_k,
)

rho_upstream = coolprop_property(
    "Dmass",
    p_upstream_abs_pa,
    t_upstream_k,
)

df["supply_density_kg_m3"] = rho_supply
df["upstream_density_kg_m3"] = rho_upstream


# ============================================================
# UPSTREAM Cp, Cv, AND GAMMA
# ============================================================

cp_upstream = coolprop_property(
    "Cpmass",
    p_upstream_abs_pa,
    t_upstream_k,
)

cv_upstream = coolprop_property(
    "Cvmass",
    p_upstream_abs_pa,
    t_upstream_k,
)

gamma_upstream = cp_upstream / cv_upstream

df["upstream_cp_J_kgK"] = cp_upstream
df["upstream_cv_J_kgK"] = cv_upstream
df["upstream_gamma"] = gamma_upstream


# ============================================================
# SUPPLY-TANK MASS
# ============================================================

tank_mass_raw = (
    rho_supply * TANK_VOLUME_M3
)

df["tank_mass_kg_raw"] = tank_mass_raw

if SMOOTHING_WINDOW_ROWS < 1:
    raise ValueError(
        "SMOOTHING_WINDOW_ROWS must be >= 1."
    )

if SMOOTHING_WINDOW_ROWS == 1:
    tank_mass_used = tank_mass_raw.copy()
else:
    tank_mass_used = (
        pd.Series(tank_mass_raw)
        .rolling(
            window=SMOOTHING_WINDOW_ROWS,
            center=True,
            min_periods=1,
        )
        .mean()
        .to_numpy()
    )

df["tank_mass_kg_used"] = tank_mass_used


# ============================================================
# BLOWDOWN MASS-FLOW RATE
# ============================================================

# Mass leaving the tank:
#
#     mdot = -dm_tank/dt
#
# np.gradient provides centered differences internally
# and one-sided differences at the ends.

dm_dt = np.gradient(
    tank_mass_used,
    time_s,
)

mass_flow_kg_s = -dm_dt

df["mass_flow_kg_s"] = mass_flow_kg_s


# ============================================================
# REAL-GAS COMPRESSIBLE IDEAL MASS FLUX
# ============================================================

real_mass_flux = np.full(len(df), np.nan)
real_p_star_pa = np.full(len(df), np.nan)
real_choked = np.zeros(len(df), dtype=bool)

for i in range(len(df)):
    try:
        (
            real_mass_flux[i],
            real_p_star_pa[i],
            real_choked[i],
        ) = real_gas_ideal_mass_flux(
            float(p_upstream_abs_pa[i]),
            float(t_upstream_k[i]),
            float(p_downstream_abs_pa[i]),
        )
    except Exception as exc:
        raise RuntimeError(
            f"Real-gas orifice calculation failed at row {i}."
        ) from exc

real_ideal_mdot = (
    ORIFICE_AREA_M2
    * real_mass_flux
)

df["real_gas_critical_pressure_psia"] = (
    real_p_star_pa / PSI_TO_PA
)
df["real_gas_choked"] = real_choked
df["real_gas_ideal_mass_flux_kg_m2_s"] = real_mass_flux
df["real_gas_ideal_mass_flow_kg_s"] = real_ideal_mdot


# ============================================================
# PRIMARY DISCHARGE COEFFICIENT: REAL-GAS COMPRESSIBLE
# ============================================================

cd_real = np.full(len(df), np.nan)

valid_real = (
    np.isfinite(mass_flow_kg_s)
    & np.isfinite(real_ideal_mdot)
    & (mass_flow_kg_s > 0.0)
    & (real_ideal_mdot > 0.0)
)

cd_real[valid_real] = (
    mass_flow_kg_s[valid_real]
    / real_ideal_mdot[valid_real]
)

df["Cd_real_gas"] = cd_real


# ============================================================
# IDEAL-GAS COMPRESSIBLE COMPARISON
# ============================================================

ideal_gas_flux = np.full(len(df), np.nan)
ideal_critical_ratio = np.full(len(df), np.nan)
ideal_gas_choked = np.zeros(len(df), dtype=bool)

for i in range(len(df)):
    (
        ideal_gas_flux[i],
        ideal_critical_ratio[i],
        ideal_gas_choked[i],
    ) = ideal_gas_mass_flux(
        float(p_upstream_abs_pa[i]),
        float(t_upstream_k[i]),
        float(p_downstream_abs_pa[i]),
        float(gamma_upstream[i]),
    )

ideal_gas_mdot = (
    ORIFICE_AREA_M2
    * ideal_gas_flux
)

cd_ideal = np.full(len(df), np.nan)

valid_ideal = (
    np.isfinite(mass_flow_kg_s)
    & np.isfinite(ideal_gas_mdot)
    & (mass_flow_kg_s > 0.0)
    & (ideal_gas_mdot > 0.0)
)

cd_ideal[valid_ideal] = (
    mass_flow_kg_s[valid_ideal]
    / ideal_gas_mdot[valid_ideal]
)

df["ideal_gas_critical_pressure_ratio"] = (
    ideal_critical_ratio
)
df["ideal_gas_choked"] = ideal_gas_choked
df["ideal_gas_ideal_mass_flux_kg_m2_s"] = ideal_gas_flux
df["ideal_gas_ideal_mass_flow_kg_s"] = ideal_gas_mdot
df["Cd_ideal_gas"] = cd_ideal


# ============================================================
# INCOMPRESSIBLE-EQUIVALENT Cd
# ============================================================

# This is NOT the preferred Cd for a sonic/choked gas orifice.
# It is included because it directly uses:
#
#     mdot = Cd * A * sqrt(2*rho*DeltaP)
#
# and can be useful as a diagnostic.

incompressible_ideal_mdot = np.full(len(df), np.nan)

valid_dp = (
    np.isfinite(delta_p_pa)
    & np.isfinite(rho_upstream)
    & (delta_p_pa > 0.0)
    & (rho_upstream > 0.0)
)

incompressible_ideal_mdot[valid_dp] = (
    ORIFICE_AREA_M2
    * np.sqrt(
        2.0
        * rho_upstream[valid_dp]
        * delta_p_pa[valid_dp]
    )
)

cd_incompressible = np.full(len(df), np.nan)

valid_incompressible = (
    valid_dp
    & np.isfinite(mass_flow_kg_s)
    & (mass_flow_kg_s > 0.0)
    & (incompressible_ideal_mdot > 0.0)
)

cd_incompressible[valid_incompressible] = (
    mass_flow_kg_s[valid_incompressible]
    / incompressible_ideal_mdot[valid_incompressible]
)

df["incompressible_ideal_mass_flow_kg_s"] = (
    incompressible_ideal_mdot
)
df["Cd_incompressible_equiv"] = cd_incompressible


# ============================================================
# CONSTANT GEOMETRY OUTPUT
# ============================================================

df["orifice_diameter_in"] = ORIFICE_DIAMETER_IN
df["orifice_diameter_m"] = ORIFICE_DIAMETER_M
df["orifice_area_m2"] = ORIFICE_AREA_M2
df["tank_volume_L"] = TANK_VOLUME_L
df["tank_volume_m3"] = TANK_VOLUME_M3


# ============================================================
# AVERAGES
# ============================================================

valid_mdot = (
    np.isfinite(mass_flow_kg_s)
    & (mass_flow_kg_s > 0.0)
)

valid_cd_real = (
    np.isfinite(cd_real)
    & (cd_real > 0.0)
)

average_mass_flow = (
    np.mean(mass_flow_kg_s[valid_mdot])
    if np.any(valid_mdot)
    else np.nan
)

average_cd_real = (
    np.mean(cd_real[valid_cd_real])
    if np.any(valid_cd_real)
    else np.nan
)

average_cd_ideal = (
    np.nanmean(cd_ideal)
    if np.any(np.isfinite(cd_ideal))
    else np.nan
)

average_cd_incompressible = (
    np.nanmean(cd_incompressible)
    if np.any(np.isfinite(cd_incompressible))
    else np.nan
)

# A useful independent check:
#
# average mdot over full run =
# total tank mass loss / total elapsed time

total_time_s = time_s[-1] - time_s[0]

overall_blowdown_mass_flow = (
    (tank_mass_used[0] - tank_mass_used[-1])
    / total_time_s
)


# ============================================================
# WRITE OUTPUT CSV
# ============================================================

output_path = Path(OUTPUT_CSV)

df.to_csv(
    output_path,
    index=False,
)


# ============================================================
# TERMINAL OUTPUT
# ============================================================

print()
print("=" * 72)
print("              BLOWDOWN SONIC-ORIFICE Cd CALIBRATION")
print("=" * 72)

print(f"CoolProp fluid              : {COOLPROP_FLUID}")
print(f"Orifice diameter            : {ORIFICE_DIAMETER_IN:.6f} in")
print(f"Orifice area                : {ORIFICE_AREA_M2:.8e} m^2")
print(f"Tank volume                 : {TANK_VOLUME_L:.6f} L")
print(f"Specific gas constant       : {GAS_CONSTANT_J_KG_K:.6f} J/(kg*K)")
print(
    "CSV pressures interpreted as: "
    + ("psig" if PRESSURES_ARE_GAUGE else "psia")
)

print()
print(f"Average mass flow rate      : {average_mass_flow:.8f} kg/s")
print(
    f"Full-run blowdown avg mdot  : "
    f"{overall_blowdown_mass_flow:.8f} kg/s"
)

print()
print(
    f"Average Cd (REAL GAS)       : "
    f"{average_cd_real:.6f}"
)
print(
    f"Average Cd (ideal gas)      : "
    f"{average_cd_ideal:.6f}"
)
print(
    f"Average Cd (incompressible) : "
    f"{average_cd_incompressible:.6f}"
)

print()
print(
    f"Real-gas choked rows        : "
    f"{np.count_nonzero(real_choked)} / {len(df)}"
)
print(
    f"Valid real-gas Cd rows      : "
    f"{np.count_nonzero(valid_cd_real)} / {len(df)}"
)

print()
print(f"Output CSV                  : {output_path.resolve()}")
print("=" * 72)
print()
