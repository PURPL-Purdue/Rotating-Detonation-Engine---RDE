import pandas as pd
from pathlib import Path
import CoolProp.CoolProp as cp
import math
import numpy as np
import matplotlib.pyplot as plt

#------------------------------------ SETUP ----------------------------------------#

# Folder containing this Python script
script_dir = Path(__file__).resolve().parent

# CSV in the same folder as the Python script
filename = script_dir / "orifice_blowdown_input.csv"

# Read CSV
data = pd.read_csv(filename)

# Conversion Factors
psi_to_Pa = 6894.75729
L_to_m3 = 0.001
in_to_m = 0.0254
P_atm_psi = 14.6959

#------------------------------------ DATA PRE-PROCESSING ----------------------------------------#

# Extract Data
time_ms = data["time_ms"].to_numpy()

supply_pressure_psi = data["supply_pressure_psi"].to_numpy()
upstream_pressure_psi = data["upstream_pressure_psi"].to_numpy()
downstream_pressure_psi = data["downstream_pressure_psi"].to_numpy()

supply_temperature_K = data["supply_temperature_K"].to_numpy()
upstream_temperature_K = data["upstream_temperature_K"].to_numpy()
downstream_temperature_K = data["downstream_temperature_K"].to_numpy()

# Convert Data
supply_pressure = (supply_pressure_psi + P_atm_psi) * psi_to_Pa
upstream_pressure = (upstream_pressure_psi + P_atm_psi) * psi_to_Pa
downstream_pressure = (downstream_pressure_psi + P_atm_psi) * psi_to_Pa

time_s = time_ms / 1000

#------------------------------------ INPUTS ----------------------------------------#

gas = "Air"
tank_volume = 50 # L
orifice_diameter = 0.24 # in

#------------------------------------ CALCULATIONS ----------------------------------------#

# Gas Properties
R_u = cp.PropsSI("GAS_CONSTANT", gas)   # J/(mol*K)
M = cp.PropsSI("MOLAR_MASS", gas)       # kg/mol
R = R_u / M                             # J/(kg*K)

# Tank Dimensions
tank_volume = tank_volume * L_to_m3

# Orifice Dimensions
orifice_radius = (orifice_diameter * in_to_m) / 2
orifice_area = math.pi * orifice_radius**2

# Tank Mass Calculations
tank_masses = []

i = 0
while i < len(supply_pressure):

    density = cp.PropsSI(
        "Dmass",
        "P", supply_pressure[i],
        "T", supply_temperature_K[i],
        gas
    )

    mass = density * tank_volume
    tank_masses.append(mass)

    i += 1

tank_masses = np.array(tank_masses)

# Actual Cumulative Mass Loss
mass_actual_cumulative = tank_masses[0] - tank_masses

# Ideal Mass Flow Calculations
m_dot_ideals = []
gammas = []

i = 0
while i < len(supply_pressure):

    # Flow Conditions
    P_0 = supply_pressure[i]
    T_0 = supply_temperature_K[i]

    # Gas Properties
    c_p = cp.PropsSI(
        "Cpmass",
        "P", P_0,
        "T", T_0,
        gas
    )

    c_v = cp.PropsSI(
        "Cvmass",
        "P", P_0,
        "T", T_0,
        gas
    )

    gamma = c_p / c_v
    gammas.append(gamma)

    # Ideal Choked Mass Flow for Cd = 1
    m_dot_ideal = (
        orifice_area
        * P_0
        / math.sqrt(T_0)
        * math.sqrt(gamma / R)
        * (2 / (gamma + 1))**(
            (gamma + 1) / (2 * (gamma - 1))
        )
    )

    m_dot_ideals.append(m_dot_ideal)

    i += 1

m_dot_ideals = np.array(m_dot_ideals)
gammas = np.array(gammas)

# Ideal Cumulative Mass Loss for Cd = 1
mass_ideal_cumulative = np.zeros(len(time_s))

i = 1
while i < len(time_s):

    dt = time_s[i] - time_s[i - 1]

    mass_ideal_cumulative[i] = (
        mass_ideal_cumulative[i - 1]
        + 0.5
        * (m_dot_ideals[i - 1] + m_dot_ideals[i])
        * dt
    )

    i += 1

# Constant Discharge Coefficient Using Least Squares
x = mass_ideal_cumulative[1:]
y = mass_actual_cumulative[1:]

Cd = np.sum(x * y) / np.sum(x**2)

# Integrated Discharge Coefficient
mass_actual_total = mass_actual_cumulative[-1]
mass_ideal_total = mass_ideal_cumulative[-1]

Cd_integrated = mass_actual_total / mass_ideal_total

# Predicted Cumulative Mass Loss
mass_predicted_cumulative = Cd * mass_ideal_cumulative

# Predicted Mass Flow
m_dot_predicted = Cd * m_dot_ideals

#------------------------------------ OUTPUT ----------------------------------------#

print(f"Calibrated Cd: {Cd:.4f}")
print(f"Integrated Cd: {Cd_integrated:.4f}")
print(f"Actual Mass Loss: {mass_actual_total:.6f} kg")
print(f"Ideal Mass Loss (Cd = 1): {mass_ideal_total:.6f} kg")

plt.figure(1, figsize=(10, 6))

plt.plot(
    time_ms,
    mass_actual_cumulative,
    linewidth=2,
    label="Measured Mass Loss"
)

plt.plot(
    time_ms,
    mass_predicted_cumulative,
    linestyle="--",
    linewidth=2,
    label=f"Model, $C_d$ = {Cd:.3f}"
)

plt.xlabel("Time (ms)", fontsize=16)
plt.ylabel("Cumulative Mass Loss (kg)", fontsize=16)
plt.title("Orifice Blowdown Calibration", fontsize=18)

plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.grid(True, alpha=0.3)
plt.legend(fontsize=13)

plt.tight_layout()
plt.show()