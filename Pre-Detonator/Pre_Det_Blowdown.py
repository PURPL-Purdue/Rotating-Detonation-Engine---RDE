import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------------------------------
# INITIALIZATION
# -------------------------------------------------------------------------
T = 300.0  # K

# -------------------------------------------------------------------------
# MIXTURE COMPOSITION AND THERMO
# -------------------------------------------------------------------------
# Mass fractions
Y_O2 = 0.88889
Y_H2 = 0.11111

# Molecular weights (kg/mol)
MW_O2 = 0.032
MW_H2 = 0.002016

# Universal gas constant
R_univ = 8.314462  # J/(mol*K)

# Species Cp (J/kg-K)
Cp_O2 = 0.918e3
Cp_H2 = 14.3e3

# --- Mass fractions -> mole fractions
nO2 = Y_O2 / MW_O2
nH2 = Y_H2 / MW_H2
n_sum = nO2 + nH2

xO2 = nO2 / n_sum
xH2 = nH2 / n_sum

# --- Mixture properties
MW_mix = Y_O2 * MW_O2 + Y_H2 * MW_H2
R_mix = R_univ / MW_mix
Cp_mix = Y_O2 * Cp_O2 + Y_H2 * Cp_H2
Cv_mix = Cp_mix - R_mix
gamma = Cp_mix / Cv_mix

# -------------------------------------------------------------------------
# MIXTURE VISCOSITY (Wilke)
# -------------------------------------------------------------------------
mu_O2 = 2.07e-5  # Pa*s
mu_H2 = 8.76e-6  # Pa*s

def phi(mu_i, mu_j, MW_i, MW_j):
    return (1 + np.sqrt(mu_i / mu_j) * (MW_j / MW_i)**0.25)**2 / \
           np.sqrt(8 * (1 + MW_i / MW_j))

phi_O2_H2 = phi(mu_O2, mu_H2, MW_O2, MW_H2)
phi_H2_O2 = phi(mu_H2, mu_O2, MW_H2, MW_O2)

mu = (
    xO2 * mu_O2 / (xO2 + xH2 * phi_O2_H2) +
    xH2 * mu_H2 / (xH2 + xO2 * phi_H2_O2)
)

# -------------------------------------------------------------------------
# GEOMETRY
# -------------------------------------------------------------------------
L = 0.1524       # m
D = 0.004572     # m
A = np.pi * (D / 2)**2
Vt = A * L

# -------------------------------------------------------------------------
# PRESSURE CONDITIONS
# -------------------------------------------------------------------------
P0 = (110 + 14.7) * 6894.76  # Pa
P_atm = 101325.0            # Pa

# -------------------------------------------------------------------------
# TIME SETTINGS
# -------------------------------------------------------------------------
dt = 1e-6
t_final = 0.005
steps = int(np.floor(t_final / dt))

time = np.zeros(steps)
P_tube = np.zeros(steps)

# -------------------------------------------------------------------------
# INITIAL CONDITIONS
# -------------------------------------------------------------------------
P = P0
m = P * Vt / (R_mix * T)
unChokeTime = np.nan

# -------------------------------------------------------------------------
# MAIN LOOP
# -------------------------------------------------------------------------
for i in range(steps):

    rho = m / Vt

    # Initial choked mass flow (bootstrap)
    if i == 0:
        m_dot_prev = (
            A * P * np.sqrt(gamma / (R_mix * T)) *
            (2 / (gamma + 1))**((gamma + 1) / (2 * (gamma - 1)))
        )

    # Velocity and Reynolds number
    V = m_dot_prev / (rho * A + 1e-12)
    Re = abs(rho * V * D / mu)

    # Discharge coefficient
    Cd = 0.5959 + 0.0312 / np.sqrt(max(Re, 1e-6))

    # Choking condition
    P_crit = P * (2 / (gamma + 1))**(gamma / (gamma - 1))

    if P_crit > P_atm:
        # Choked
        m_dot = (
            Cd * A * P * np.sqrt(gamma / (R_mix * T)) *
            (2 / (gamma + 1))**((gamma + 1) / (2 * (gamma - 1)))
        )
    else:
        # Unchoked
        if np.isnan(unChokeTime):
            unChokeTime = i * dt
        m_dot = Cd * A * np.sqrt(2 * rho * (P - P_atm))

    m_dot_prev = m_dot

    # Update mass and pressure
    m -= m_dot * dt
    m = max(m, 0.0)
    P = m * R_mix * T / Vt

    time[i] = (i + 1) * dt
    P_tube[i] = P

    if P <= P_atm:
        lastIndex = i + 1
        break

# -------------------------------------------------------------------------
# TRIM ARRAYS
# -------------------------------------------------------------------------
time = time[:lastIndex]
P_tube = P_tube[:lastIndex]

# -------------------------------------------------------------------------
# PLOT (pressure in psia)
# -------------------------------------------------------------------------
plt.figure()
plt.plot(time * 1e3, P_tube / 6894.76, linewidth=2)
plt.xlabel("Time (ms)")
plt.ylabel("Pressure (psia)")
plt.title("Tube Pressure Decay – Choked → Unchoked")
plt.grid(True)

if not np.isnan(unChokeTime):
    plt.axvline(unChokeTime * 1e3, linestyle="--", linewidth=1.5, label="Unchoked")

plt.axhline(14.7, linestyle="--", linewidth=1, label="Atmospheric")
plt.legend()
plt.show()
