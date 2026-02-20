import numpy as np
import math
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

# -------------------------------------------------------------------------
# INITIALIZATION
# -------------------------------------------------------------------------
T = 293.15  # K
in_to_m = 0.0254
psi_to_Pa = 6894.76

# -------------------------------------------------------------------------
# MIXTURE COMPOSITION
# -------------------------------------------------------------------------
Y_O2 = 0.88889
Y_H2 = 0.11111

MW_O2 = 0.032
MW_H2 = 0.002016
R_univ = 8.314462

p1_O2 = 400 * psi_to_Pa
p1_H2 = 200 * psi_to_Pa

Cp_O2 = PropsSI('CPMASS', 'T', T, 'P', p1_O2, 'Oxygen')
Cp_H2 = PropsSI('CPMASS', 'T', T, 'P', p1_H2, 'Hydrogen')

# Mole fractions
nO2 = Y_O2 / MW_O2
nH2 = Y_H2 / MW_H2
n_sum = nO2 + nH2

xO2 = nO2 / n_sum
xH2 = nH2 / n_sum

MW_mix = xO2 * MW_O2 + xH2 * MW_H2
R_mix = R_univ / MW_mix
Cp_mix = xO2 * Cp_O2 + xH2 * Cp_H2
Cv_mix = Cp_mix - R_mix
gamma = Cp_mix / Cv_mix

# -------------------------------------------------------------------------
# VISCOSITY (Wilke)
# -------------------------------------------------------------------------
mu_O2 = 2.07e-5
mu_H2 = 8.76e-6

def phi(mu_i, mu_j, MW_i, MW_j):
    return (1 + np.sqrt(mu_i/mu_j)*(MW_j/MW_i)**0.25)**2 / \
           np.sqrt(8*(1 + MW_i/MW_j))

phi_O2_H2 = phi(mu_O2, mu_H2, MW_O2, MW_H2)
phi_H2_O2 = phi(mu_H2, mu_O2, MW_H2, MW_O2)

mu = (
    xO2 * mu_O2 / (xO2 + xH2 * phi_O2_H2) +
    xH2 * mu_H2 / (xH2 + xO2 * phi_H2_O2)
)

# -------------------------------------------------------------------------
# GEOMETRY (CONVERTED TO METERS)
# -------------------------------------------------------------------------
thick1 = 0.049 * in_to_m
thick2 = 0.065 * in_to_m
do1 = 0.25 * in_to_m
do2 = 0.375 * in_to_m

di1 = do1 - 2*thick1
di2 = do2 - 2*thick2

a1 = math.pi * (di1/2)**2
a2 = math.pi * (di2/2)**2

l1 = 1.5 * in_to_m
l2 = (4.3955 + 0.74) * in_to_m

spring_volume = 0.0216 * in_to_m**3
manifold_volume = math.pi * ((0.3873*in_to_m)/2)**2 * (1.19 * in_to_m)

tube_volume = (a1*l1) + (a2*l2) - spring_volume
Vt = manifold_volume + tube_volume

# -------------------------------------------------------------------------
# PRESSURES
# -------------------------------------------------------------------------
P0 = (xH2*p1_H2 + xO2*p1_O2)*0.5283
P_atm = 101325.0

# -------------------------------------------------------------------------
# TIME
# -------------------------------------------------------------------------
dt = 1e-6
t_final = 0.1
steps = int(t_final/dt)

time = np.zeros(steps)
P_tube = np.zeros(steps)

# -------------------------------------------------------------------------
# INITIAL CONDITIONS
# -------------------------------------------------------------------------
P = P0
m = P*Vt/(R_mix*T)
lastIndex = steps
unChokeTime = np.nan

# -------------------------------------------------------------------------
# MAIN LOOP
# -------------------------------------------------------------------------
for i in range(steps):

    rho = m/Vt

    # Bootstrap mass flow
    if i == 0:
        m_dot_prev = (
            a1 * P * np.sqrt(gamma/(R_mix*T)) *
            (2/(gamma+1))**((gamma+1)/(2*(gamma-1)))
        )

    V1 = m_dot_prev/(rho*a1 + 1e-12)
    Re = abs(rho*V1*di1/mu)

    Cd = 0.5959 + 0.0312/np.sqrt(max(Re,1e-6))

    P_crit = P*(2/(gamma+1))**(gamma/(gamma-1))

    if P_crit > P_atm:
        # Choked
        m_dot = (
            Cd*a1*P*np.sqrt(gamma/(R_mix*T)) *
            (2/(gamma+1))**((gamma+1)/(2*(gamma-1)))
        )
    else:
        # Unchoked
        if np.isnan(unChokeTime):
            unChokeTime = i*dt
        deltaP = max(P - P_atm, 0)
        m_dot = Cd*a1*np.sqrt(2*rho*deltaP)

    m_dot_prev = m_dot

    m -= m_dot*dt
    m = max(m, 0)

    P = m*R_mix*T/Vt

    time[i] = (i+1)*dt
    P_tube[i] = P

    if P <= P_atm:
        lastIndex = i+1
        break

# -------------------------------------------------------------------------
# TRIM
# -------------------------------------------------------------------------
time = time[:lastIndex]
P_tube = P_tube[:lastIndex]

# -------------------------------------------------------------------------
# PLOT
# -------------------------------------------------------------------------
plt.figure()
plt.plot(time*1e3, P_tube/psi_to_Pa, linewidth=2)
plt.xlabel("Time (ms)")
plt.ylabel("Pressure (psia)")
plt.title("Tube Pressure Decay – Choked → Unchoked")
plt.grid(True)

if not np.isnan(unChokeTime):
    plt.axvline(unChokeTime*1e3, linestyle="--", linewidth=1.5, label="Unchoked")

plt.axhline(14.7, linestyle="--", linewidth=1, label="Atmospheric")
plt.legend()
plt.show()