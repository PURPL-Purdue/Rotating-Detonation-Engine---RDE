import numpy as np
import math
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI

# -------------------------------------------------------------------------
# INITIALIZATION
# -------------------------------------------------------------------------
T = 293.15                  # K
psi_to_Pa = 6894.76         # psi to pa conversion
in_to_m = 0.0254            # inches to meters conversion
R_univ = 8.314462           # J/mol-K
P_atm = 101325.0            # Pa
Cd = 0.61                   # discharge coefficient

# -------------------------------------------------------------------------
# MIXTURE COMPOSITION
# -------------------------------------------------------------------------
Y_O2 = 0.88889 # mass fraction of o2
Y_H2 = 0.11111 # mass fraction of h2

MW_O2 = 0.032 # molar weight of o2 (g/mol)
MW_H2 = 0.002016 # molar weight of h2 (g/mol)

# Convert mass fractions → mole fractions
nO2 = Y_O2 / MW_O2 # moles of o2
nH2 = Y_H2 / MW_H2 # moles of h2
n_sum = nO2 + nH2 # total moles

xO2 = nO2 / n_sum # mole fraction of o2
xH2 = nH2 / n_sum # mole fraction of h2

MW_mix = xO2*MW_O2 + xH2*MW_H2 # molar weight of mixture
R_mix = R_univ / MW_mix # gas constant of mixture

# -------------------------------------------------------------------------
# INITIAL FEED PRESSURES
# -------------------------------------------------------------------------
p1_O2 = 400 * psi_to_Pa # feed pressure of o2
p1_H2 = 200 * psi_to_Pa # feed pressure of h2

# Assume chamber initially at choked pressure (isentropic)
cp_GOx = PropsSI('CPMASS', 'T', T, 'P', p1_O2, 'Oxygen')      # J/kg-K
cp_GH2 = PropsSI('CPMASS', 'T', T, 'P', p1_H2, 'Hydrogen')    # J/kg-K

cv_GOx = PropsSI('CVMASS', 'T', T, 'P', p1_O2, 'Oxygen')      # J/kg-K
cv_GH2 = PropsSI('CVMASS', 'T', T, 'P', p1_H2, 'Hydrogen')    # J/kg-K

g_GOx = cp_GOx / cv_GOx # gamma of o2
g_GH2 = cp_GH2 / cv_GH2 # gamma of h2

P0_O2 = p1_O2 * (2 / (g_GOx + 1))**(g_GOx / (g_GOx - 1)) # choked pressure of o2
P0_H2 = p1_H2 * (2 / (g_GH2 + 1))**(g_GH2 / (g_GH2 - 1)) # choked pressure of h2
P0 = xO2*P0_O2 + xH2*P0_H2 # total choked pressure

# -------------------------------------------------------------------------
# THERMODYNAMIC PROPERTIES (evaluated at chamber conditions)
# -------------------------------------------------------------------------
Cp_O2 = PropsSI('CPMASS','T',T,'P', P0_O2 ,'Oxygen') # new cp of o2
Cp_H2 = PropsSI('CPMASS','T',T,'P', P0_H2 ,'Hydrogen') # new cp of h2

# Convert to molar, mix, convert back
Cp_O2_m = Cp_O2 * MW_O2 # molar cp of o2
Cp_H2_m = Cp_H2 * MW_H2 # molar cp of h2
Cp_mix_m = xO2*Cp_O2_m + xH2*Cp_H2_m # molar cp of mixture
Cp_mix = Cp_mix_m / MW_mix # cp of mixture

Cv_mix = Cp_mix - R_mix # cv of mixture
gamma = Cp_mix / Cv_mix # gamma of mixture

# Viscosity from CoolProp
mu_O2 = PropsSI('VISCOSITY','T',T,'P',P0_O2,'Oxygen') # viscosity of o2
mu_H2 = PropsSI('VISCOSITY','T',T,'P',P0_H2,'Hydrogen') # viscosity of h2

# Wilke mixing rule
def phi(mu_i, mu_j, MW_i, MW_j):
    return (1 + np.sqrt(mu_i/mu_j)*(MW_j/MW_i)**0.25)**2 / \
           np.sqrt(8*(1 + MW_i/MW_j))

phi_O2_H2 = phi(mu_O2, mu_H2, MW_O2, MW_H2)
phi_H2_O2 = phi(mu_H2, mu_O2, MW_H2, MW_O2)

mu_mix = (
    xO2*mu_O2/(xO2 + xH2*phi_O2_H2) +
    xH2*mu_H2/(xH2 + xO2*phi_H2_O2)
)

# -------------------------------------------------------------------------
# GEOMETRY
# -------------------------------------------------------------------------
thick1 = 0.049 * in_to_m # thickness of tube 1
thick2 = 0.065 * in_to_m # thickness of tube 2
do1 = 0.25 * in_to_m # outer diameter of tube 1
do2 = 0.375 * in_to_m # outer diameter of tube 2

di1 = do1 - 2*thick1 # inner diameter of tube 1
di2 = do2 - 2*thick2 # inner diameter of tube 2

a1 = math.pi*(di1/2)**2 # inner area of tube 1
a2 = math.pi*(di2/2)**2 # inner area of tube 2

l1 = (4.3955 + 0.74)*in_to_m # length of tube 1
l2 = 1.5*in_to_m # length of tube 2

Vt = 0.1407 * in_to_m**3 # volume of pre-detonator

# -------------------------------------------------------------------------
# TIME SETUP
# -------------------------------------------------------------------------
dt = 1e-6 # time step
t_final = 0.1 # total time
steps = int(t_final/dt) # number of steps

time = np.zeros(steps) # time vector
P_tube = np.zeros(steps) # pressure vector

# Initial mass
P = P0 # initial pressure
m = P*Vt/(R_mix*T) # mass in volume
unChokeTime = np.nan # time when flow unchokes
lastIndex = steps

# Critical pressure ratio
critical_ratio = (2/(gamma+1))**(gamma/(gamma-1)) # critical pressure ratio

# -------------------------------------------------------------------------
# MAIN LOOP
# -------------------------------------------------------------------------
for i in range(steps):

    rho = m/Vt

    if P_atm/P <= critical_ratio:
        # Choked
        m_dot = Cd*a1*P*np.sqrt(gamma/(R_mix*T)) * \
                (2/(gamma+1))**((gamma+1)/(2*(gamma-1)))
    else:
        # Unchoked (compressible)
        if np.isnan(unChokeTime):
            unChokeTime = i*dt

        pr = P_atm/P
        m_dot = Cd*a1*P*np.sqrt(
            (2*gamma/(R_mix*T*(gamma-1))) *
            (pr**(2/gamma) - pr**((gamma+1)/gamma))
        )

    m -= m_dot*dt
    m = max(m,0)

    P = m*R_mix*T/Vt

    time[i] = (i+1)*dt
    P_tube[i] = P

    if P <= P_atm:
        lastIndex = i+1
        break

# Trim arrays
time = time[:lastIndex]
P_tube = P_tube[:lastIndex]

# -------------------------------------------------------------------------
# PLOT
# -------------------------------------------------------------------------
plt.figure()
plt.plot(time*1e3, P_tube/psi_to_Pa, linewidth=2)
plt.xlabel("Time (ms)")
plt.ylabel("Pressure (psia)")
plt.title("Pre-Detonator Pressure Decay")
plt.grid(True)

if not np.isnan(unChokeTime):
    plt.axvline(unChokeTime*1e3, linestyle="--", label="Unchoked")

plt.axhline(14.7, linestyle="--", label="Atmospheric")
plt.legend()
plt.show()