import math as math
import CoolProp.CoolProp as cp
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve

#--------------------------------- CONVERSIONS -------------------------------#
psi_to_Pa = 6894.76
bar_to_psia = 14.5
in_to_m = 0.0254
choke_factor = 0.5283

#--------------------------------- MAIN CODE ---------------------------------#

# Inputs
p1_GOx = 300 * psi_to_Pa
p1_GH2 = 150 * psi_to_Pa
d_GOx = 0.01 * in_to_m
d_GH2 = 0.01 * in_to_m

# Constants
Rbar = 8.3144 # J / K*mol
MW_GOx = 31.9988/1000        # Oxygen Molecular Weight (kg/mol)
MW_GH2 = 2.01588/1000        # Hydrogen Molecular Weight (kg/mol)
o_f = 8                     # O/F Ratio, reference (phi = 1)
Tt = 20 + 273.15

# Calculations
cp_GOx = PropsSI('CPMASS', 'T', Tt, 'P', p1_GOx, 'Oxygen') #J/kg-K
cv_GOx = PropsSI('CVMASS', 'T', Tt, 'P', p1_GOx, 'Oxygen') #J/kg-K
g_GOx = cp_GOx / cv_GOx
cp_GH2 = PropsSI('CPMASS', 'T', Tt, 'P', p1_GH2, 'Hydrogen') #J/kg-K
cv_GH2 = PropsSI('CVMASS', 'T', Tt, 'P', p1_GH2, 'Hydrogen') #J/kg-K
g_GH2 = cp_GH2 / cv_GH2
R_GOx = Rbar / MW_GOx
R_GH2 = Rbar / MW_GH2
A_GOx = math.pi * (d_GOx / 2)**2
A_GH2 = math.pi * (d_GH2 / 2)**2
rho_GOx = p1_GOx / (R_GOx * Tt)
rho_GH2 = p1_GH2 / (R_GH2 * Tt)

md_GOx = (A_GOx * p1_GOx / math.sqrt(Tt)) * math.sqrt(g_GOx / R_GOx) * ((g_GOx + 1) / 2)**(-1 * ((g_GOx + 1) / (2 * (g_GOx - 1))))
md_GH2 = (A_GH2 * p1_GH2 / math.sqrt(Tt)) * math.sqrt(g_GH2 / R_GH2) * ((g_GH2 + 1) / 2)**(-1 * ((g_GH2 + 1) / (2 * (g_GH2 - 1))))

v_GOx = md_GOx / (rho_GOx * A_GOx)
v_GH2 = md_GH2 / (rho_GH2 * A_GH2)

#------------------------------ PRINT COMMANDS ------------------------------#
print(f"O2 Mass Flow Rate: {md_GOx*1000:.3f} g/s")
print(f"H2 Mass Flow Rate: {md_GH2*1000:.3f} g/s")
print(f"Mass Ratio: {md_GOx / md_GH2:.6f}")
print(f"O2 Flow Velocity: {v_GOx:.3f} m/s")
print(f"H2 Flow Velocity: {v_GH2:.3f} m/s")