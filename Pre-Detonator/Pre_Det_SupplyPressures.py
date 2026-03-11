import math as math
import CoolProp.CoolProp as cp
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve

#--------------------------------- CONVERSIONS -------------------------------#
psi_to_Pa = 6894.76
bar_to_psia = 14.5
in_to_m = 0.0254
m_to_ft = 3.28084

#--------------------------------- MAIN CODE ---------------------------------#

# Inputs
p1_GOx = 400 * psi_to_Pa
p1_GH2 = 200 * psi_to_Pa
d_GOx = 0.01 * in_to_m
d_GH2 = 0.01 * in_to_m
dl_GOx = 0.159 * in_to_m
dl_GH2 = 0.159 * in_to_m

# Constants
Rbar = 8.3144               # J/K-mol
MW_GOx = 31.9988/1000       # Oxygen Molecular Weight (kg/mol)
MW_GH2 = 2.01588/1000       # Hydrogen Molecular Weight (kg/mol)
o_f = 8                     # O/F Ratio, reference (phi = 1)
Tt = 10 + 273.15            # Static Temperature
Cd = 0.61                   # Discharge Coefficient

# Mass Flow Calculations

cp_GOx = PropsSI('CPMASS', 'T', Tt, 'P', p1_GOx, 'Oxygen')      # J/kg-K
cp_GH2 = PropsSI('CPMASS', 'T', Tt, 'P', p1_GH2, 'Hydrogen')    # J/kg-K

cv_GOx = PropsSI('CVMASS', 'T', Tt, 'P', p1_GOx, 'Oxygen')      # J/kg-K
cv_GH2 = PropsSI('CVMASS', 'T', Tt, 'P', p1_GH2, 'Hydrogen')    # J/kg-K

g_GOx = cp_GOx / cv_GOx
g_GH2 = cp_GH2 / cv_GH2

cf_GOx = (2 / (g_GOx + 1))**(g_GOx/(g_GOx-1))
cf_GH2 = (2 / (g_GH2 + 1))**(g_GH2/(g_GH2-1))

R_GOx = Rbar / MW_GOx
R_GH2 = Rbar / MW_GH2

A_GOx = math.pi * (d_GOx / 2)**2
A_GH2 = math.pi * (d_GH2 / 2)**2

Tc_GOx = PropsSI('TCRIT', 'T', Tt, 'P', p1_GOx, 'O2')
Tc_GH2 = PropsSI('TCRIT', 'T', Tt, 'P', p1_GH2, 'H2')

md_GOx = Cd * (A_GOx * p1_GOx / math.sqrt(Tt)) * math.sqrt(g_GOx / R_GOx) * ((g_GOx + 1) / 2)**(-1 * ((g_GOx + 1) / (2 * (g_GOx - 1))))
md_GH2 = Cd * (A_GH2 * p1_GH2 / math.sqrt(Tt)) * math.sqrt(g_GH2 / R_GH2) * ((g_GH2 + 1) / 2)**(-1 * ((g_GH2 + 1) / (2 * (g_GH2 - 1))))

v_GOx = math.sqrt(g_GOx * R_GOx * Tc_GOx)
v_GH2 = math.sqrt(g_GH2 * R_GH2 * Tc_GH2)

# Line Velocity Calculations

d_line = 0.250 * in_to_m     # inches
t_line = 0.049 * in_to_m     # inches
flow_area = math.pi * ((d_line - t_line*2)/2)**2    # m^2

rho1_GOx = PropsSI('DMASS', 'T', Tt, 'P', p1_GOx, 'O2')
rho1_GH2 = PropsSI('DMASS', 'T', Tt, 'P', p1_GH2, 'H2')

v_line_GOx = md_GOx / (rho1_GOx * flow_area)         # m/s
v_line_GH2 = md_GH2 / (rho1_GH2 * flow_area)         # m/s

# Fill Time Calculations

p2_GOx = p1_GOx * (2 / (g_GOx + 1))**(g_GOx / (g_GOx - 1))
p2_GH2 = p1_GH2 * (2 / (g_GH2 + 1))**(g_GH2 / (g_GH2 - 1))

rho2_GOx = PropsSI('DMASS', 'T', Tc_GOx, 'P', p2_GOx, 'O2')
rho2_GH2 = PropsSI('DMASS', 'T', Tc_GH2, 'P', p2_GH2, 'H2')

q_GOx = md_GOx / rho2_GOx
q_GH2 = md_GH2 / rho2_GH2
q_tot = q_GOx + q_GH2

vol_t = (0.1407 * in_to_m**3) # man_vol + tube_vol # Total Volume, m^3

t_fill = vol_t / q_tot

#------------------------------ PRINT COMMANDS ------------------------------#
print()
print(f"Pre-Detonator Fill Time: {t_fill * 1000:.3f} ms")
print(f"Mass Ratio: {md_GOx / md_GH2:.6f}")
print(f"O2 Mass Flow Rate: {md_GOx*1000:.3f} g/s")
print(f"H2 Mass Flow Rate: {md_GH2*1000:.3f} g/s")
print(f"O2 Choked Velocity: {v_GOx * m_to_ft:.3f} ft/s")
print(f"O2 Line Velocity: {v_line_GOx * m_to_ft:.3f} ft/s")
print(f"H2 Choked Velocity: {v_GH2 * m_to_ft:.3f} ft/s")
print(f"H2 Line Velocity: {v_line_GH2 * m_to_ft:.3f} ft/s")
print()