import math
import CoolProp.CoolProp as cp
from scipy.optimize import fsolve

#--------------------------------- CONVERSIONS -------------------------------#
psi_to_Pa = 6894.76
bar_to_psia = 14.5
choke_factor = 0.5283

#--------------------------------- MAIN CODE ---------------------------------#

# Inputs
o_f = 8                     # O/F Ratio
h = 2966.93 # J/g
H_t = 0.1 # J
det_v = 2921.3 # m/s
length = 0.159385 # m

# Constants
P_man = 4.5 * 10**5         # Manifold Pressure (Pa)
MW_o2 = 31.9988/1000        # Oxygen Molecular Weight (kg/mol)
MW_h2 = 2.01588/1000        # Hydrogen Molecular Weight (kg/mol)

# Calculations
t = length / (det_v/2) # s
m = H_t / h # g
mdot = m / t # g/s
mdo2 = mdot * (o_f / (o_f + 1)) 
mdh2 = mdot * (1 / (o_f + 1))

no2 = mdo2 / (MW_o2*1000)
nh2 = mdh2 / (MW_h2*1000)
xo2 = no2 / (no2 + nh2)
xh2 = nh2 / (no2 + nh2)
po2 = xo2 * P_man
ph2 = xh2 * P_man

pox_feed = 1.893 * po2
ph2_feed = 1.899 * ph2

cd = 0.61
gamma_o2 = 1.4
gamma_h2 = 1.41
R_o2 = 8.314/MW_o2
R_h2 = 8.314/MW_h2
temp = 293

#area_o2_WRONG = mdo2 / (1000 * cd * math.sqrt(gamma_o2 * (po2 / R_o2 * temp) * po2 * (2/(gamma_o2+1)**((gamma_o2+1)/(gamma_o2-1)))))
area_o2 = (mdo2/1000) * math.sqrt(temp) / pox_feed / math.sqrt(gamma_o2/R_o2) * ((gamma_o2+1)/2)**((gamma_o2+1)/(2*(gamma_o2-1)))
area_h2 = (mdh2/1000) * math.sqrt(temp) / ph2_feed / math.sqrt(gamma_h2/R_h2) * ((gamma_h2+1)/2)**((gamma_h2+1)/(2*(gamma_h2-1)))
rad_o2 = math.sqrt(area_o2 / math.pi)
rad_h2 = math.sqrt(area_h2 / math.pi)
dia_o2 = rad_o2 * 2 * 1000 # mm
dia_h2 = rad_h2 * 2 * 1000 # mm

#def f(m):)
#    return P_man - (((1/(o_f+1)) * m/MW_h2) * P_h2) + (((o_f/(o_f+1)) * m/MW_o2) * P_o2) / ((1/(o_f+1)*m)/MW_h2 + (o_f/(o_f+1)*m)/MW_o2)
#m_dot_total = fsolve(f, 0.1)
#m_dot_h2 = 1/(o_f+1) * m_dot_total
#m_dot_o2 = o_f/(o_f+1) * m_dot_total

#------------------------------ PRINT COMMANDS ------------------------------#
print(f"Mass Flow Rate: {mdot:.3f} g/s")
print(f"O2 P3: {po2 / psi_to_Pa:.3f} psia")
print(f"H2 P3: {ph2 / psi_to_Pa:.3f} psia")
print(f"O2 Mass Flow Rate: {mdo2:.3f} g/s")
print(f"H2 Mass Flow Rate: {mdh2:.3f} g/s")
print(f"Diameter of O2: {dia_o2/25.4:.6f} in")
print(f"Diameter of H2: {dia_h2/25.4:.6f} in")