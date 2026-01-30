import math
import CoolProp.CoolProp as cp
from scipy.optimize import fsolve

#--------------------------------- CONVERSIONS -------------------------------#
psi_to_Pa = 6894.76
bar_to_psia = 14.5
choke_factor = 0.5283

#--------------------------------- MAIN CODE ---------------------------------#

# Inputs
P_o2 = 100 * psi_to_Pa      # Feed pressure of O2 (psia)
P_h2 = 100 * psi_to_Pa      # Feed pressure of H2 (psia)
o_f = 8                     # O/F Ratio

# Constants
P_man = 4.5 * 10^5          # Manifold Pressure (Pa)
MW_o2 = 15.999/1000         # Oxygen Molecular Weight (g/mol)
MW_h2 = 2.01588/1000        # Hydrogen Molecular Weight (g/mol)

# Calculations
def f(m):
    return P_man - (((1/(o_f+1)) * m/MW_h2) * P_h2) + (((o_f/(o_f+1)) * m/MW_o2) * P_o2) / ((1/(o_f+1)*m)/MW_h2 + (o_f/(o_f+1)*m)/MW_o2)

m_dot_total = fsolve(f, 0.1)
m_dot_h2 = 1/(o_f+1) * m_dot_total
m_dot_o2 = o_f/(o_f+1) * m_dot_total

#------------------------------ PRINT COMMANDS ------------------------------#
print(f"O2 Supply Pressure: {P_o2 / psi_to_Pa:.3f} psia")
print(f"O2 Mass Flow Rate: {m_dot_o2:.3f} kg/s")
print(f"H2 Supply Pressure: {P_h2 / psi_to_Pa:.3f} psia")
print(f"H2 Mass Flow Rate: {m_dot_h2:.3f} kg/s")