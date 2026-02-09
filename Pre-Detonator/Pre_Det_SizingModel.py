import math
import CoolProp.CoolProp as cp

#--------------------------------- CONVERSIONS -------------------------------#
psi_to_Pa = 6894.76
bar_to_psia = 14.5038
in_to_m = 0.0254

#--------------------------------- MAIN CODE ---------------------------------#

# Feed Conditions
P_h2 = 200 * psi_to_Pa                      # H2 pressure (Pa)
P_o2 = 300 * psi_to_Pa                      # O2 pressure (Pa)
T_h2 = 10 + 273.15                          # H2 temperature (K)
T_o2 = 10 + 273.15                          # O2 temperature (K)

gamma_h2 = cp.PropsSI("CPMOLAR", "T", T_h2, "P", P_h2, "H2") / cp.PropsSI("CVMOLAR", "T", T_h2, "P", P_h2, "H2")
gamma_o2 = cp.PropsSI("CPMOLAR", "T", T_o2, "P", P_o2, "O2") / cp.PropsSI("CVMOLAR", "T", T_o2, "P", P_o2, "O2")

R = 8314                                    # Universal gas constant (J/kg-K)
MW_h2 = 2.016                               # MW of H2 (kg/kmol)
R_h2 = R / MW_h2                            # Gas constant for H2
MW_o2 = 32                                  # MW of O2 (kg/kmol)
R_o2 = R / MW_o2                            # Gas constant for O2

# Orifice Sizes
D_h2 = 0.040 * in_to_m                      # H2 orifice diameter (m)
D_o2 = 0.040 * in_to_m                      # O2 orifice diameter (m)
A_h2 = math.pi * (D_h2/2)**2                # H2 area (m^2)
A_o2 = math.pi * (D_o2/2)**2                # O2 area (m^2)

# Mass flow rates
m_dot_h2 = (A_h2 * P_h2 / T_h2) * math.sqrt(gamma_h2/R_h2) * (((gamma_h2 + 1)/2) ** ((-gamma_h2-1)/(2 * (gamma_h2-1))))     # kg/s
m_dot_o2 = (A_o2 * P_o2 / T_o2) * math.sqrt(gamma_o2/R_o2) * (((gamma_o2 + 1)/2) ** ((-gamma_o2-1)/(2 * (gamma_o2-1))))     # kg/s
o_f = m_dot_o2 / m_dot_h2

# Print statements
print()
print("H2 Mass Flow Rate:", str(m_dot_h2*1000), "g/s")
print("O2 Mass Flow Rate:", str(m_dot_o2*1000), "g/s")
print("O/F Ratio:", str(o_f))
print()