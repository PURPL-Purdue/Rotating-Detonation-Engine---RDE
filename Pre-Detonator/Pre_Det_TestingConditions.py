import math
import CoolProp.CoolProp as cp

#--------------------------------- CONVERSIONS -------------------------------#
in_to_m = 0.0254
psi_to_Pa = 6894.76
bar_to_psia = 14.5
atm_to_psia = 14.7
kg_to_lbm = 2.20462

#--------------------------------- MAIN CODE ---------------------------------#

# Pre-Det Dimensions
quarter_length = (4.3955 + 0.74) * in_to_m # tube 1 length, m
threeeight_length = 1.5 * in_to_m # tube 2 length, m
wall_thickness = 0.049 # tube 1 thickness, in
wall_thickness2 = 0.065 # tube 2 thickness, in

# Shchelkin Spiral Calculations
OD = 0.3 # inches
ID = 0.21 # inches
BR = (OD**2 - ID**2) / (OD**2)
print("Blockage Ratio:", str(BR))

# Hoop Stress Calculations
P_i = 1360.14 # psia
P_o = 14.7 # psia
r_o = 0.25/2 # inches
r2_o = 0.375/2 # inches
r2_i = r2_o - wall_thickness2 # inches
r_i = r_o - wall_thickness # inches
hoop_stress = ((P_i * r_i**2 - P_o * r_o**2)/(r_o**2 - r_i**2)) + (r_i**2 * r_o**2 * (P_i - P_o) / (r_i**2 * (r_o**2 - r_i**2))) # psia
hoop_stress2 = ((P_i * r2_i**2 - P_o * r2_o**2)/(r2_o**2 - r2_i**2)) + (r2_i**2 * r2_o**2 * (P_i - P_o) / (r2_i**2 * (r2_o**2 - r2_i**2))) # psia
print("Hoop Stress 1/4 in:", str(hoop_stress), "psia")
print("Hoop Stress 3/8 in:", str(hoop_stress2), "psia")
hoop_stress_safety = 7514.7 # psia
hoop_stress_safety2 = 6514.7 # psia
print("Pressure Safety Factor 1/4 in:", str(hoop_stress_safety/hoop_stress))
print("Pressure Safety Factor 3/8 in:", str(hoop_stress_safety2/hoop_stress2))