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
l1 = 1.5 * in_to_m # tube 2 length, m
l2 = 2.5 * in_to_m # tube 1 length, m
t1 = 0.065 # tube 1 thickness, in
t2 = 0.049 # tube 2 thickness, in
r1_o = 0.375/2 # inches
r2_o = 0.25/2 # inches

# Shchelkin Spiral Calculations
OD = 0.3 # inches
ID = 0.21 # inches
BR = (OD**2 - ID**2) / (OD**2)
print("Blockage Ratio:", str(BR))

# Hoop Stress Calculations
P_i = 3004 # psia
P_o = 14.7 # psia
r1_i = r1_o - t1 # inches
r2_i = r2_o - t2 # inches
sigma1_hoop = ((r1_i**2 * P_i - r1_o**2 * P_o) / (r1_o**2 - r1_i**2)) + (r1_i**2 * r1_o**2 * (P_i - P_o) / (r1_i**2 * (r1_o**2 - r1_i**2))) # psia
sigma2_hoop = ((r2_i**2 * P_i - r2_o**2 * P_o) / (r2_o**2 - r2_i**2)) + (r2_i**2 * r2_o**2 * (P_i - P_o) / (r2_i**2 * (r2_o**2 - r2_i**2))) # psia
print("Hoop Stress 3/8 in:", str(sigma1_hoop), "psia")
print("Hoop Stress 1/4 in:", str(sigma2_hoop), "psia")
hoop_stress_safety1 = 6500 # psia
hoop_stress_safety2 = 7500 # psia
print("Pressure Safety Factor 3/8 in:", str(hoop_stress_safety1/sigma1_hoop))
print("Pressure Safety Factor 1/4 in:", str(hoop_stress_safety2/sigma2_hoop))