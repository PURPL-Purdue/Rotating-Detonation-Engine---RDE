import math

#--------------------------------- CONVERSIONS -------------------------------#
in_to_m = 0.0254
psi_to_Pa = 6894.76
bar_to_psia = 14.5
atm_to_psia = 14.7
kg_to_lbm = 2.20462

#--------------------------------- MAIN CODE ---------------------------------#
# Pre-Det Dimensions
l = 6.275 * in_to_m # tube length, ms
t = 0.049 # tube thickness, in
r_o = 0.25/2 # inches
r_i = r_o - t # inches

print("Tube Length:", 6.275, "in")
print("Tube Outer Diameter:", 2*r_o, "in")
print("Tube Inner Diameter:", 2*r_i, "in")

# Hoop Stress Calculations
P_i = 2300 # set to keep Pressure Safety Factor > 1.5
P_o = 14.7 # psia
sigma_hoop = ((r_i**2 * P_i - r_o**2 * P_o) / (r_o**2 - r_i**2)) + (r_i**2 * r_o**2 * (P_i - P_o) / (r_i**2 * (r_o**2 - r_i**2))) # psi
print("Hoop Stress 1/4 in:", str(sigma_hoop), "psi")
hoop_stress_safety = 7500 # psi
print("Pressure Safety Factor 1/4 in:", str(hoop_stress_safety/sigma_hoop))
