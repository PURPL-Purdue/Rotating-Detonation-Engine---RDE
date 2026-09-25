import math
import numpy as np
import matplotlib.pyplot as plt
from CoolProp.CoolProp import PropsSI
from scipy.optimize import brentq

#--------------------------------- CONVERSIONS -------------------------------#
psi_to_Pa = 6894.76     # psi to pascals
in_to_m = 0.0254        # inches to meters

#--------------------------------- INPUTS ------------------------------------#
d_GOx = 0.01 * in_to_m          # GOx sonic orifice diameter (m)
d_GH2 = 0.01 * in_to_m          # GH2 sonic orifice diameter (m)
Cd_GOx = 0.61                   # GOx discharge coefficient
Cd_GH2 = 0.61                   # GH2 discharge coefficient
Tt = 10 + 273.15                # Supply (stagnation) temperature (K)
o_f_target = 8                  # Target O/F ratio

p_GOx_min = 200                 # psia
p_GOx_max = 1000                # psia
n_points = 81

#--------------------------------- CONSTANTS ---------------------------------#
Rbar = 8.3144                   # J/mol-K
MW_GOx = 31.9988 / 1000         # kg/mol
MW_GH2 = 2.01588 / 1000         # kg/mol
R_GOx = Rbar / MW_GOx
R_GH2 = Rbar / MW_GH2

A_GOx = math.pi * (d_GOx / 2)**2
A_GH2 = math.pi * (d_GH2 / 2)**2

#--------------------------------- FUNCTIONS ---------------------------------#
def choked_mdot(p1, fluid, R, A, Cd):
    """Choked mass flow (kg/s) through a sonic orifice. p1 in Pa."""
    cp = PropsSI('CPMASS', 'T', Tt, 'P', p1, fluid)
    cv = PropsSI('CVMASS', 'T', Tt, 'P', p1, fluid)
    g = cp / cv
    return Cd * (A * p1 / math.sqrt(Tt)) * math.sqrt(g / R) \
        * ((g + 1) / 2)**(-(g + 1) / (2 * (g - 1)))

def required_GH2_pressure(p_GOx_psia):
    """Returns GH2 supply pressure (psia) that gives O/F = o_f_target."""
    md_GOx = choked_mdot(p_GOx_psia * psi_to_Pa, 'Oxygen', R_GOx, A_GOx, Cd_GOx)
    md_GH2_target = md_GOx / o_f_target

    def residual(p_GH2_psia):
        md_GH2 = choked_mdot(p_GH2_psia * psi_to_Pa, 'Hydrogen', R_GH2, A_GH2, Cd_GH2)
        return md_GH2 - md_GH2_target

    return brentq(residual, 1, 10000)

#--------------------------------- SWEEP -------------------------------------#
p_GOx_range = np.linspace(p_GOx_min, p_GOx_max, n_points)
p_GH2_req = np.array([required_GH2_pressure(p) for p in p_GOx_range])

plt.figure(figsize=(8, 5))
plt.plot(p_GOx_range, p_GH2_req, 'b-', linewidth=2)
plt.xlabel('GOx Supply Pressure (psia)')
plt.ylabel('Required GH2 Supply Pressure (psia)')
plt.title(f'GH2 Supply Pressure for O/F = {o_f_target}')
plt.grid(True, alpha=0.4)
plt.xlim(p_GOx_min, p_GOx_max)
plt.tight_layout()
plt.savefig('GH2_vs_GOx_pressure.png', dpi=200)
plt.show()

#--------------------------------- USER PROMPT -------------------------------#
while True:
    entry = input("\nEnter GOx supply pressure in psia (or 'q' to quit): ").strip()
    if entry.lower() == 'q':
        break
    try:
        p_in = float(entry)
        if p_in <= 0:
            print("Pressure has to be positive.")
            continue
        p_out = required_GH2_pressure(p_in)
        md_O = choked_mdot(p_in * psi_to_Pa, 'Oxygen', R_GOx, A_GOx, Cd_GOx)
        md_H = choked_mdot(p_out * psi_to_Pa, 'Hydrogen', R_GH2, A_GH2, Cd_GH2)
        x_GOx = md_O / MW_GOx
        x_GH2 = md_H / MW_GH2
        n_GOx = x_GOx / (x_GH2 + x_GOx)
        n_GH2 = x_GH2 / (x_GH2 + x_GOx)
        cph2 = int(PropsSI('CPMASS', 'T', Tt, 'P', p_out, "Hydrogen"))
        cvh2 = int(PropsSI('CVMASS', 'T', Tt, 'P', p_out, "Hydrogen"))
        cpo2 = int(PropsSI('CPMASS', 'T', Tt, 'P', p_in, "Oxygen"))
        cvo2 = int(PropsSI('CVMASS', 'T', Tt, 'P', p_in, "Oxygen"))
        gh2 = cph2 / cvh2
        go2 = cpo2 / cvo2
        p2_GH2 = p_out / ((1 + (gh2 - 1) / 2) ** (gh2 / (gh2 - 1)))
        p2_GOx = p_in / ((1 + (go2 - 1) / 2) ** (go2 / (go2 - 1)))
        p_total = (p2_GH2*n_GH2) + (p2_GOx*n_GOx)
        print(f"Total pressure: {p_total:.2f} psia")
        print(f"Required GH2 pressure: {p_out:.2f} psia")
        print(f"Check -> O2: {md_O*1000:.5f} g/s, H2: {md_H*1000:.5f} g/s, O/F = {md_O/md_H:.4f}")
        False
    except ValueError:
        print("Invalid input, enter a number.")
