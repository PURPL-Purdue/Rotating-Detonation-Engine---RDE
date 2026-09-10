"""
HADES_combined_heatmap.py
--------------------------
Produces two figures:
    Fig 1 — Air plenum pressure heatmap (single panel)
    Fig 2 — H2 plenum pressure heatmaps (one panel per phi = 0.8–1.2)

Both figures:
    X axis : discharge coefficient Cd  (0.5 – 0.9)
    Y axis : inlet stagnation temp T_0 (15 – 35 deg C)
    Color  : plenum pressure P_0 [bar]
    White band = 19.5 – 20.5 bar target window

gamma and R pulled from CoolProp at each T_0 for the respective fluid.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import CoolProp.CoolProp as CP

# Fixed inputs
M_DOT_AIR  = 2.000 * 0.453592     # kg/s  (2 lbm/s converted)
A_AIR      = 0.472 * 6.4516e-4    # m2    (0.472 in2 converted, CDR p.73)
A_H2       = 25.13e-6              # m2    (32 x 1mm orifices, CDR p.67)
P_REF      = 101325.0              # Pa    (1 atm, for CoolProp lookup)
AFR_STOICH = 34.5                  # stoichiometric air-to-fuel ratio for H2/air

# Sweep parameters 
PHI_VALS = [0.8, 0.9, 1.0, 1.1, 1.2]
CD_GRID  = np.linspace(0.50, 0.90, 200)
T0_GRID  = np.linspace(15.0, 35.0, 200)

# Highlight band and CDR design-point references 
P_LO, P_HI = 19.5, 20.5     # bar  target plenum pressure window
CDR_CD_AIR = 0.62            # CDR CFD-validated air Cd  (CDR p.75)
CDR_CD_H2  = 0.71            # CDR CFD-validated H2  Cd  (CDR p.69)
CDR_T0     = 20.0            # deg C  CDR inlet temperature


# Physics helpers 
def get_props(T_0_C, fluid):
    """Pull gamma and R = Cp - Cv from CoolProp at T_0_C, 1 atm."""
    T_K = T_0_C + 273.15
    Cp  = CP.PropsSI('Cpmass', 'T', T_K, 'P', P_REF, fluid)
    Cv  = CP.PropsSI('Cvmass', 'T', T_K, 'P', P_REF, fluid)
    return Cp / Cv, Cp - Cv

def calc_MFP(gamma):
    """Choked-flow mass flow parameter: m_dot = Cd * A * P0 * MFP / sqrt(R * T0)."""
    return np.sqrt(gamma) * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))

def calc_Pp_air(T_0_C, Cd):
    """Air plenum pressure [bar] from inverted choked flow equation."""
    T_K  = T_0_C + 273.15
    g, R = get_props(T_0_C, 'Air')
    return M_DOT_AIR * np.sqrt(R * T_K) / (Cd * A_AIR * calc_MFP(g)) / 1e5

def calc_Pp_H2(T_0_C, Cd, phi):
    """H2 plenum pressure [bar]; m_dot_H2 scales with phi via AFR_stoich."""
    T_K  = T_0_C + 273.15
    g, R = get_props(T_0_C, 'Hydrogen')
    mdot = M_DOT_AIR * phi / AFR_STOICH
    return mdot * np.sqrt(R * T_K) / (Cd * A_H2 * calc_MFP(g)) / 1e5


# Build 2D pressure grids 

Z_air = np.array([[calc_Pp_air(T, cd) for cd in CD_GRID] for T in T0_GRID])
print(f"  Air:  {Z_air.min():.1f}-{Z_air.max():.1f} bar")

Z_H2 = {}
for phi in PHI_VALS:
    Z_H2[phi] = np.array([[calc_Pp_H2(T, cd, phi) for cd in CD_GRID] for T in T0_GRID])
    print(f"  H2  phi={phi:.1f}:  {Z_H2[phi].min():.1f}-{Z_H2[phi].max():.1f} bar")


# Shared panel drawing function
def draw_panel(ax, Z, Cd_star, title):
    """Draw a heatmap panel with white 19.5-20.5 bar band and CDR marker."""
    z_min, z_max = Z.min(), Z.max()

    pcm = ax.pcolormesh(CD_GRID, T0_GRID, Z,
                        cmap='coolwarm', shading='auto',
                        vmin=z_min, vmax=z_max)

    # White highlight band — only draw if visible within this panel's range
    if z_min <= P_HI and z_max >= P_LO:
        ax.contourf(CD_GRID, T0_GRID, Z,
                    levels=[P_LO, P_HI], colors=['white'], alpha=0.25)
        cs = ax.contour(CD_GRID, T0_GRID, Z,
                        levels=[P_LO, P_HI],
                        colors=['white'], linewidths=1.8, linestyles='--')
        ax.clabel(cs, fmt='%.1f bar', fontsize=7.5, inline=True)
    else:
        ax.text(0.5, 0.5,
                f'19.5-20.5 bar\nnot in range\n(min = {z_min:.1f} bar)',
                transform=ax.transAxes, ha='center', va='center',
                fontsize=9, color='white', fontweight='medium',
                bbox=dict(boxstyle='round,pad=0.4', fc='#444', alpha=0.7))

    # Subtle 1-bar iso-pressure contours across full range
    ax.contour(CD_GRID, T0_GRID, Z,
               levels=np.arange(np.floor(z_min), np.ceil(z_max) + 1, 1),
               colors=['#cccccc'], linewidths=0.4, alpha=0.5)

    # CDR design-point star
    ax.scatter(Cd_star, CDR_T0, color='black', zorder=8, s=100, marker='*',
               label=f'CDR  ($C_d$ = {Cd_star})')

    ax.set_title(title, fontsize=11, fontweight='medium', pad=5)
    ax.set_xlabel('$C_d$', fontsize=10)
    ax.set_ylabel('$T_0$ (°C)', fontsize=10)
    ax.set_xlim(CD_GRID[0], CD_GRID[-1])
    ax.set_ylim(T0_GRID[0], T0_GRID[-1])
    ax.xaxis.set_major_locator(ticker.MultipleLocator(0.1))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.05))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))
    ax.tick_params(labelsize=9)
    ax.legend(fontsize=8, loc='upper right', framealpha=0.8)

    return pcm


# Figure 1 — Air 
fig1, ax1 = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
fig1.suptitle(
    'Air plenum pressure $P_0$ [bar]  --  white band = 19.5-20.5 bar\n'
    f'$\\dot{{m}}_{{air}}$ = {M_DOT_AIR:.4f} kg/s,  '
    f'$A_{{air}}$ = {A_AIR*1e6:.1f} mm$^2$',
    fontsize=12, fontweight='medium'
)
pcm1 = draw_panel(ax1, Z_air, CDR_CD_AIR, 'Air')
cbar1 = fig1.colorbar(pcm1, ax=ax1, pad=0.02, fraction=0.035)
cbar1.set_label('$P_0$ (bar)', fontsize=11)
cbar1.ax.tick_params(labelsize=9)


# Figure 2 — H2 (one subplot per phi) 
fig2, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True)
fig2.suptitle(
    'H$_2$ plenum pressure $P_0$ [bar]  --  white band = 19.5-20.5 bar\n'
    f'$\\dot{{m}}_{{air}}$ = {M_DOT_AIR:.4f} kg/s,  '
    f'$A_{{H_2}}$ = {A_H2*1e6:.2f} mm$^2$,  '
    f'AFR$_{{stoich}}$ = {AFR_STOICH}',
    fontsize=12, fontweight='medium'
)

for i, phi in enumerate(PHI_VALS):
    ax  = axes.flat[i]
    pcm = draw_panel(ax, Z_H2[phi], CDR_CD_H2, f'H$_2$  --  $\\phi$ = {phi:.1f}')
    cbar = fig2.colorbar(pcm, ax=ax, pad=0.02, fraction=0.046)
    cbar.set_label('$P_0$ (bar)', fontsize=9)
    cbar.ax.tick_params(labelsize=8)

# Hide unused 6th subplot (2x3 grid, 5 phi values)
axes.flat[5].set_visible(False)

plt.show()