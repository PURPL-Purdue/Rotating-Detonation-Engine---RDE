"""Standard volumetric flow (SCFM) of the air and H2 supplies (cleaned-up ``scfm.py``).

Converts the design mass flows (2 lbm/s air; H2 at phi = 1.2) at 5.171 MPa
(750 psi) and ambient temperature into actual and standard cubic feet per
minute: ``SCFM = ACFM * (P_psia/14.7) * (519 R / T_R)``.

Run as a script::

    uv run python -m HADES.data_analysis.scfm
"""

from __future__ import annotations

# Global
max_phi = 1.2
stoic_OF = 34.3
T_amb = 293  # K
T_amb_F = 67.73  # degF

M3S_TO_CFM = 2118.8799727597  # (m^3/s) -> (ft^3/min)
PA_TO_PSI = 0.000145038

# Air
m_dot_air = 0.907185  # kg/s
p_air = 5.171e6  # Pa
R_air = 287  # J/kg-K

# H2
m_dot_h2 = m_dot_air / stoic_OF * max_phi  # kg/s
p_h2 = 5.171e6  # Pa
R_h2 = 4126  # J/kg-K


def scfm_from_mass_flow(m_dot, p, R, T=T_amb, T_F=T_amb_F):
    """Standard cubic feet per minute for ``m_dot`` [kg/s] stored at ``p`` [Pa], ``T`` [K]."""
    rho = p / (R * T)  # kg/m^3
    Q = m_dot / rho  # m^3/s
    Q_acfm = Q * M3S_TO_CFM  # acfm
    return Q_acfm * ((p * PA_TO_PSI) / 14.7) * (519 / (460 + (T_F)))


def scfm():
    """Returns ``(Q_air_scfm, Q_h2_scfm)``."""
    return scfm_from_mass_flow(m_dot_air, p_air, R_air), scfm_from_mass_flow(m_dot_h2, p_h2, R_h2)


def main() -> None:
    Q_air_scfm, Q_h2_scfm = scfm()
    print("Air:", str(Q_air_scfm), "scfm")
    print("H2:", str(Q_h2_scfm), "scfm")


if __name__ == "__main__":
    main()
