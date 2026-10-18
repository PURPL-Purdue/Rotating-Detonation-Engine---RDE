"""Port of ``HADES_size_P0_v1.m`` (initial chamber pressure estimate)."""

from __future__ import annotations

__all__ = ["HADES_size_P0_v1"]


def HADES_size_P0_v1():
    """Initial chamber pressure [bar] from the injector-area-ratio relation.

    All inputs are hard-coded as in the MATLAB: Ai/Aw = 0.1135, M0 = 0.65,
    plenum pressure 20 bar, gamma = 1.4.
    """
    Ai_Aw = 0.1135  # unitless
    M_0 = 0.65  # unitless
    P_p = 20  # bar
    gamma_unburned = 1.4
    P_c = (2 / (gamma_unburned + 1)) ** (gamma_unburned / (gamma_unburned - 1)) * P_p  # bar
    P_0 = (
        Ai_Aw
        * (1 / M_0)
        * ((gamma_unburned + 1) / (2 + (gamma_unburned - 1) * M_0**2)) ** (1 / 2)
        * P_c
    )  # bar
    print("Initial Chamber Pressure %.2f bar" % P_0)
    return P_0
