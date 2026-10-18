"""Dataclasses replacing the MATLAB structs of the RDE MoC model (same field names)."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class ChamberDim:
    """Annular chamber dimensions (MATLAB ``ChamberDim``), all in [mm]."""

    outer_chamber_diameter: float  # [mm]
    inner_chamber_diameter: float  # [mm]
    chamber_length: float  # [mm] axial length from injector face to exit


@dataclass
class InjV:
    """Injector inputs for ``MoC_Injection_Velocity`` (MATLAB ``InjV``)."""

    injector_area: float  # [mm^2] area of a single set of injector orifices
    N_holes: int  # number of injectors (1 for a plenum)
    P3: float  # [Pa] manifold pressure before the orifice
    T3: float  # [K] gas temperature
    mdot: float  # [kg/s]
    gamma1: float  # unreacted gas gamma
    Pa: float  # [Pa] chamber static pressure
    R: float  # [J/kg-K] gas constant
    Cp: float  # [J/kg-K] specific heat at constant pressure
    A4: float | None = None  # [m^2] annulus area (set locally by MoC_Injection_Velocity)


@dataclass
class ITripleParam:
    """Triple-point inputs (MATLAB ``iTripleParam``)."""

    P1: float  # [Pa] injection static pressure (= bounding gas initial pressure)
    P2: float  # [Pa] post-detonation static pressure
    T2: float  # [K] post-detonation static temperature
    Vcj: float  # [m/s] Chapman-Jouguet detonation velocity
    gamma2: float  # specific heat ratio of post-detonation gas
    R1: float  # [J/kg-K] specific gas constant of bounding gas


@dataclass
class InjectionVelocityOut:
    """Output struct ``out`` of ``MoC_Injection_Velocity``.

    ``V4``, ``P4``, ``T4``: both roots of the momentum/energy/continuity
    quadratic.  ``Vwave``: as built by the MATLAB (``out.Vwave(invalid) = NaN``
    on an undefined field: zeros up to the last invalid root, NaN at invalid
    roots, empty if none).  ``M4``: Mach number V4/sqrt(gamma1*R*T4) (NaN for
    T4 <= 0), the value printed in the 5th column (see PORTING_NOTES.md).
    """

    V4: np.ndarray
    P4: np.ndarray
    T4: np.ndarray
    Vwave: np.ndarray = field(default_factory=lambda: np.zeros(0))
    M4: np.ndarray = field(default_factory=lambda: np.zeros(0))


@dataclass
class PostProc:
    """Plot-domain extents (MATLAB ``PostProc``) in [mm]."""

    xMax: float
    yMax: float
