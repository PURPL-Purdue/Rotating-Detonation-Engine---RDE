"""Dataclasses replacing the MATLAB structs of the MoC nozzle solver.

Field names are identical to the MATLAB struct fields (``geom.yt``,
``params.gamma``, ``plots.patches_data``, ``chocked.mach`` ...).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class Geom:
    """Nozzle geometry and discretisation (MATLAB ``geom``)."""

    delta: int  # [0/1] 0: planar nozzle, 1: axisymmetric nozzle
    yt: float  # [m] throat radius (reference length)
    rhou: float  # [m] throat upstream radius of circular arc (IVLINE)
    rhod: float  # [m] throat downstream radius of circular arc
    xe: float  # [m] nozzle length
    xplume: float  # [m] axial length of the plume after the exit lip
    ta: float  # [deg] attachment angle between circular arc and wall contour
    te: float  # [deg] exit lip angle
    NI: int  # number of points on the initial-value line
    circdownTheta: np.ndarray  # [deg] discretisation of the downstream arc (e.g. 1:1:ta)
    NIexpansion: int  # number of Prandtl-Meyer expansion waves at the lip
    circdownX: np.ndarray | None = None  # [m] set by the driver
    circdownY: np.ndarray | None = None  # [m] set by the driver


@dataclass
class Params:
    """Flow parameters (MATLAB ``params``)."""

    gamma: float  # specific heat ratio Cp/Cv
    R: float  # [J/kg-K] gas constant
    P: float  # [Pa] stagnation pressure
    PRatio: float  # Pstatic_lip / Pstatic_ambient
    T: float  # [K] stagnation temperature
    Pstatic: float | None = None  # [Pa] ambient static pressure, set in extend_plume


@dataclass
class Plots:
    """Post-processing switches (MATLAB ``plots``)."""

    patches: int  # [0/1] plot the 2-D characteristic patches
    patches_data: int  # 0 outlines only, 1 Mach, 2 p, 3 T, 4 rho, 5 |V|, 6 flow angle
    patches_xlim: float  # abscissa above which patches are not plotted


@dataclass
class Chocked:
    """1-D isentropic choked-nozzle solution (MATLAB ``chocked``)."""

    mach: np.ndarray = field(default_factory=lambda: np.zeros(0))
    pressure: np.ndarray = field(default_factory=lambda: np.zeros(0))  # p / P_stag
