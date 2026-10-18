"""Driver: supersonic isentropic flow in a nozzle of known shape by the MoC.

Python translation of the team's ``MOC_2D_steady_irrotational_main.m``
(MODELING_size_Nozzle_MOC_supersonic-master) from xavierdechamps/Nozzle_MOC_supersonic
(GPL-3.0); see MODELING/moc/README.md.

Run as a script (``uv run python -m MODELING.moc.MOC_2D_steady_irrotational_main``)
to reproduce the MATLAB script: it solves the team's default case and shows
the three figures.  ``run()`` is the reusable, non-plotting entry point.
"""

from __future__ import annotations

import dataclasses
import warnings
from dataclasses import dataclass

import numpy as np

from ._matlab import colon, cosd, sind
from .MOC_2D_steady_irrotational_extend_expansion import MOC_2D_steady_irrotational_extend_expansion
from .MOC_2D_steady_irrotational_extend_initial_expansion import (
    MOC_2D_steady_irrotational_extend_initial_expansion,
)
from .MOC_2D_steady_irrotational_extend_initial_line import MOC_2D_steady_irrotational_extend_initial_line
from .MOC_2D_steady_irrotational_extend_plume import MOC_2D_steady_irrotational_extend_plume
from .MOC_2D_steady_irrotational_get_thermo import MOC_2D_steady_irrotational_get_thermo
from .MOC_2D_steady_irrotational_IVLINE import MOC_2D_steady_irrotational_IVLINE
from .structs import Geom, Params, Plots


@dataclass
class MocSolution:
    """Workspace of the MATLAB script after the solve (same variable names).

    ``X, Y, U, V`` are (row J, column I) node arrays, zero padded like MATLAB;
    column I is one C- characteristic, row 1 the wall / jet boundary, row
    ``LENG_INDI[I-1]`` the axis.  ``patch_*`` hold 1-based MATLAB indices.
    """

    geom: Geom
    params: Params
    ythroat: np.ndarray
    xsonic: np.ndarray
    xvnull: np.ndarray
    uvnull: np.ndarray
    vvnull: float
    indI: int
    X: np.ndarray
    Y: np.ndarray
    U: np.ndarray
    V: np.ndarray
    LENG_INDI: np.ndarray
    ind_patch_tri: int
    ind_patch_quad: int
    patch_i_tri: np.ndarray
    patch_j_tri: np.ndarray
    patch_i_quad: np.ndarray
    patch_j_quad: np.ndarray
    LOC_PATCHES: np.ndarray
    fig1: object = None
    fig2: object = None
    fig3: object = None

    def thermo(self):
        """``(Mach, pressure, density, temperature, sound)`` at every node."""
        return MOC_2D_steady_irrotational_get_thermo(self.U, self.V, self.params)


def default_inputs():
    """Input section of the team's ``MOC_2D_steady_irrotational_main.m``.

    Returns ``(geom, params, plots)``.
    """
    # Parameters for the flow field
    params = Params(
        gamma=1.2,  # Specific heat ratio = Cp / Cv
        R=320,  # [J/kg-K] Gas constant
        P=70.0e5,  # [Pa] Stagnation pressure
        PRatio=2,  # Pstatic_lip / Pstatic_ambiant [>1 to have a Prandtl-Meyer expansion]
        T=3000,  # [K] Stagnation temperature
    )

    # Parameters for the geometry
    yt = 1.0  # [m] Throat radius -> reference length for the whole nozzle
    ta = 15  # [deg] Attachment angle between circular arc and line
    geom = Geom(
        delta=0,  # [0/1] 0: planar nozzle, 1: axisymmetric nozzle
        yt=yt,
        rhou=1 * yt,  # [m] Throat upstream radius of circular arc (IVLINE)
        rhod=1e-6 * yt,  # [m] Throat downstream radius of circular arc
        xe=10 * yt,  # [m] Nozzle length
        xplume=10 * yt,  # [m] Axial length of the plume from the exit lip point
        ta=ta,
        te=15,  # [deg] Exit lip point angle
        # Parameters for the discretization
        NI=11,  # Number of points on the initial-value line
        circdownTheta=colon(1, 1, ta),  # Discretization of the downstream circular arc
        NIexpansion=5,  # Number of expansion waves, if any exists
    )

    # Parameters for the post-processing of the results
    plots = Plots(
        patches=1,  # [0/1] Plot the 2D patches (characteristics + solution)
        patches_data=0,  # 0 outlines, 1 Mach, 2 p, 3 T, 4 rho, 5 |V|, 6 flow angle [deg]
        patches_xlim=40 * yt,  # Abscissa above which the patch plot is cut
    )
    return geom, params, plots


def run(geom: Geom, params: Params, plots: Plots | None = None) -> MocSolution:
    """Solve the nozzle + plume (body of the MATLAB script after the inputs).

    ``plots=None`` skips post-processing; otherwise the three figures are
    created (not shown) and stored on the returned :class:`MocSolution`.
    """
    # Geometry of the downstream circular arc
    theta = np.asarray(geom.circdownTheta, dtype=float)
    geom = dataclasses.replace(
        geom,
        circdownX=geom.rhod * sind(theta),
        circdownY=geom.yt + geom.rhod * (1 - cosd(theta)),
    )

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        # The initial-value line is chosen as the line where V = 0
        ythroat = colon(0, geom.yt / (geom.NI - 1), geom.yt)
        # xsonic: abscissae of the sonic line (not used); xvnull: abscissae of V = 0;
        # uvnull: value of U on the V = 0 line
        xsonic, xvnull, uvnull = MOC_2D_steady_irrotational_IVLINE(geom, params, ythroat)
        vvnull = 0.0

        ythroat = ythroat.copy()
        uvnull = uvnull.copy()
        ythroat[0] = ythroat[0] + 1.0e-6  # To avoid singularity on axis
        uvnull[0] = uvnull[0] + 1.0e-6  # To avoid singularity on axis

        X = np.array([[xvnull[0]]])
        Y = np.array([[ythroat[0]]])
        U = np.array([[uvnull[0]]])
        V = np.array([[vvnull]])

        LOC_PATCHES = np.zeros((4, 2), dtype=int)

        print("Computing the extent from the initial-value line...")
        out = MOC_2D_steady_irrotational_extend_initial_line(
            geom, params, xvnull, ythroat, uvnull, vvnull, X, Y, U, V
        )
        LOC_PATCHES[0, :] = [out[6], out[7]]

        print("Computing the extent of the flow field determined by the initial-expansion contour...")
        out = MOC_2D_steady_irrotational_extend_initial_expansion(geom, params, out[0], out[5], *out[1:5], *out[6:])
        LOC_PATCHES[1, :] = [out[6], out[7]]

        print("Computing the extent of the flow field determined by the wall...")
        out = MOC_2D_steady_irrotational_extend_expansion(geom, params, out[0], out[5], *out[1:5], *out[6:])
        LOC_PATCHES[2, :] = [out[6], out[7]]

        out = MOC_2D_steady_irrotational_extend_plume(geom, params, out[0], out[5], *out[1:5], *out[6:])
        LOC_PATCHES[3, :] = [out[6], out[7]]

    indI, X, Y, U, V, LENG_INDI, ind_tri, ind_quad, i_tri, j_tri, i_quad, j_quad = out
    if np.iscomplexobj(X) or np.iscomplexobj(U):
        bad = np.flatnonzero(np.any((np.imag(X) != 0) | (np.imag(Y) != 0) | (np.imag(U) != 0) | (np.imag(V) != 0), axis=0))
        warnings.warn(
            f"MoC produced complex-valued nodes from column {bad[0] + 1} on (no real free-pressure "
            "solution in the plume); the original MATLAB propagates them silently. Those nodes are "
            "not physical.",
            RuntimeWarning,
            stacklevel=2,
        )
    sol = MocSolution(
        geom=geom, params=params, ythroat=ythroat, xsonic=xsonic, xvnull=xvnull, uvnull=uvnull,
        vvnull=vvnull, indI=indI, X=X, Y=Y, U=U, V=V, LENG_INDI=LENG_INDI,
        ind_patch_tri=ind_tri, ind_patch_quad=ind_quad, patch_i_tri=i_tri, patch_j_tri=j_tri,
        patch_i_quad=i_quad, patch_j_quad=j_quad, LOC_PATCHES=LOC_PATCHES,
    )

    if plots is not None:
        from .MOC_2D_steady_irrotational_postprocess import MOC_2D_steady_irrotational_postprocess

        print("Plotting the results...")
        sol.fig1, sol.fig2, sol.fig3 = MOC_2D_steady_irrotational_postprocess(
            geom, params, plots, indI, LENG_INDI, X, Y, U, V,
            ind_tri, ind_quad, i_tri, j_tri, i_quad, j_quad, LOC_PATCHES,
        )
    return sol


def main():
    """Equivalent of running ``MOC_2D_steady_irrotational_main.m``."""
    import matplotlib.pyplot as plt

    plt.close("all")
    geom, params, plots = default_inputs()
    sol = run(geom, params, plots)
    # sol.fig1.savefig("Characteristics.pdf"); sol.fig2.savefig("Pressure.pdf")
    # sol.fig3.savefig("Mach_number.pdf")
    plt.show()
    return sol


if __name__ == "__main__":
    main()
