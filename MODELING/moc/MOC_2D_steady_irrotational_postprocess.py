"""Post-processing of the MoC nozzle solution (characteristic net, wall/axis p and M).

Python translation of ``MOC_2D_steady_irrotational_postprocess.m`` (team version in
``src2/``, including the team's fix that colours the characteristic lines) from
xavierdechamps/Nozzle_MOC_supersonic (GPL-3.0); see MODELING/moc/README.md.
MATLAB ``patch`` calls are reproduced with matplotlib collections.
"""

from __future__ import annotations

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.tri import Triangulation

from ._matlab import atand
from .get_Laval_theory import get_Laval_theory
from .MOC_2D_steady_irrotational_get_thermo import MOC_2D_steady_irrotational_get_thermo


@dataclass
class PatchData:
    """Patch arrays assembled exactly as in the MATLAB post-processing (0-based faces)."""

    nodes_tri: np.ndarray
    faces_tri: np.ndarray
    nodesC_tri: np.ndarray
    nodesRGB_tri: np.ndarray
    facesC_tri: np.ndarray
    nodes_quad: np.ndarray
    faces_quad: np.ndarray
    nodesC_quad: np.ndarray
    nodesRGB_quad: np.ndarray
    facesC_quad: np.ndarray


def _output_patch(patches_data, Mach, pressure, temperature, density, U, V):
    if patches_data == 0:
        return False, Mach
    if patches_data == 1:
        return True, Mach
    if patches_data == 2:
        return True, pressure
    if patches_data == 3:
        return True, temperature
    if patches_data == 4:
        return True, density
    if patches_data == 5:
        return True, np.sqrt(U**2 + V**2)
    if patches_data == 6:
        # atand(V./U); complex-capable (plume nodes can be complex, see PORTING_NOTES)
        with np.errstate(divide="ignore", invalid="ignore"):
            return True, atand(V / U)
    raise ValueError("Unknown type of data to display")


def build_patches(X, Y, output_patch, ind_tri, ind_quad, i_tri, j_tri, i_quad, j_quad, LOC_PATCHES, patches_xlim):
    """Assemble the patch vertex/face/colour arrays of the MATLAB post-processing.

    ``i_*``/``j_*`` are the 1-based (column, row) indices returned by the solver.
    Faces are only kept when all their vertices satisfy ``x < patches_xlim``.
    The region colour index (black/red/blue/magenta for initial line / initial
    expansion / nozzle / plume) follows ``LOC_PATCHES`` like the MATLAB code.
    """
    colours = np.array([[0, 0, 0], [1, 0, 0], [0, 0, 1], [1, 0, 1]], dtype=float)
    LOC = np.asarray(LOC_PATCHES, dtype=int)

    def assemble(n_patch, nv, i_idx, j_idx, loc_col):
        nodes = np.zeros((nv * n_patch, 2))
        nodesC = np.full(nv * n_patch, np.nan)
        nodesRGB = np.zeros((nv * n_patch, 3))
        faces = []
        facesC = []
        indP2 = 0
        for I in range(1, n_patch + 1):
            rows = j_idx[:, I - 1] - 1
            cols = i_idx[:, I - 1] - 1
            base = nv * (I - 1)
            nodes[base : base + nv, 0] = np.real(X[rows, cols])
            nodes[base : base + nv, 1] = np.real(Y[rows, cols])
            if I == LOC[indP2, loc_col] + 1:
                indP2 += 1  # change of region [initial line; initial expansion; nozzle; plume]
            nodesRGB[base : base + nv, :] = colours[indP2]
            if np.sum(nodes[base : base + nv, 0] < patches_xlim) == nv:
                faces.append(np.arange(base, base + nv))
                facesC.append(colours[indP2])
                nodesC[base : base + nv] = np.real(output_patch[rows, cols])
        faces = np.asarray(faces, dtype=int).reshape(-1, nv)
        facesC = np.asarray(facesC, dtype=float).reshape(-1, 3)
        return nodes, faces, nodesC, nodesRGB, facesC

    nt, ft, ct, rgbt, fct = assemble(int(ind_tri), 3, np.asarray(i_tri), np.asarray(j_tri), 0)
    nq, fq, cq, rgbq, fcq = assemble(int(ind_quad), 4, np.asarray(i_quad), np.asarray(j_quad), 1)
    return PatchData(nt, ft, ct, rgbt, fct, nq, fq, cq, rgbq, fcq)


def axis_nodes(U):
    """Boolean mask of the last nonzero entry of each column of U (axis nodes)."""
    U = np.asarray(U)
    condplot = np.zeros(U.shape, dtype=bool)
    for I in range(U.shape[1]):
        nz = np.flatnonzero(U[:, I])
        if nz.size:
            condplot[nz[-1], I] = True
    return condplot


def MOC_2D_steady_irrotational_postprocess(
    geom, params, plots, indI, LENG_INDI, X, Y, U, V,
    ind_tri, ind_quad, i_tri, j_tri, i_quad, j_quad, LOC_PATCHES,
):
    """Plot the characteristic net, the wall/axis pressure and Mach number.

    Returns ``(fig1, fig2, fig3)`` matplotlib figures (MATLAB figures 1, 2, 3).
    Figures are not shown; call ``plt.show()`` to display them.
    """
    print("Postprocessing the results...")
    X = np.asarray(X)
    Y = np.asarray(Y)
    U = np.asarray(U)
    V = np.asarray(V)
    Mach, pressure, density, temperature, _ = MOC_2D_steady_irrotational_get_thermo(U, V, params)
    Xr, Yr = np.real(X), np.real(Y)

    # ---- Figure 1: characteristics ------------------------------------------------
    fig1 = plt.figure(1)
    fig1.clf()
    ax1 = fig1.add_subplot(111)
    ax1.plot(Xr[0, :indI] / geom.yt, Yr[0, :indI] / geom.yt, "k", linewidth=2)
    ax1.grid(True)
    ax1.set_xlabel("Axial location x / y_t")
    ax1.set_ylabel("Radial location y / y_t")

    if plots.patches > 0:
        plot_colours, output_patch = _output_patch(plots.patches_data, Mach, pressure, temperature, density, U, V)
        pd = build_patches(
            X, Y, output_patch, ind_tri, ind_quad, i_tri, j_tri, i_quad, j_quad, LOC_PATCHES, plots.patches_xlim
        )
        if plot_colours:
            # 'FaceColor','interp' -> Gouraud-shaded triangles (quads split in two)
            pts = np.vstack([pd.nodes_tri, pd.nodes_quad]) / geom.yt
            cvals = np.concatenate([pd.nodesC_tri, pd.nodesC_quad])
            off = pd.nodes_tri.shape[0]
            tris = [pd.faces_tri]
            if pd.faces_quad.size:
                fq = pd.faces_quad + off
                tris += [fq[:, [0, 1, 2]], fq[:, [0, 2, 3]]]
            tris = np.vstack(tris)
            # keep only the vertices of plotted faces (others carry NaN colour data)
            used, tris = np.unique(tris, return_inverse=True)
            tris = tris.reshape(-1, 3)
            cused = cvals[used]
            vmin, vmax = np.nanmin(cused), np.nanmax(cused)
            tpc = ax1.tripcolor(
                Triangulation(pts[used, 0], pts[used, 1], tris), cused, shading="gouraud",
                cmap="jet", vmin=vmin, vmax=vmax,
            )
            fig1.colorbar(tpc, ax=ax1)
        else:
            # Characteristic lines coloured by region ('EdgeColor','flat','FaceColor','none')
            for nodes, faces, facesC in (
                (pd.nodes_tri, pd.faces_tri, pd.facesC_tri),
                (pd.nodes_quad, pd.faces_quad, pd.facesC_quad),
            ):
                if faces.size:
                    verts = [nodes[f] / geom.yt for f in faces]
                    ax1.add_collection(
                        PolyCollection(verts, facecolors="none", edgecolors=facesC, linewidths=1)
                    )
        ax1.autoscale_view()
    ax1.set_aspect("equal")

    # ---- Figure 2: static pressure on the wall / axis vs 1-D theory -----------------
    fig2 = plt.figure(2)
    fig2.clf()
    ax2 = fig2.add_subplot(111)
    condplot = axis_nodes(U)
    wallnodesJ = np.arange(geom.NI + 1, X.shape[1] + 1) - 1  # 0-based columns
    wallnodesJ = wallnodesJ[Xr[0, wallnodesJ] < geom.xe]
    # MATLAB X(condplot) is column-major: one axis node per column, ordered by column
    ax_rows, ax_cols = np.nonzero(condplot.T)[::-1]
    chocked = get_Laval_theory(geom, params, Xr[0, wallnodesJ], Yr[0, wallnodesJ])
    p = np.real(pressure)
    ax2.semilogy(Xr[0, wallnodesJ] / geom.yt, p[0, wallnodesJ] / params.P, "r", linewidth=2)
    ax2.semilogy(Xr[ax_rows, ax_cols] / geom.yt, p[ax_rows, ax_cols] / params.P, "b", linewidth=2)
    ax2.semilogy(Xr[0, wallnodesJ] / geom.yt, chocked.pressure, "k", linewidth=2)
    ax2.grid(True, which="both")
    ax2.set_ylabel("Static pressure / Stagnation pressure [-]")
    ax2.set_xlabel("Axial location x / y_t")
    ax2.legend(["Wall", "Axis", "1D"])
    ax2.axis([0, geom.xe, 0.01, 1])
    ax2.set_xticks(np.arange(0, 11, 1))

    # ---- Figure 3: Mach number on the wall / axis vs 1-D theory ---------------------
    fig3 = plt.figure(3)
    fig3.clf()
    ax3 = fig3.add_subplot(111)
    M = np.real(Mach)
    ax3.plot(Xr[0, wallnodesJ] / geom.yt, M[0, wallnodesJ], "r", linewidth=2)
    ax3.plot(Xr[ax_rows, ax_cols] / geom.yt, M[ax_rows, ax_cols], "b", linewidth=2)
    ax3.plot(Xr[0, wallnodesJ] / geom.yt, chocked.mach, "k", linewidth=2)
    ax3.grid(True)
    ax3.set_ylabel("Mach number [-]")
    ax3.set_xlabel("Axial location x / y_t")
    ax3.axis([0, geom.xe, 1, 4])
    ax3.set_xticks(np.arange(0, 11, 1))
    ax3.legend(["Wall", "Axis", "1D"], loc="lower right")

    return fig1, fig2, fig3
