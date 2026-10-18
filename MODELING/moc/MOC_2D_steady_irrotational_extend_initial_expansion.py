"""Flow field determined by the initial-expansion contour (downstream circular arc).

Python translation of ``MOC_2D_steady_irrotational_extend_initial_expansion.m``
from xavierdechamps/Nozzle_MOC_supersonic (GPL-3.0); see MODELING/moc/README.md.
All row/column indices below are MATLAB 1-based indices.
"""

import numpy as np

from ._matlab import mge
from ._net import Net
from .MOC_2D_steady_irrotational_internal_point import MOC_2D_steady_irrotational_internal_point
from .MOC_2D_steady_irrotational_wall_inverse import MOC_2D_steady_irrotational_wall_inverse


def MOC_2D_steady_irrotational_extend_initial_expansion(
    geom, params, indIin, LENG_INDIin, Xin, Yin, Uin, Vin,
    ind_triin, ind_quadin, i_triin, j_triin, i_quadin, j_quadin,
):
    """Inverse wall-point method along the prescribed points of the downstream arc.

    Same input/output convention as the MATLAB function; returns
    ``(indI, X, Y, U, V, LENG_INDI, ind_tri, ind_quad, i_tri, j_tri, i_quad, j_quad)``.
    """
    indI = int(indIin)
    net = Net(Xin, Yin, Uin, Vin, LENG_INDIin, i_triin, j_triin, i_quadin, j_quadin)
    LENG_INDI = net.LENG_INDI
    ind_tri = int(ind_triin)
    ind_quad = int(ind_quadin)
    internal = MOC_2D_steady_irrotational_internal_point
    circdownX = np.ravel(geom.circdownX)
    circdownY = np.ravel(geom.circdownY)
    circdownTheta = np.ravel(geom.circdownTheta)

    for I in range(1, circdownX.size + 1):
        indI = indI + 1
        theta4 = circdownTheta[I - 1]
        LENG_INDI[indI] = 1
        J = 1
        # Prescribed point on the circular arc -> get the values of U and V
        xtmp = [net.X[J + 1, indI - 1], 0.0, net.X[J, indI - 1], circdownX[I - 1]]
        ytmp = [net.Y[J + 1, indI - 1], 0.0, net.Y[J, indI - 1], circdownY[I - 1]]
        utmp = [net.U[J + 1, indI - 1], 0.0, net.U[J, indI - 1], 0.0]
        vtmp = [net.V[J + 1, indI - 1], 0.0, net.V[J, indI - 1], 0.0]
        u4, v4, _, _, _, _ = MOC_2D_steady_irrotational_wall_inverse(
            xtmp, ytmp, utmp, vtmp, theta4, geom, params
        )
        net.U[J, indI] = u4
        net.V[J, indI] = v4
        net.X[J, indI] = circdownX[I - 1]
        net.Y[J, indI] = circdownY[I - 1]

        J = 2
        JJ = None  # MATLAB: empty loop leaves JJ = [] (and the j_tri assignment would fail)
        for JJ in range(2, int(LENG_INDI[indI - 1]) + 1):  # Internal points
            # Point 1: right-running C-  ---> ( J-1 , indI   )
            # Point 2: left -running C+  ---> ( JJ  , indI-1 )
            xt, yt_, ut, vt = internal(*net.point(JJ, indI - 1), *net.point(J - 1, indI), geom, params)
            if mge(yt_, net.Y[1, indI]):
                # The C+ characteristic exits the nozzle: delete the last computed point
                print("... Initial expansion: deleting a point outside of the nozzle")
                continue
            net.set_point(J, indI, xt, yt_, ut, vt)
            LENG_INDI[indI] = LENG_INDI[indI] + 1
            ind_quad += 1
            net.i_quad.set(ind_quad, [indI, indI - 1, indI - 1, indI])
            net.j_quad.set(ind_quad, [J - 1, JJ - 1, JJ, J])
            J = J + 1
        # Point on the axis of symmetry
        net.set_point(J, indI, *internal(*net.mirror(J - 1, indI), *net.point(J - 1, indI), geom, params))
        net.Y[J, indI] = net.Y[J, indI] + 1.0e-6  # To avoid singularity on axis
        LENG_INDI[indI] = LENG_INDI[indI] + 1
        ind_tri += 1
        if JJ is None:
            raise ValueError("Dimensions of arrays being concatenated are not consistent (JJ = []).")
        net.i_tri.set(ind_tri, [indI, indI - 1, indI])
        net.j_tri.set(ind_tri, [J - 1, JJ, J])

    return net.outputs(indI, ind_tri, ind_quad)
