"""Extend the solution from the initial-value line.

Python translation of ``MOC_2D_steady_irrotational_extend_initial_line.m`` from
xavierdechamps/Nozzle_MOC_supersonic (GPL-3.0); see MODELING/moc/README.md.
All row/column indices below are MATLAB 1-based indices.
"""

import numpy as np

from ._net import Net
from .MOC_2D_steady_irrotational_internal_point import MOC_2D_steady_irrotational_internal_point


def MOC_2D_steady_irrotational_extend_initial_line(
    geom, params, xvnull, ythroat, uvnull, vvnull, Xin, Yin, Uin, Vin
):
    """Fill X, Y, U, V from the initial-value (v = 0) line.

    Returns ``(indI, X, Y, U, V, LENG_INDI, ind_tri, ind_quad, i_tri, j_tri,
    i_quad, j_quad)``; ``i_*``/``j_*`` hold 1-based (column, row) node indices
    of the triangular/quadrilateral patches.
    """
    xvnull = np.ravel(xvnull)
    ythroat = np.ravel(ythroat)
    uvnull = np.ravel(uvnull)

    net = Net(Xin, Yin, Uin, Vin)
    LENG_INDI = net.LENG_INDI
    LENG_INDI[1] = 1
    ind_tri = 0
    ind_quad = 0
    internal = MOC_2D_steady_irrotational_internal_point

    I = 1
    for I in range(2, geom.NI + 1):
        net.set_point(1, I, xvnull[I - 1], ythroat[I - 1], uvnull[I - 1], vvnull)
        LENG_INDI[I] = 1
        J = 1
        for J in range(2, 2 * (I - 1) + 1):  # Internal points
            # Point 1: right-running C-  ---> ( J-1 , I   )
            # Point 2: left -running C+  ---> ( J-1 , I-1 )
            net.set_point(J, I, *internal(*net.point(J - 1, I - 1), *net.point(J - 1, I), geom, params))
            LENG_INDI[I] = LENG_INDI[I] + 1
            if J == 2:
                ind_tri += 1
                net.i_tri.set(ind_tri, [I, I - 1, I])
                net.j_tri.set(ind_tri, [J - 1, J - 1, J])
            else:
                ind_quad += 1
                net.i_quad.set(ind_quad, [I, I - 1, I - 1, I])
                net.j_quad.set(ind_quad, [J - 1, J - 2, J - 1, J])
        # Point on the axis of symmetry
        J = J + 1
        net.set_point(J, I, *internal(*net.mirror(J - 1, I), *net.point(J - 1, I), geom, params))
        net.Y[J, I] = net.Y[J, I] + 1.0e-6  # To avoid singularity on axis
        LENG_INDI[I] = LENG_INDI[I] + 1
        ind_tri += 1
        net.i_tri.set(ind_tri, [I, I - 1, I])
        net.j_tri.set(ind_tri, [J - 1, J - 2, J])

    indI = I
    return net.outputs(indI, ind_tri, ind_quad)
