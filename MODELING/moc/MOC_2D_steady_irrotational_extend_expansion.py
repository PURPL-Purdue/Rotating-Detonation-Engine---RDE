"""Flow field determined by the wall downstream of the circular arc (direct wall points).

Python translation of ``MOC_2D_steady_irrotational_extend_expansion.m`` from
xavierdechamps/Nozzle_MOC_supersonic (GPL-3.0); see MODELING/moc/README.md.
All row/column indices below are MATLAB 1-based indices.
"""

from ._matlab import atand, mgt, mle, mlt
from ._net import Net
from .MOC_2D_steady_irrotational_get_geometry import MOC_2D_steady_irrotational_get_geometry
from .MOC_2D_steady_irrotational_internal_point import MOC_2D_steady_irrotational_internal_point
from .MOC_2D_steady_irrotational_wall import MOC_2D_steady_irrotational_wall
from .MOC_2D_steady_irrotational_wall_inverse import MOC_2D_steady_irrotational_wall_inverse

# Safety stop for the open-ended ``while 1`` marching loop (MATLAB would hang).
MAX_COLUMNS = 100000


def MOC_2D_steady_irrotational_extend_expansion(
    geom, params, indIin, LENG_INDIin, Xin, Yin, Uin, Vin,
    ind_triin, ind_quadin, i_triin, j_triin, i_quadin, j_quadin,
):
    """Direct wall-point method up to the exit lip, then an inverse wall point at ``xe``.

    Same input/output convention as the MATLAB function; returns
    ``(indI, X, Y, U, V, LENG_INDI, ind_tri, ind_quad, i_tri, j_tri, i_quad, j_quad)``.
    """
    indI = int(indIin)
    net = Net(Xin, Yin, Uin, Vin, LENG_INDIin, i_triin, j_triin, i_quadin, j_quadin)
    X, Y, U, V = net.X, net.Y, net.U, net.V
    LENG_INDI = net.LENG_INDI
    ind_tri = int(ind_triin)
    ind_quad = int(ind_quadin)
    internal = MOC_2D_steady_irrotational_internal_point
    intersection = None  # MATLAB variable persists across iterations of the while loop

    while True:
        indI = indI + 1
        if indI > MAX_COLUMNS:
            raise RuntimeError("extend_expansion: the wall never reached geom.xe")
        # Wall point
        J = 1
        LENG_INDI[indI] = 1
        # Point 2: left-running C+  ---> ( J+1 , indI-1 )
        net.set_point(J, indI, *MOC_2D_steady_irrotational_wall(*net.point(J + 1, indI - 1), geom, params))
        if mgt(X[J, indI], geom.xe):
            break
        ind_tri += 1
        net.i_tri.set(ind_tri, [indI - 1, indI - 1, indI])
        net.j_tri.set(ind_tri, [J, J + 1, J])

        J = 2
        for JJ in range(2, int(LENG_INDI[indI - 1]) - 1 + 1):  # Internal points
            intersection = False
            # Point 1: right-running C-  ---> ( J-1  , indI   )
            # Point 2: left -running C+  ---> ( JJ+1 , indI-1 )
            try:
                xtmp, ytmp, utmp, vtmp = internal(
                    *net.point(JJ + 1, indI - 1), *net.point(J - 1, indI), geom, params
                )
            except Exception:
                # Possible intersection of characteristics -> shock -> not supported by the MoC.
                # Stop the C- characteristic and copy the values from the upstream C- characteristic
                xtmp, ytmp, utmp, vtmp = net.point(JJ + 1, indI - 1)

            if mlt(xtmp, X[J - 1, indI]):
                print("Intersection of left-running characteristics inside the diverging region -> weak shock ")
                ind_tri += 1
                net.i_tri.set(ind_tri, [indI, indI - 1, indI - 1])
                net.j_tri.set(ind_tri, [J - 1, JJ, J + 1])
                continue  # Delete the current point and continue with the next one
            net.set_point(J, indI, xtmp, ytmp, utmp, vtmp)

            test1 = (X[J, indI] - X[J - 1, indI]) * (Y[J, indI - 1] - Y[J - 1, indI]) - (
                Y[J, indI] - Y[J - 1, indI]
            ) * (X[J, indI - 1] - X[J - 1, indI])
            test2 = (X[J, indI] - X[J - 1, indI]) * (Y[J + 1, indI - 1] - Y[J - 1, indI]) - (
                Y[J, indI] - Y[J - 1, indI]
            ) * (X[J + 1, indI - 1] - X[J - 1, indI])
            if mle(test1 * test2, 0.0):  # Intersection of characteristics
                intersection = True
                break

            LENG_INDI[indI] = LENG_INDI[indI] + 1
            ind_quad += 1
            net.i_quad.set(ind_quad, [indI, indI - 1, indI - 1, indI])
            net.j_quad.set(ind_quad, [J - 1, JJ, JJ + 1, J])
            J = J + 1

        if intersection is None:
            raise NameError("Undefined function or variable 'intersection'.")
        if intersection:
            # Copy the rest of the C- characteristic from the upstream C- characteristic
            loclen = int(LENG_INDI[indI - 1]) - J
            for k in range(loclen):
                net.set_point(J + k, indI, *net.point(J + 1 + k, indI - 1))
            LENG_INDI[indI] = LENG_INDI[indI] + loclen

            print("Intersection of right-running characteristics inside the diverging region -> weak shock")

            ind_tri += 1
            net.i_tri.set(ind_tri, [indI, indI - 1, indI])
            net.j_tri.set(ind_tri, [J - 1, J, J])
        else:
            # Point on the axis of symmetry
            net.set_point(J, indI, *internal(*net.mirror(J - 1, indI), *net.point(J - 1, indI), geom, params))
            Y[J, indI] = Y[J, indI] + 1.0e-6  # To avoid singularity on axis
            LENG_INDI[indI] = LENG_INDI[indI] + 1
            ind_tri += 1
            net.i_tri.set(ind_tri, [indI, indI - 1, indI])
            net.j_tri.set(ind_tri, [J - 1, int(LENG_INDI[indI - 1]), J])

    # The last point computed above lies downstream of the nozzle exit lip point.
    # Compute the C- from the exit lip point with the inverse wall method.
    _, y4, theta4 = MOC_2D_steady_irrotational_get_geometry(geom.xe, 2, geom)
    J = 1
    theta4 = atand(theta4)
    xtmp = [X[J + 1, indI - 1], 0.0, X[J, indI - 1], geom.xe]
    ytmp = [Y[J + 1, indI - 1], 0.0, Y[J, indI - 1], y4]
    utmp = [U[J + 1, indI - 1], 0.0, U[J, indI - 1], 0.0]
    vtmp = [V[J + 1, indI - 1], 0.0, V[J, indI - 1], 0.0]
    u4, v4, _, _, _, _ = MOC_2D_steady_irrotational_wall_inverse(xtmp, ytmp, utmp, vtmp, theta4, geom, params)
    U[J, indI] = u4
    V[J, indI] = v4
    X[J, indI] = geom.xe
    Y[J, indI] = y4

    J = None  # MATLAB: J = [] if the loop below is empty
    for J in range(2, int(LENG_INDI[indI - 1]) + 1):
        # Point 1: right-running C-  ---> ( J-1 , indI   )
        # Point 2: left -running C+  ---> ( J   , indI-1 )
        net.set_point(J, indI, *internal(*net.point(J, indI - 1), *net.point(J - 1, indI), geom, params))
        LENG_INDI[indI] = LENG_INDI[indI] + 1
        ind_quad += 1
        net.i_quad.set(ind_quad, [indI, indI - 1, indI - 1, indI])
        net.j_quad.set(ind_quad, [J - 1, J - 1, J, J])
    if J is None:
        raise IndexError("extend_expansion: empty upstream characteristic")
    # Point on axis of symmetry
    J = J + 1
    LENG_INDI[indI] = LENG_INDI[indI] + 1
    net.set_point(J, indI, *internal(*net.mirror(J - 1, indI), *net.point(J - 1, indI), geom, params))
    Y[J, indI] = Y[J, indI] + 1.0e-6  # To avoid singularity on axis
    ind_tri += 1
    net.i_tri.set(ind_tri, [indI, indI - 1, indI])
    net.j_tri.set(ind_tri, [J - 1, int(LENG_INDI[indI - 1]), J])

    return net.outputs(indI, ind_tri, ind_quad)
