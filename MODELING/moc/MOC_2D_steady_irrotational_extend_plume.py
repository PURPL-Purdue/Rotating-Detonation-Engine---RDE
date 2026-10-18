"""Exit-lip Prandtl-Meyer expansion (or weak shock) and free-pressure plume.

Python translation of ``MOC_2D_steady_irrotational_extend_plume.m`` from
xavierdechamps/Nozzle_MOC_supersonic (GPL-3.0); see MODELING/moc/README.md.
All row/column indices below are MATLAB 1-based indices.
"""

import dataclasses

from ._matlab import cosd, fzero, num2str, mgt, mle, mlt, sind
from ._net import Net
from .FinalWaveExpansion import FinalWaveExpansion
from .get_Mach_from_nu import get_Mach_from_nu
from .get_prandtl_meyer_function import get_prandtl_meyer_function
from .MOC_2D_steady_irrotational_free_pressure import MOC_2D_steady_irrotational_free_pressure
from .MOC_2D_steady_irrotational_get_thermo import MOC_2D_steady_irrotational_get_thermo
from .MOC_2D_steady_irrotational_get_thermo2 import MOC_2D_steady_irrotational_get_thermo2
from .MOC_2D_steady_irrotational_internal_point import MOC_2D_steady_irrotational_internal_point

# Safety stop for the open-ended ``while 1`` marching loop (MATLAB would hang).
MAX_COLUMNS = 100000


def MOC_2D_steady_irrotational_extend_plume(
    geom, params, indIin, LENG_INDIin, Xin, Yin, Uin, Vin,
    ind_triin, ind_quadin, i_triin, j_triin, i_quadin, j_quadin,
):
    """Plume computation downstream of the exit lip.

    If ``params.PRatio`` (p_lip / p_ambient) >= 1 a Prandtl-Meyer fan of
    ``geom.NIexpansion`` waves is centred on the lip; 0.5 < PRatio < 1 is
    treated as a weak oblique shock (no special treatment); PRatio <= 0.5
    stops the computation (outputs returned unchanged).  The jet boundary is
    then marched with the free-pressure point method up to ``xe + xplume``.
    ``params`` is not modified (MATLAB value semantics); the ambient pressure
    is set on a local copy.
    """
    indI = int(indIin)
    net = Net(Xin, Yin, Uin, Vin, LENG_INDIin, i_triin, j_triin, i_quadin, j_quadin)
    X, Y, U, V = net.X, net.Y, net.U, net.V
    LENG_INDI = net.LENG_INDI
    ind_tri = int(ind_triin)
    ind_quad = int(ind_quadin)
    internal = MOC_2D_steady_irrotational_internal_point

    # Thermodynamic state at the exit lip point (row 1, column indI)
    Mach_lip, pressure_lip, _, _, _ = MOC_2D_steady_irrotational_get_thermo(U[1, indI], V[1, indI], params)
    Mach_lip = Mach_lip[()]
    pressure_lip = pressure_lip[()]
    plume_case = 1  # [1/2] 1: weak shock at exit lip point, 2: Prandtl-Meyer expansion
    if params.PRatio < 1:
        print("Oblique shock wave emanating from the exit lip point")
        if params.PRatio > 0.5:
            print("Weak shock wave -> MOC can continue in the plume without special treatment")
        else:
            print("Strong shock wave at the lip point. The MOC cannot handle such problems.")
            print("Stopping the computation now...")
            return net.outputs(indI, ind_tri, ind_quad)
    else:
        print("Prandtl-Meyer expansion wave at the exit lip point")
        print("Need to iterate to find the final angle that matches the static pressure outside the nozzle.")
        plume_case = 2

    params = dataclasses.replace(params, Pstatic=pressure_lip / params.PRatio)

    if plume_case == 2:  # Prandtl-Meyer expansion
        # Iterate on the final angle so that the static pressure matches the ambient one
        nu1, _ = get_prandtl_meyer_function(Mach_lip, params.gamma)  # before the expansion
        deflection = 30  # deflection angle [deg], nu2 = nu1 + deflection
        print("    Mach number before the expansion: " + num2str(Mach_lip))
        deflection = fzero(
            lambda delta12: FinalWaveExpansion(delta12, params.gamma, nu1, params.P, params.Pstatic),
            (0.0, float(deflection)),
            tolx=1e-12,
        )
        Mach2 = get_Mach_from_nu(nu1 + deflection, params.gamma)
        print("    Mach number after the expansion: " + num2str(Mach2))
        print("    Deflection angle (degrees)       " + num2str(deflection))

        incrdefl = deflection / geom.NIexpansion
        for I in range(1, geom.NIexpansion - 1 + 1):  # 1st to penultimate expansion wave
            # Data for the exit lip point for the current expansion wave
            Mach2 = get_Mach_from_nu(nu1 + incrdefl * I, params.gamma)
            _, _, _, a2 = MOC_2D_steady_irrotational_get_thermo2(Mach2, params)
            J = 1
            indI = indI + 1
            LENG_INDI[indI] = 1
            X[J, indI] = X[J, indI - 1]
            Y[J, indI] = Y[J, indI - 1]
            U[J, indI] = Mach2 * a2 * cosd(geom.te + incrdefl * I)
            V[J, indI] = Mach2 * a2 * sind(geom.te + incrdefl * I)

            # Internal points
            J = None
            for J in range(2, int(LENG_INDI[indI - 1]) + 1):
                # Point 1: right-running C-  ---> ( J-1 , indI   )
                # Point 2: left -running C+  ---> ( J   , indI-1 )
                net.set_point(J, indI, *internal(*net.point(J, indI - 1), *net.point(J - 1, indI), geom, params))
                LENG_INDI[indI] = LENG_INDI[indI] + 1
                ind_quad += 1
                net.i_quad.set(ind_quad, [indI, indI - 1, indI - 1, indI])
                net.j_quad.set(ind_quad, [J - 1, J - 1, J, J])
            if J is None:
                raise IndexError("extend_plume: empty upstream characteristic")
            # Point on the axis of symmetry
            J = J + 1
            net.set_point(J, indI, *internal(*net.mirror(J - 1, indI), *net.point(J - 1, indI), geom, params))
            Y[J, indI] = Y[J, indI] + 1.0e-6  # To avoid singularity on axis
            LENG_INDI[indI] = LENG_INDI[indI] + 1
            ind_tri += 1
            net.i_tri.set(ind_tri, [indI, indI - 1, indI])
            net.j_tri.set(ind_tri, [J - 1, int(LENG_INDI[indI - 1]), J])

    # Final stage of the Prandtl-Meyer expansion, or just after the weak shock
    intersection = None  # MATLAB variable persists across iterations of the while loop
    while True:
        indI = indI + 1
        if indI > MAX_COLUMNS:
            raise RuntimeError("extend_plume: the jet boundary never reached xe + xplume")
        # Jet point
        J = 1
        LENG_INDI[indI] = 1
        # Point 2: left-running C+  ---> ( J+1 , indI-1 )
        net.set_point(
            J, indI,
            *MOC_2D_steady_irrotational_free_pressure(
                *net.point(J, indI - 1), *net.point(J + 1, indI - 1), geom, params
            ),
        )
        if mgt(X[J, indI], geom.xe + geom.xplume):
            indI = indI - 1
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
                print("Intersection of left-running characteristics inside the plume -> weak shock ")
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

            print("Intersection of right-running characteristics inside the plume -> weak shock")

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

    return net.outputs(indI, ind_tri, ind_quad)
