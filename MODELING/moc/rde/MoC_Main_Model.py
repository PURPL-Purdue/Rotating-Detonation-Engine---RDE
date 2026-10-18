"""RDE method-of-characteristics model: main driver.

Python translation of ``RDE_Toolbox_2/MoC_Main_Model.m``.  Run with
``uv run python -m MODELING.moc.rde.MoC_Main_Model``.

Description (from the MATLAB): in addition to running the entire program this
main file initializes by computing the 2D injector velocity and CJ speed. It
constructs the first velocity triangle to determine the angle and height of
the detonation wave, before interpolating starting points along the
detonation wave and drawing characteristic lines from them.
Definitions: IVLine = initial value line = detonation wave front.

Logic process (still to be implemented in the MATLAB as well): draw first
velocity triangle -> draw main lines (detonation line, slip line, injector
fill line, shock line) -> define conditions in each zone -> draw
characteristic lines -> draw more characteristic lines in the slip line and
oblique shock.
"""

from __future__ import annotations

from .MoC_Calculate_Domain_Size import MoC_Calculate_Domain_Size
from .MoC_Injection_Velocity import MoC_Injection_Velocity
from .MoC_Resolve_Triple_Point import MoC_Resolve_Triple_Point
from .post_processing.Post_Process_Main import Post_Process_Main
from .structs import ChamberDim, InjV, ITripleParam, PostProc


def default_inputs():
    """Input section of ``MoC_Main_Model.m``; returns ``(ChamberDim, InjV, iTripleParam)``."""
    # Cell sizer from CEA should go here to render the unwrapping hand inputs obsolete.
    # Annulus mean diameter / circumference for X and Y max length (unwrapping RDE).
    chamber = ChamberDim(
        outer_chamber_diameter=52.243,  # [mm]
        inner_chamber_diameter=37.167,  # [mm]
        chamber_length=29.21,  # [mm] injector face to exit
    )
    # Triple point inputs
    triple = ITripleParam(
        P1=250000,  # injection static pressure (= bounding gas initial pressure) [Pa]
        P2=8e6,  # post-detonation static pressure [Pa]
        T2=3905,  # post-detonation static temperature [K]
        Vcj=2563,  # Chapman-Jouguet detonation wave velocity [m/s]
        gamma2=1.1405,  # specific heat ratio of post-detonation gas
        R1=430.47,  # specific gas constant of bounding gas [J/kg-K]
    )
    # Delta and beta angles are relative to phi angle. delta + beta = phi angle.
    inj = InjV(
        injector_area=3.7138,  # [mm^2] one set of injector orifices (doublet: sum of both)
        N_holes=36,  # number of injectors (1 if using a plenum)
        P3=1e6,  # [Pa] initial manifold pressure before entering orifice
        T3=398,  # [K] initial gas temperature
        mdot=0.59,  # [kg/s]
        gamma1=1.3569,  # unreacted gas gamma
        Pa=1e6,  # [Pa] static pressure of chamber
        R=315,  # gas constant
        Cp=315,  # specific heat at constant pressure (as in the MATLAB input)
    )
    return chamber, inj, triple


def run(chamber: ChamberDim, inj: InjV, triple: ITripleParam, plot: bool = True):
    """Run the model; returns ``(out, (beta, delta, P_match, M3, M1_p), PostProc)``."""
    print("\n--- RUNNING INJECTION VELOCITY CALCULATIONS ---")
    out = MoC_Injection_Velocity(inj, chamber)

    print("\n--- RUNNING TRIPLE POINT RESOLUTION ---")
    beta, delta, P_match, M3, M1_p = MoC_Resolve_Triple_Point(triple)

    # Single-line formatted printout of the triple-point outputs
    print(
        "Beta Shock: %.2f° | Slip Angle: %.2f° | P Matched: %.2f Pa | M3: %.2f | M1': %.2f\n"
        % (beta, delta, P_match, M3, M1_p)
    )

    # Post-processing: graph domain
    print("Graphing Results...", end="")
    xMax, yMax = MoC_Calculate_Domain_Size(chamber)
    post = PostProc(xMax=xMax, yMax=yMax)
    if plot:
        Post_Process_Main(post)
    return out, (beta, delta, P_match, M3, M1_p), post


def main():
    """Equivalent of running ``MoC_Main_Model.m``."""
    import matplotlib.pyplot as plt

    plt.close("all")
    results = run(*default_inputs())
    plt.show()
    return results


if __name__ == "__main__":
    main()
