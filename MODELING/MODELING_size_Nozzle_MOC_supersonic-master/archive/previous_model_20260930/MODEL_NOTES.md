# Straight-boundary RDC model

Run `addpath('RDE_Toolbox_2/src'); r=MoC_Main_Model;` from the project directory. Inputs remain at the top of `MoC_Main_Model.m`. Chamber/port dimensions are mm and mm²; solver quantities are SI and radians.

## Geometry and characteristic net

The slip line and oblique shock are straight again. Their slopes come from the triple-point shock/expansion compatibility calculation; pressure matching is enforced at the triple point only. `diagnostics.slipPressureMismatch` reports the subsequent pressure mismatch under this approximation. The curved pressure-matched solver is preserved in `archive/MoC_Coupled_Field.m`.

The detonation is tilted from the axial direction by `zeta=asin(Vinjection/Vcj)` (Grunenwald Eq. 2.21), with foot offset `waveHeight*tan(zeta)`. It is not forced vertical. Default injection velocity is low compared with CJ speed, so the calculated tilt is only about 1.6 degrees. Enlarging the angle requires a different physical injection/wave-speed ratio, not a cosmetic change to the wave geometry.

`numerics.displayPhase=0.58` moves the displayed detonation into the domain interior, similar to thesis Figure 15. This wraps one calculated pitch for display; it is not a periodic-flow closure. Field data and solver coordinates remain unchanged. Discontinuous wrap segments are split rather than connected across the plot.

The characteristic-only panel and `results/characteristic_net.png` draw continuous red C- and blue C+ paths through the interpolated wave-fixed direction field `theta +/- asin(1/Mwave)`. Traces stop at region boundaries. These are visualization paths, not new solver points or fan seeds. The actual characteristic nodes and edges remain in `r.mesh`. Traces in the bounding region use its prescribed uniform state. No characteristics are drawn in refill.

## Three different Mach quantities

1. `geometry.detonationMach = Vcj/a_reactants`: propagation Mach of the detonation relative to the reactants, about 5.5 for the default inputs.
2. `field.Mwave`: local **gas** Mach in the wave-fixed frame. Immediately behind a CJ detonation it is near one; the approved numerical seed is 1.02. This is not the propagation Mach of the wave.
3. `field.Mlab`: local laboratory-frame gas Mach. Refill is subsonic (about 0.153). The same refill gas is high-Mach in the wave-fixed frame because the wave moves rapidly relative to it.

Changing refill `Mwave` to 0.153, or changing the post-CJ gas seed to Mach 5 to represent the wave, would break the velocity triangle and gas dynamics. The plots now avoid those ambiguities: the Mach contour explicitly uses `Mlab`; the characteristic panel has no Mach colormap and labels the refill `Mlab`. Propagation Mach and CJ gas Mach are labeled separately in the figure title. The true `Mwave` remains available in the output.

## Solver and thermodynamics

The injector calculation conserves mass and adiabatic total enthalpy under the supplied momentum/area-change closure. It derives `Cp=gamma*R/(gamma-1)` and selects a physical subsonic root. Port pressure/temperature are static. Optional round-port diameter overrides port-set area. `N_holes` counts port sets.

CEA receives the computed injection pressure/temperature and returns CJ chemistry. It runs in a temporary directory without changing archived user working files. `chemistry.useCEA=false` selects explicit manual product properties. Subsequent regions use constant effective gamma and gas constant, rather than resolving chemistry or the ZND zone.

Product characteristics start at interpolated points along the tilted detonation, with `seedMach=1.02` avoiding the sonic degeneracy. The planar compatibility/intersection approach comes from the archived `src2` nozzle solver. The product mesh remeshes on successive lines and traces interpolated characteristic feet. The independent post-shock mesh starts from interpolated shock points. Under the straight-shock/uniform-upstream closure its constant state can be extended exactly outside its characteristic hull; such grid cells are marked `analyticClosure`.

Interpolation is separate for refill (region 1), products (2), shocked products (3), and bounding products (4). Product cells outside the computed hull remain NaN rather than receiving arbitrary extrapolation. Refill primitive pressure, temperature, density and velocity are assigned directly, and `nu`, `Kplus`, and `Kminus` remain NaN there. `nu` is the Prandtl-Meyer coordinate, not a scalar velocity potential.

## Limitations and verification

`refillFraction` prescribes the open injection fraction and nominal wave height. Injection `mdot` is the open-port flow, not a time-averaged blocked-injection mass flow. Bounding products are a prescribed uniform expanded state. The full periodic Mach/refill-height iteration from the thesis is not implemented. Restoring straight boundaries also restores the downstream contact-pressure mismatch; it is reported rather than hidden. Smooth contours and full grid coverage do not imply a validated pressure-matched flow field.

Run `addpath('RDE_Toolbox_2/tests'); report=run_tests;` for CEA parsing, mass/enthalpy conservation, equation of state, reference-frame conversion, refill NaNs, characteristic compatibility, straight/tilted geometry, invalid-input handling, and mesh-refinement checks. Previous curved-boundary validation numbers do not apply to this restored solver.

Outputs: `results/MoC_result.mat` (structure `saved`), `results/unwrapped_flow.png`, and `results/characteristic_net.png`. Node columns are `[x y theta Mwave P0wave]`; edge columns are `[parent child family]`. `geometry.curves` contains marching coordinate and the normal coordinates of the straight slip/shock boundaries. Global boundary coordinates are `slipXY` and `shockXY`.

Reference: John Andrew (Jack) Grunenwald (2023), *Investigation of Rotating Detonation Physics and Design of a Mixer for a Rotating Detonation Engine*, Chapter 2, Eqs. 2.7-2.28 and Figure 15. The supplied thesis and original nozzle-source license are retained.

The restored straight-boundary regressions passed in MATLAB R2026a. The 11-to-21-seed manual-chemistry pressure refinement check changed the common product-region mean-normalized pressure by about 3.85%. The default CEA run covered the entire grid; its propagation Mach was 5.499, refill laboratory Mach about 0.153, wave tilt 1.596 degrees, and downstream slip-pressure mismatch about 92.5%. These are numerical checks of the reduced model, not experimental validation.
