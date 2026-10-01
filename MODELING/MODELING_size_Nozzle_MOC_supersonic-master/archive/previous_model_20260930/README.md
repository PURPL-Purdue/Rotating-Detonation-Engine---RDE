# Unwrapped rotating-detonation MoC model

The active entry point is `RDE_Toolbox_2/src/MoC_Main_Model.m`. From this directory:

```matlab
addpath('RDE_Toolbox_2/src');
r = MoC_Main_Model;
```

Inputs are grouped at the top of that function. Chamber dimensions and port area/diameter are in mm and mm²; the solver uses SI units and radians. The returned structure contains inputs, injector state, chemistry, boundary geometry, characteristic meshes, gridded primitive variables, and diagnostics. Plots and a MAT result are written under `RDE_Toolbox_2/results`.

The original nozzle solver and unrelated HADES sizing programs have been preserved under `archive/`, outside the active dependency chain. Existing CEA working files were archived rather than overwritten. Do not add the archive recursively to MATLAB's path.

Read `RDE_Toolbox_2/MODEL_NOTES.md` before interpreting a result as a physical prediction. Grid coverage is not a physical validation metric.

Reference: John Andrew (Jack) Grunenwald, *Investigation of Rotating Detonation Physics and Design of a Mixer for a Rotating Detonation Engine*, Purdue University, 2023, Chapter 2, especially Eqs. 2.7–2.28 and §2.3.1. The supplied thesis PDF is retained in this directory. [Purdue thesis record](https://hammer.purdue.edu/articles/thesis/INVESTIGATION_OF_ROTATING_DETONATION_PHYSICS_AND_DESIGN_OF_A_MIXER_FOR_A_ROTATING_DETONATION_ENGINE/24747450).

Original nozzle characteristic intersection/predictor-corrector logic: Xavier Dechamps, `Nozzle_MOC_supersonic`, retained in `archive/src2`. The existing repository license is retained.
## Running without plots or without CEA

```matlab
r = MoC_Main_Model(struct('numerics',struct('plot',false,'save',false)));
r = MoC_Main_Model(struct('chemistry',struct('useCEA',false)));
```

For regression checks:

```matlab
addpath('RDE_Toolbox_2/tests');
report = run_tests;
```

The program rejects invalid geometry, a sonic initial seed, nonphysical injection roots, detached/subsonic shock states, failed nonlinear compatibility, and insufficient marching extent.

The current default uses straight shock/slip boundaries. The characteristic-only figure is saved as `results/characteristic_net.png` and shows continuous C+/C- traces. `numerics.displayPhase` changes only the displayed phase of the wave. Gas Mach contours use the laboratory frame; detonation propagation Mach is reported separately. See the model notes for the downstream pressure mismatch inherent in the straight-boundary closure.
