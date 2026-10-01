# RDE model framework

This project is a starting framework for a hand-written MATLAB model. **Only the NASA CEA call runs.** Injection, velocity triangles, wave geometry, characteristics, interpolation, and plots are intentionally unimplemented.

Open `RDE_Toolbox_2/src/MoC_Main_Model.m` and press Run, or from this project directory:

```matlab
run('RDE_Toolbox_2/src/MoC_Main_Model.m');
```

Main is a readable script again. It restores the original named inputs (`ChamberDim`, `InjV`, `iTripleParam`, fuel/oxidizer, pressure, and temperature), calls CEA, and leaves `ceaOut` in the workspace. Paths are resolved from main's location, so it also works when MATLAB's current folder is elsewhere.

The existing source filenames remain as short TODO placeholders. Their calls in main are commented out. Calling a placeholder raises `MoC:NotImplemented`; it does not return fabricated states or plots. Implement and verify one stage at a time, then enable its call in main.

The CEA executable and libraries remain in `RDE_Toolbox_2/src/CEA`. Each main run saves readable CEA input/output under `RDE_Toolbox_2/results/CEA`. See [MODEL_NOTES.md](RDE_Toolbox_2/MODEL_NOTES.md) for how this interface works.

The previous implemented model, tests, notes, and generated results are preserved in `archive/previous_model_20260930`. Older nozzle/HADES references remain elsewhere in `archive`. Do not add the archive recursively to MATLAB's path.

To check the CEA interface only:

```matlab
addpath('RDE_Toolbox_2/tests');
run_tests;
```
