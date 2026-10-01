# Framework and CEA interface

The active program is now a blank modeling framework. It restores the original main-file inputs and section sequence, with paths adapted to the current structure. The previous flow solver is archived; its approximations and verification claims do not apply to this framework.

## What currently runs

`src/MoC_Main_Model.m` is a script, not the previous options/result function. Running it clears the workspace, defines the original inputs, locates CEA, runs one chemistry calculation, and leaves the inputs and `ceaOut` in the workspace. All proposed model calls below the CEA section are commented out.

CEA uses the original independent conditions: CH4/O2, phi 1.3, 2.5 bar, and 283 K. It does not use a calculated injection state, and it does not overwrite the hand-entered `iTripleParam` values. Optional assignments from CEA to the model inputs appear as comments for you to review and enable.

`InjV.Cp = 315` is restored from the original main for review. It is unused; together with `InjV.R = 315` it is not a consistent ideal-gas property pair. Resolve the intended gas properties and the static/stagnation meaning of the injection inputs before using them. No replacement injection equations have been supplied.

## What FCEA2 is

`FCEA2.exe` is the bundled NASA Chemical Equilibrium with Applications executable. It performs the requested thermochemical calculation. It is not a MATLAB function, and it does not calculate the two-dimensional RDE characteristic field.

The MATLAB interface consists of two files:

- `getCEAPath.m` builds the full path to `src/CEA/FCEA2.exe` using its own file location. It does not depend on MATLAB's current directory.
- `HADES_size_ceaDet.m` writes the input, starts the executable, and parses its text output into a MATLAB structure.

The wrapper remains the existing CEA interface, with comments and optional retained run files added for inspection. The executable itself is unchanged.

## One call, step by step

1. **Read the arguments.** Main passes fuel, oxidizer, exactly one mixture definition, initial reactant pressure and temperature, the executable path, and the folder for retained output.
2. **Convert units.** The wrapper converts pressure to psia and temperature to K. The legacy argument names `P0` and `T0` refer here to the initial unburned reactant state; they do not automatically denote stagnation conditions.
3. **Prepare a working folder.** It creates a fresh temporary folder and copies `FCEA2.exe`, `thermo.lib`, and `trans.lib` into it. These are the executable, thermodynamic database, and transport database. A fresh directory avoids accidentally reading output from an earlier run.
4. **Write `cea_det.inp`.** The `problem` / `det` lines request a detonation calculation. The file specifies initial pressure/temperature, mixture ratio, and reactants. Output options request a short report with transport properties.
5. **Run FCEA2.** The executable asks for the input base name. `run.txt` contains `cea_det`, without the `.inp` extension. The wrapper uses a shell command equivalent to:

   ```text
   "full-path-to-FCEA2.exe" < run.txt > cea_console.txt
   ```

   `<` sends the base name to the executable's standard input. `>` saves its console messages. CEA reads `cea_det.inp` and produces `cea_det.out`.
6. **Read the answer.** The wrapper checks the exit status and output file, then searches labeled lines in `cea_det.out` for the requested values. MATLAB is extracting numbers, not repeating the chemistry calculation.
7. **Keep the audit files.** Main supplies `outputDir`, so the wrapper copies the input, output, and console log to `results/CEA`. Those filenames are overwritten by the next run. It returns to the original working directory and removes only its temporary run folder.

The returned `inputFile` and `outputFile` fields point to the retained files. If the wrapper is called without `outputDir`, these two fields are empty and temporary files are removed after parsing.

## Inspecting the result

After running main:

```matlab
ceaOut.cjVel
ceaOut.detMach
ceaOut.P_burned_bar
ceaOut.T_cj
edit(ceaOut.inputFile)
edit(ceaOut.outputFile)
```

| MATLAB field | Meaning | Units |
| --- | --- | --- |
| `cjVel` | CJ detonation propagation velocity | m/s |
| `detMach` | CEA detonation Mach number | dimensionless |
| `P_burned_bar` | Burned-gas pressure | bar |
| `T_cj` | Burned-gas temperature | K |
| `R_specific` | Burned-gas specific gas constant, from molecular weight | J/(kg K) |
| `gamma_unburned`, `gamma_burned` | CEA-reported gas gammas | dimensionless |
| `Son_speed_unburned`, `Son_speed_burned` | CEA-reported sound speeds | m/s |
| `P_ratio`, `T_ratio`, `rho_ratio` | Burned/unburned ratios | dimensionless |
| `Cp_eq`, `Cp_frozen` | Equilibrium/frozen heat capacity | kJ/(kg K) |

The parser retains the units printed by CEA for transport quantities; inspect their headings in the raw report before using them. It converts molecular weight to `R_specific` using the universal gas constant. Missing parsed values remain NaN; they are not substitute physical states.

The original CH4/O2 example returns roughly 2565 m/s CJ velocity and 3904 K burned temperature. Those are results for that chemistry and initial state, not injected gas velocity or an RDE flow-field solution.

## Where to hand-code the model

Keep the main sections as a review sequence: dimensions, injection/velocity triangle, detonation/triple point, IVLine seeds, product characteristics, shock-region characteristics, refill/interpolation, and plots. Each existing `MoC_*.m` helper is now a short placeholder stating its purpose. Its current signature is a starting suggestion, not a fixed data schema. Write your derivation next to each implementation and define units and reference frames explicitly.

No wave-height closure, startup Mach offset, Prandtl-Meyer formula, shock-polar solver, mesh marching, interpolation, or display-phase adjustment runs in the framework. A placeholder raises `MoC:NotImplemented` if called accidentally. The old numerical tests have been replaced with a small CEA-interface smoke check; no model validation is claimed.

The original nozzle program remains in `archive/src2`. The Grunenwald thesis remains at the project root. The most recent implemented model is preserved in `archive/previous_model_20260930`, including its old documentation and results. Its CEA binaries are not duplicated; if you ever run that archived copy independently, provide the active bundled CEA path explicitly.
