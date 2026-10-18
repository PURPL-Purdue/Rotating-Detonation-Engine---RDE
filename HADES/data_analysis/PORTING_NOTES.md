# HADES data analysis: MATLAB to Python porting notes

Source: `HADES/DATA ANALYSIS/` (removed from the working tree after the port; `mkdir -p /tmp/rde_src && git archive 2a30c95 | tar -x -C /tmp/rde_src` restores it). Target: this package (`HADES.data_analysis`).
Tests: `tests/hades_data_analysis/`, run with `uv run pytest tests/hades_data_analysis`.

## File map

| Original | Python | Notes |
|---|---|---|
| `BOS_Background_Generation.m` | `bos_background_generation.py` (`BOS_Background_Generation`, `main`) | Writes `./BOS Backgrounds/` in the current directory, as MATLAB did with `pwd`. Override with `output_folder=` or `--output-dir` |
| `BOS Backgrounds/BOS_random_dot_background.{svg,png}` | `bos_print_assets/` (copied, for printing) and `tests/hades_data_analysis/reference/bos/` (read-only test oracle, pinned by SHA-256) | MATLAB output for the current inputs. The code never writes to either folder |
| `BOS_random_dot_background.svg` (top level) | `bos_print_assets/BOS_random_dot_background_FOV60mm.svg` (copied, renamed). The test pins it by SHA-256 only | An older run with `FOV_width_mm = FOV_height_mm = 60` (231483 dots). Byte-identical when regenerated with those inputs |
| `HeatFluxCalc.m` | `heat_flux_calc.py` (`HeatFluxCalc`, `main`) | Derivative bug fixed. No default model: `equations="cubic"` or `"legacy"`, plus `depths_m=`, are both required. See bug 5 and open question 1 |
| `sample_thermocouple_data.csv` | `sample_thermocouple_data.csv` (copied) | |
| `PURPL_test_fake_ITP_data_v2.m` | `purpl_fake_data.PURPL_test_fake_ITP_data_v2` | Canonical ITP generator (Sept 2026) |
| `PURPL_ITP_data_fft_new.m` | superseded | Identical to v2 apart from the function name. Uploaded 3 minutes before v2 |
| `PURPL_Fake_ITP_Pressure_Data.m` | `purpl_fake_data.PURPL_Fake_ITP_Pressure_Data` | April generator: 1 to 8 kHz tones, 3-sample delay, spikes |
| `PURPL_Fake_Plenum_pressure.m` | `purpl_fake_data.PURPL_Fake_Plenum_pressure` | Returns `[time; h2; air]` as 3 x N |
| `PURPL_test_FFT_v2.m` | `purpl_fft.PURPL_test_FFT_v2` | Canonical pipeline: averaged one-sided amplitude spectrum |
| `PURPL_FFT_new.m` | superseded | Identical to v2, but calls `PURPL_ITP_data_fft`, which does not exist, so it never ran |
| `ITP_FFT_for_PURPL_FakeData.m` | `purpl_fft.ITP_FFT_for_PURPL_FakeData` | Raw-magnitude pipeline. Frequency axis made general |
| `PURPLE_MassFlow_from_PlenumPressure.m` | `purpl_mass_flow_from_plenum_pressure.py` (`PURPLE_MassFlow_from_PlenumPressure`, `plenum_mass_flow`, `main`) | Row-swap bug fixed |
| `cD_analysis.m` | `choked_nozzle_analysis.cD_analysis` | |
| `temp_analysis.m` | `choked_nozzle_analysis.temp_analysis` | Shares constants with `cD_analysis` |
| `CDRConditions.py` | `cdr_conditions.py` | Same code. No import-time side effects. `plt.show()` moved to `main` |
| `UpstreamConditions.py` | `upstream_conditions.py` | Same numbers. CoolProp is cached per temperature (about 200x fewer calls) |
| `scfm.py` | `scfm.py` | Wrapped in functions. Output unchanged |
| `test.m` | dropped | Contains only `% whaddup`. It also shadows Octave's built-in `test` |
| (new) | `choked_flow.py` | Choked-flow equation shared by the scripts above |
| (new) | `rdelib.matlab_compat` | `matlab_round`, MATLAB-compatible random streams (`matlab_rng`, `as_rng`, `rand`, `randn`, `randi`); first written here as `_matlab_compat.py`, moved verbatim |
| (new) | `purpl_fft.single_sided_amplitude_spectrum`, `dominant_frequencies` | Reusable spectrum and peak picking for real ITP data |

Every module has a `main()` and runs with `uv run python -m HADES.data_analysis.<module>`.
Files are only ever written to the current directory or to an explicit output option:
* BOS: `./BOS Backgrounds/`, or `--output-dir`.
* v2 fake ITP data and the v2 FFT: `./synthetic_pressure_data.csv` (as MATLAB did),
  `--output-dir`, `--no-csv`, or `csv_path=`.
Functions take `plot=False` so they can run without figures. `rng=` sets the random
stream: `None` gives a fresh stream, an `int` seed gives MATLAB's `rng(seed)` uniform
stream, or pass a numpy Generator.

## Validation

* **Octave 11.3.0 with signal 1.4.8.** The original `.m` files ran unmodified. The
  harness in `tests/hades_data_analysis/_octave_harness.py` adds shims in a temp
  directory:
  * `rand`, `randn` and `randi` replay number streams generated in Python, so both
    sides draw identical random numbers.
  * `plot` records its data.
  * The other graphics calls are no-ops.
  * `readmatrix`, `writematrix`, `subtitle` and `sgtitle` are filled in.

  Results:
  * Fake ITP data v1 and v2, and fake plenum data: bit-identical.
  * v2 spectrum: within 3e-16. v1 raw spectrum: within 5e-16 relative. v1 frequency
    axis: within 7e-15 kHz.
  * Mass flow, using MATLAB's row assignment: bit-identical.
  * `cD_analysis` and `temp_analysis`: within 5e-16 relative.
  * `HeatFluxCalc`: scipy `savgol_filter(mode="interp")` matches Octave `sgolayfilt`
    within 4e-12 K. Eq. 3/4, fed MATLAB's derivative, match within 1e-14 relative.

  Golden values from these runs are hard-coded in the tests. The live Octave tests
  read the originals from `$RDE_MATLAB_SOURCE_ROOT/HADES/DATA ANALYSIS` and skip,
  with the recovery command, when Octave or the original `.m` files are missing:
  `RDE_MATLAB_SOURCE_ROOT=/tmp/rde_src uv run pytest tests/hades_data_analysis`.
* **MATLAB output.** The committed BOS SVG is byte-identical and the PNG
  (11811 x 8268) is pixel-identical. The oracles are read-only fixtures in
  `tests/hades_data_analysis/reference/bos/` plus SHA-256 hashes of the original
  files, and the code never writes there. numpy `RandomState(seed)` reproduces MATLAB's
  `rng(seed)` twister stream. MATLAB maps seed 0 to 5489. The stray 60 mm-FOV SVG is
  also byte-identical.
* **Original Python scripts.** The stdout of `CDRConditions.py`,
  `UpstreamConditions.py` and `scfm.py`, captured before clean-up, is stored in
  `tests/hades_data_analysis/reference/`. The cleaned-up output must match it exactly.
  The arrays also match the originals' namespaces within 1e-10 (CoolProp).
* **Hand and analytic checks**:
  * one-sided amplitude convention: a pure tone of amplitude A shows height A;
  * the FFT recovers injected tones within 0.5 Hz for several seeds, for both
    pipelines, including the delay and the spikes;
  * heat flux, "cubic" model: exact against heat-equation solutions, both directly
    and through the CSV / Savitzky-Golay pipeline. Example:
    `T = 300 - 20x + 50x^2 + 100*alpha*t` gives q0 = 7800 W/m^2;
  * heat flux, "legacy" model: steady-state and Eq. 3/4 coefficient checks;
  * one-sided spectrum: odd-N and even-N Nyquist cases;
  * the SG derivative is exact for quadratic data;
  * choked-flow hand calculations;
  * BOS dot count and Poisson coverage: 0.3028 measured, 0.3033 expected.

## MATLAB bugs fixed (before -> after)

1. **`HeatFluxCalc.m`: the "temperature derivative" was `T_s / dt`.** That is
   temperature divided by the time step, about 3e5 "K/s". It is now the Savitzky-Golay
   first derivative with the same order and window: scipy
   `savgol_filter(T, 21, 2, deriv=1, delta=dt)`. Results on `sample_thermocouple_data.csv`,
   using `equations="legacy"` and `depths_m=MATLAB_SENSOR_DEPTHS_M`:

   | t [s] | dT1/dt before | after [K/s] | T0 before | T0 after [K] | q0 before | q0 after [W/m^2] |
   |---|---|---|---|---|---|---|
   | 0.000 | 300453 | 39.09 | 306.0033 | 301.4365 | -3821.75 | 1528.37 |
   | 0.250 | 412280 | 296.90 | 428.6971 | 422.5817 | 8835.62 | 16056.43 |
   | 1.000 | 492883 | 12.49 | 500.4213 | 492.8992 | -8776.08 | 24.88 |
   | 1.999 | 499507 | -72.89 | 506.7397 | 499.1103 | -9545.56 | -620.84 |
   | mean q0 | | | | | -5083.46 | 3266.98 |

   The hard-coded data path `/Users/Abby/Desktop/RDE/...` now points at the copy of
   the CSV next to the module.
2. **`PURPLE_MassFlow_from_PlenumPressure.m`: plenum rows were swapped.**
   `PURPL_Fake_Plenum_pressure` returns `[time; h2; air]`, but the script read row 2
   as air and row 3 as H2. Both plenums sit at 20 bar, so only the oscillation
   pattern moved between the gases. Results on the replayed stream (seed 13), mean
   mass flow in lb/s:
   * air: 1.972392 before, 1.969316 after;
   * H2: 0.0491239 before, 0.0492007 after.

   The mass-flow figure was titled "Static Plenum Pressure Over Time". It is now
   "Mass Flow Rate Over Time".
3. **`ITP_FFT_for_PURPL_FakeData.m`: the frequency axis was
   `linspace(1, fs, fs+1)/1000`.** That is only correct when the record is exactly
   one second long (N == fs). It is now `k*fs/N`, k = 1..N/2. The default
   1 s / 100 kHz data are unchanged (difference 7e-15 kHz). On a 0.25 s record, a
   4 kHz tone would have been labelled about 1 kHz. Test:
   `test_v1_axis_fix_for_records_not_one_second_long`.
4. **`PURPL_FFT_new.m` calls the undefined function `PURPL_ITP_data_fft`.** It is
   superseded by `PURPL_test_FFT_v2`.
5. **`HeatFluxCalc.m` Eq. 3/4 are invalid, so there is no default model any more.**
   `HeatFluxCalc(data_file, *, equations, depths_m, ...)` requires both keywords. On
   the command line, `--equations {legacy,cubic}` and either `--depths-m X1 X2` or
   `--matlab-depths` are required.
   * `"cubic"`: the equations derived from the heat equation (open question 1).
   * `"legacy"`: the original Eq. 3/4 (`surface_*_legacy`). The docstring marks them
     known-invalid and they emit `LegacyEquationWarning`.
   * Depths are in metres. The MATLAB pair is kept as `MATLAB_SENSOR_DEPTHS_M =
     (0.25, 0.575)`, and using it emits `SensorDepthUnitsWarning` (units unconfirmed,
     probably mm).

   On the exact solution `T = 300 - 20x + 50x^2 + 100*alpha*t`:
   * legacy gives q0 = -8287.50 W/m^2;
   * cubic gives +7800.00 W/m^2, both directly and through the CSV pipeline;
   * the true value is +7800 W/m^2.
6. **Inherited bug in `PURPL_test_FFT_v2.m` / `PURPL_FFT_new.m`: odd-length
   records.** `p_single_sided(2:end-1) = 2*...` skips the last kept bin for every N.
   For odd N there is no Nyquist bin, so the highest positive-frequency amplitude
   was halved. Example: 9 samples of `3*cos(2*pi*4*n/9)` give 1.5 at bin 4 in Octave
   with the original lines. It now gives 3.0. The fix doubles `[1:]` for odd N and
   `[1:-1]` for even N. Even-length results, including all the 100000-sample fake
   data, are bit-identical to before (tests `test_odd_length_last_bin_is_doubled`,
   `test_even_length_nyquist_not_doubled_and_unchanged`).

## Open questions (the team should decide)

1. **`HeatFluxCalc.m` Eq. 3/4 look wrong.**
   * The transient terms `alpha*(...)*T_dot` have units K m^3/s^2, but the
     `(T2-T1)/(x2-x1)` term has K/m. They need `1/alpha`.
   * The transient terms also have the wrong sign.
   * The equations fail an exact solution of the heat equation. Uniform heating,
     `T = 300 + 2*alpha*c*t + b*x + c*x^2`, gives q0 = -8288 W/m^2 where the exact
     value is +7800 W/m^2 (T0 is 292.8 K instead of 300.0 K).
   * With x1 = 0, Eq. 3 gives T0 != T1.

   I derived replacements from `dT/dt = alpha*T_xx` with `dT/dt` linear in depth:
   `surface_temperature_cubic` and `surface_heat_flux_cubic`, selected with
   `equations="cubic"`. They are exact for the analytic solutions in the tests.
   Following the orchestrator's decision there is no default model (bug 5).

   Still open: the sensor depth units.
   * With `depths_m=(0.25, 0.575)` (MATLAB values read as metres), the cubic
     equations give unphysical values on the sample data: mean q0 = 1.4e8 W/m^2 and
     |T0| up to 4e5 K.
   * With `depths_m=(0.25e-3, 0.575e-3)`, they give mean q0 = 3.41 MW/m^2, and T0 is
     within 0.4 K of the legacy result.

   Please check the source paper's Eq. 3/4 and confirm the depths.
2. **`PURPL_Fake_Plenum_pressure.m`, kept as written:**
   * `sin(time .* f)` treats f as rad/s, not Hz. With `2*pi`, the 600 to 1200 Hz
     terms would alias at 1 kHz sampling.
   * The third H2 term reuses `f2`, so `f3` is drawn but never used.
   * The oscillation sum is rectified with `abs`.
3. **`PURPL_Fake_ITP_Pressure_Data.m`:** the comment says "scale to PSI (30-70
   range)", but the code gives 30 to 200 (sensor 1) and 30 to 218 (sensor 2), and
   the y-axis label says kPa. Kept as written.
4. **`PURPLE_MassFlow_from_PlenumPressure.m`** uses the static plenum pressure as
   the stagnation pressure. `stag_density_h2` and `stag_density_air` are defined
   but never used.
5. **Constants disagree between scripts** (each port keeps its script's own values):

   | Quantity | Values |
   |---|---|
   | H2/air stoichiometry | FAR 0.0292 (`cD`/`temp`); O/F 34.3 (`CDRConditions`, `scfm`); AFR 34.5 (`UpstreamConditions`) |
   | R_H2 [J/kg-K] | 4157 (`cD`/`temp`); 4124.2 (mass flow); 4126 (`scfm`); CoolProp (CDR, Upstream) |
   | Air injector area [mm^2] | 304.77 (mass flow); 304.52 = 0.472 in^2 (`cD`/`temp`/Upstream); 316.04 annulus (CDR) |
   | H2 injector area [mm^2] | 25.00 = 0.03875 in^2 (`cD`/`temp`); 25.13 (CDR, Upstream, mass flow) |
6. **`CDRConditions.py`:** the air Reynolds number puts the whole 0.907 kg/s
   through one 2 mm hydraulic diameter (Re = 3.2e7). The Fanno friction factor and
   pressures depend on this. It may need the per-passage flow and area.

## Shared helpers (`rdelib`)

* Done: `matlab_round`, `matlab_rng`, `as_rng` and the `rand`/`randn`/`randi`
  adapters moved verbatim to `rdelib.matlab_compat` (fake data and BOS outputs
  bit-identical; the BOS SVG/PNG are still pinned by SHA-256).
* Done: the friction factor in `cdr_conditions.fanno_pdrop` is
  `rdelib.gasdynamics.tkachenko_friction_factor` (the same expression, so the
  CDR results and printout are bit-identical).
* Not done: `choked_flow.py`. Its three forms are algebraically equal to
  `rdelib.gasdynamics.choked_flow_factor` / `choked_mass_flow` but use a different
  operation order, and the scripts are pinned to Octave at 5e-16 / bit-identical.
* Not done: the Mach solve in `cdr_conditions.fanno_pdrop` (it squares a Python
  float, `rdelib.gasdynamics.fanno_fLD` a numpy one, which is not guaranteed to
  round identically) and the inline unit factors (2.20462 lb/kg,
  0.453592 kg/lbm, 6.4516e-4 m^2/in^2, 2118.88 cfm per m^3/s, 0.000145038 psi/Pa,
  25.4 mm/in), which are kept as written in each script.

## Dependencies

PNG writing uses Pillow (`PIL`), declared in `pyproject.toml`. The SVG is written
with LF line endings on every OS, so it stays byte-identical on Windows.
