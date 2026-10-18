# HADES sizing model: MATLAB to Python porting notes

Source: `HADES/SIZING MODEL/` (MATLAB; removed from the working tree after the port, `mkdir -p /tmp/rde_src && git archive 2a30c95 | tar -x -C /tmp/rde_src` restores it). Target: `HADES/sizing_model/` + shared helpers in `rdelib/`.
Tests: `tests/hades_sizing/` (115 tests). Run everything with
`uv run python -m HADES.sizing_model.HADES_size_main`.

## 1. File mapping

| Original | Python | Notes |
|---|---|---|
| `HADES_size_main_v1.m` | `HADES_size_main.py` (`HADES_size_main()`, `main()`) | Inputs are keyword args with the MATLAB values as defaults. Returns `MainResults`. |
| `HADES_size_ceaDet.m` | `HADES_size_ceaDet.py` (+ `rdelib/cea.py`) | Name/value pairs become kwargs (MATLAB-style positional pairs also accepted). Names are matched like `inputParser`: case-insensitive, exact match first, then a unique prefix (`'fu'` -> `fuel`); ambiguous names raise. Returns `CeaDetResults` with the struct's field names. |
| `HADES_size_ceaRocket.m` | `HADES_size_ceaRocket.py` (+ `rdelib/cea.py`) | Legacy fields keep the MATLAB (throat) meaning; both stations are also exposed under unambiguous names, see 4a. Returns `CeaRocketResults`. |
| `HADES_size_convert.m` | `HADES_size_convert.py` -> `rdelib.units.convert` | **Bug fix** for `'atm'`, see 4b. |
| `HADES_size_geometry.m` | `HADES_size_geometry.py` | |
| `HADES_size_chamberPressure.m` | `HADES_size_chamberPressure.py` | `0:0.001:T` via `rdelib.matlab_compat.colon`. |
| `HADES_size_annulusgap_fillheight.m` | `HADES_size_annulusgap_fillheight.py` | `readmatrix(xlsx)` -> openpyxl. The data file is copied to `HADES/sizing_model/HADES_size_cellsizedata.xlsx` (identical md5). `fit(...,'a./x')` -> closed-form least squares (plot only). |
| `HADES_size_HoopStressTemps.m` | `HADES_size_HoopStressTemps.py` | `fzero(f, x0)` via `rdelib.matlab_compat.fzero` (same bracket search as MATLAB). Returns `(Fail_temp_C, Fail_temp_F)`. |
| `HADES_size_P0_v1.m` | `HADES_size_P0_v1.py` | |
| `HADES_size_h2Inj.m` | `HADES_size_h2Inj.py` | `"txt" + number` printing reproduced with MATLAB `num2str` formatting. |
| `HADES_size_airInj.m` | `HADES_size_airInj.py` | Returns `(air_area, Area_Annulus)`. |
| `HADES_size_max_flux_range.m` | `HADES_size_max_flux_range.py` | |
| `HADES_size_P_plenum_fanno.m` | `HADES_size_P_plenum_fanno.py` | Fanno `fsolve` -> `scipy.optimize.fsolve` (same start point); Colebrook `fsolve` from 0.02 -> bracketed `rdelib.gasdynamics.colebrook` (unique root, see 3); `fzero` on bracket -> `brentq` on `rdelib.gasdynamics.subcritical_orifice_mass_flow`. |
| `HADES_size_pdrop_fanno.m` | `HADES_size_pdrop_fanno.py` | Returns `(P_inlet, P_exit, Ma_inlet, f)`. |
| `HADES_size_throat.m` | `HADES_size_throat.py` | |
| `HADES_size_thrust.m` | `HADES_size_thrust.py` | Default = MATLAB formula. Opt-in `station="exit"` and `exit_thrust_from_cea(...)` include the pressure term, see 4a. |
| `HADES_size_chambCooltime.m` (v1) | `HADES_size_chambCooltime.py` | Kept: constant-h model and the heat-loss breakdown are unique to v1. `ode45` -> `rdelib.ode.ode45`. |
| `HADES_size_chambCooltime_v3.m` | `HADES_size_chambCooltime_v3.py` | Four scenarios. |
| `HADES_size_chambCooltime_v2.m` | superseded by v3 | Same physics, renamed variables. Octave gives identical numbers for v2 and v3 (asserted in the tests). |
| `HADES_size_thermal_v1.m` | `HADES_size_thermal_v1.py` | `pdepe` -> `rdelib.ode.pdepe`. `phiIndex` stays 1-based and must be an integer in 1..4 (MATLAB would error on 0 or negative values; Python indexing would silently wrap). |
| `HADES_size_sonicnozzlearea_air.m` | `HADES_size_sonicnozzlearea_air.py` | **Bug fix** `R = 8.134 -> 8.314`, see 4c. |
| `HADES_size_sonicnozzlearea.m` | superseded by `_air` | Differs only in `A_feed`, which is unused. |
| `HADES_size_sonicnozzlearea_h2.m` | superseded by `_air` | Byte-identical copy of the air script; it still uses air properties. |
| `untitled.m` | dropped | Copy of `HADES_size_sonicnozzlearea.m` with git merge-conflict markers. |
| `cell_size2.m` (+ `.asv`) | superseded by `HADES_size_annulusgap_fillheight.py` | Its unique calculation (cell width from the Caltech-data fit at phi) is kept as `cellwidth_from_caltech_fit(phi)` and validated against Octave. |
| `testing.m` | dropped | Scratch; calls `HoopStressTemps` with the old 4-argument signature. |
| `TEST IGNORE.m`, `list of programs` | dropped | Junk text and an empty file. |
| `*.asv` (convert, airInj, h2Inj, HoopStressTemps, annulusgap_fillheight, throat, cell_size2) | dropped | MATLAB autosaves of older versions. |
| `CEA/FCEA2.exe`, `thermo.*`, `trans.*`, `cea2.inp`, `run.txt` | replaced by the `cea` package | It ships its own thermo/trans databases. |
| `CEA/cea_det.*`, `cea_run.*`, `h2oPHI1.inp`, `h2ophi1.out` | `tests/hades_sizing/reference/` | Kept as validation fixtures. |

## 2. CEA replacement

**Choice: NASA's official open-source CEA re-implementation, `cea` on PyPI (v3.3.4, Apache-2.0,
https://github.com/nasa/cea).** It has wheels for macOS (arm64/x86_64), Linux and Windows (CPython 3.11-3.14), and
solves CJ detonation, rocket (IAC/FAC) and fixed-state equilibrium problems in-process, with the same
NASA Glenn thermo and transport databases. It was added to `pyproject.toml` as `cea>=3.3.4`.

Options considered:
- `rocketcea`: wraps the CEA2 Fortran, but PyPI has only Windows wheels (elsewhere it needs a gfortran build), and it has no detonation problem.
- Cantera plus an SDToolbox-style CJ solver: pip-installable, but it would mean writing and maintaining our own CJ solver and NASA-9 data. Agreement with CEA2 would be approximate rather than exact.
- `pycea`: deprecated in favour of the official `cea` package.
- The CLI bundled in the `cea` wheel cannot run the HADES decks as written. It rejects the empty `only` dataset, and reactants with no `t=` make it return NaN. It also prints in a different format, so the MATLAB parser would not work on its output.

`rdelib/cea.py` wraps the Python API as `cj_detonation`, `rocket` and `equilibrium` (hp/tp). The HADES wrappers reproduce
the MATLAB problem exactly:
- **Unit conversions:** the MATLAB's own conversions are kept.
- **Deck formatting:** values are rounded to the 6 decimals of the deck's `%f` format.
- **Pressure units:** CEA2's internal psia->bar factor (1.01325/14.696006) is applied.
- **Reactant enthalpy:** zero for the rocket problem, because the decks gave no reactant temperature.
- **Output rounding:** outputs are rounded to the precision CEA2 printed in the `.out` file, which is what the MATLAB parser returned (`cea2_rounding=True`, the default). Full-precision results are always in `.raw`.

The decimals follow CEA2's formats as observed in the committed outputs; quantities that vary by magnitude use CEA2's
VARFMT rule. P/P1 and T/T1 use fixed 3 decimals; the review confirmed this against the legacy CEA2 source (`cea2.f`, line 611).
The plain `rdelib.cea` functions always return full precision; only the legacy HADES wrappers round.

`cea2_input_deck(...)` in both wrappers renders the exact deck the MATLAB wrote, which is useful for anyone still running CEA2.

**Agreement achieved** (`tests/hades_sizing/test_rdelib_cea.py`): every value printed in `cea_det.out`, `cea_run.out`
and `h2ophi1.out` is reproduced to the printed precision (within half a unit of the last printed digit). This covers
unburned/burned states, transport properties, detonation parameters, all three rocket stations, the performance
parameters and the mole fractions. Exceptions:
- Detonation equilibrium conductivity: 5.21016 vs printed 5.2103 (2.7e-5 relative).
- Rocket throat T: 2094.91497 K vs printed 2094.92. The value sits on the rounding boundary, so the difference is at least 3e-5 K (about 1e-8 relative).
- Rocket chamber H: 1.4e-4 kJ/kg vs exactly 0. This is solver residual.
- `h2ophi1` (hp, products `only H2O`): CEA3 refuses a product list shorter than the reactant list when reactants are attached (it prints CRITICAL and then segfaults). CEA2 reduced the components instead. `rdelib.cea.equilibrium` now handles this in the public API, and the reference test goes through it. It solves the physically equivalent problem: product amounts with the same elemental composition, and the reactant enthalpy as the hp target. T = 4930.56 K, and every printed property matches.

**Guarding the native solver.** Bad restricted product sets can crash or silently mislead the native library. When `only` or `omit` is given, `rdelib.cea` validates the elemental feasibility in Python *before* any native call:
- Element compositions are parsed from CEA species names. Reactant element totals come from a standard full-product equilibrium, since elements are conserved.
- It raises `ValueError` if a product contains an element absent from the reactants (native: CRITICAL, then crash); if a reactant element appears in no product (native: silently "converges" to 298 K garbage); or if the product set cannot reproduce the elemental proportions (e.g. `only=['H2O']` off stoichiometric).
- It raises `NotImplementedError` for detonation or rocket problems with fewer products than reactants.
- Unknown species names raise `RuntimeError` from `cea.Mixture`.
- All of this is tested in a subprocess, so a regression cannot take down the test run.
- Sanity checks: stoichiometric H2-air CJ at 1 atm / 298 K gives about 1968 m/s, and stoichiometric H2/O2 hp gives 3050-3110 K.

### MODELING copy of `HADES_size_ceaDet.m`
`MODELING/.../RDE_Toolbox_2/src/HADES_size_ceaDet.m` differs from the HADES copy in one line only: it uses `fopen(inputFile)`
(read-only) instead of `fopen(inputFile,'w')`. That copy cannot write the CEA deck, so it is a regression rather than an
improvement, and nothing was carried over. Every other file of `HADES/SIZING MODEL` (including the CEA folder) exists in
that folder byte-identically. The MoC port should import `HADES.sizing_model.*`.

## 3. Numerical-method equivalences

- **`ode45`:** `rdelib.ode.ode45` is a line-by-line reimplementation of MATLAB's algorithm: Dormand-Prince coefficients, initial-step heuristic, infinity-norm error test with threshold AbsTol/RelTol, step update, and Refine=4 output with the ntrp45 interpolant. It exists for MATLAB sample-grid compatibility: the cooldown scripts report the *first output sample* below a target, so the output grid matters. For new calculations, prefer `scipy.integrate.solve_ivp` with `events`.
  - **Input checks:** finite, strictly monotonic `tspan`; finite `y0`; positive tolerances, `max_step` and `initial_step`; correct derivative shape.
  - **Minimum step:** `eps(t)` is `np.spacing(abs(t))`. MATLAB's `eps` is always positive, while `np.spacing` of a negative t is negative and previously hung.
  - **Failures raise instead of hanging or returning NaN/inf:** non-finite derivatives, error estimates, or stage, accepted and interpolated states (e.g. overflow to inf with a finite derivative); step-size underflow; and failure to meet the tolerances. MATLAB would carry on silently.
- **`pdepe`:** `rdelib.ode.pdepe` implements the Skeel-Berzins discretisation MATLAB uses (m = 0 or 1, scalar PDE, flux/Robin boundaries) and integrates in time with scipy BDF.
  - Default tolerances are rtol 1e-6 / atol 1e-8 rather than MATLAB's 1e-3 / 1e-6. scipy BDF at 1e-3 is 0.06-0.08 % off the converged time-to-limit, while Octave's `ode15s` at 1e-3 is within 0.02 %, so the tighter setting lands closer to MATLAB's answer.
  - Validated against the analytic steady state (cylinder, Robin-Robin), the analytic convective semi-infinite solid (m = 0), mesh refinement, and an independent Octave `ode15s` implementation.
- **`fzero(f, x0)`:** reproduces MATLAB's bracket search (dx = x0/50, growing by sqrt(2)), which decides *which* root is found, and then refines with Brent's method. Bracketed `fzero` maps to `brentq`.
- **`fsolve`:** `scipy.optimize.fsolve` from the same start points. Differences from MATLAB/Octave are around 1e-8 relative (the test tolerance for `P_plenum_fanno` is 1e-7).
- **Colebrook:** `rdelib.gasdynamics.colebrook` used `fsolve` from f = 0.02 until the final cleanup; it now shares DEIMOS's bracketed solve in x = 1/sqrt(f) (`colebrook_darcy`). The Colebrook root is unique, so this only removes fsolve's stopping error: on the main inputs f changes by at most 1.6e-15 relative and `P_plenum_req` by 1 ulp (1954432.758810425 -> 1954432.7588104247 Pa, -1.2e-16); the printed report and the second Octave case are unchanged. On a 1200-point grid, `Re = np.geomspace(3e3, 1e8, 200)` x `eps/D` in {1e-6, 1e-4, 1.25e-3, 2.7e-3, 0.0167, 0.05}, the old fsolve stalled at its 0.02 start point with a RuntimeWarning ("invalid value encountered in sqrt") in 54 cases; elsewhere the two agree to 6e-14, and the new root's Colebrook residual is below 1e-12 everywhere. (An independent review on a 40 x 30 log grid found 79/1200 stalls, max new residual 5.3e-15.)
- **Other mappings:** `a:b:c` -> `colon` (MATLAB element count and symmetric generation), `polyfit`/`polyval` -> numpy (same least squares), `interp1` -> `np.interp` (`T_inner` is monotonic).
- **Plots:** every figure is reproduced with matplotlib. Functions take `plot=True` and never call `plt.show()`; only `main()` does.

## 4. Station ambiguity (4a) and MATLAB bugs fixed (4b, 4c; before -> after)

**a. Rocket Isp station: ambiguity flagged, legacy preserved, exit option added (not treated as a bug).**
The MATLAB parser takes the first number on each PERFORMANCE PARAMETERS line, i.e. the **throat** column: `AeAt` = 1,
`isp` = 927.8 m/s, `cf` = 0.6884, `ivac` = 1679.1 m/s. The deck's `supar=10` is only the wrapper's fallback, and main supplies no
nozzle geometry. Using the exit column would also be wrong on its own: mdot x Isp_exit leaves out the (pe - pa)Ae term, and here pe is
0.0984 bar against an ambient of 1.013 bar. Per the orchestrator's decision:
- `HADES_size_ceaRocket` keeps the legacy fields `AeAt`, `cf`, `isp`, `ivac` with the MATLAB (throat) meaning.
- It also exposes both stations unambiguously: `*_throat` (same as legacy), `AeAt_exit`, `cf_exit`, `isp_exit`, `ivac_exit`, `pe_exit_bar` and `pc_bar`.
- `HADES_size_main` is unchanged: thrust 863.770 N (194.176 lbf). The printed report is identical to the MATLAB, line for line.
- Exit-plane thrust is an explicit opt-in: `HADES_size_thrust(Ve, mdot_lb, station="exit", exit_pressure=..., ambient_pressure=..., exit_area=...)`, which applies the NASA rocket thrust equation F = mdot Ve + (pe - p0) Ae, or the convenience wrapper `exit_thrust_from_cea(ceaRock_results, mdot_lb, p0)`.
  - Validated by hand calculation from the printed `cea_run.out` exit column (within 2e-4).
  - Validated against the CEA identity Ivac = Ve + pe Ae/mdot: with p0 = 0, F = mdot x Ivac to 1e-6.
  - For the main case at sea level (eps = 10, over-expanded) it gives about 843.6 N, versus 2010 N from mdot x Isp_exit.

**b. `HADES_size_convert` 'atm' factor was inverted.** The table holds "units per bar", but 'atm' had 1.01325 where
0.986923 belongs. Before: `convert(1,'atm','bar')` = 0.98692. After: 1.01325. Not used by the main chain, so no change to the main output.

**c. `HADES_size_sonicnozzlearea_air`: `R = 8.134` is a typo for 8.314 J/(mol K).** Choke diameter rises 0.55 %:
- 20 bar: 16.0885 -> 16.1768 mm
- 40 bar: 11.3763 -> 11.4387 mm
- 60 bar: 9.2887 -> 9.3397 mm
- The 11 mm minimum-diameter line is crossed at 43.26 bar instead of 42.79 bar.

`HADES_size_sonicnozzlearea_air(R=8.134)` reproduces the original.

## 5. Questionable but left unchanged (physics as written)
- **`HADES_size_max_flux_range`:** the label says kg/(s*mm^2), but main passes lbm/s, so the value is in lbm/(s mm^2).
- **`HADES_size_main` fill-height inputs:** CJ speed 1931.5 m/s, a2 1075.6 m/s and P2 32.33 are hand inputs. They are not the CEA results (1929.7 m/s, 1074.2 m/s, 28.927 bar).
- **`HADES_size_annulusgap_fillheight`:** returns the fixed design gap of 9.525 mm. The 1.4 x cell-width correlation (8.917 mm) is only printed, and the data fit only plotted.
- **`HADES_size_HoopStressTemps`:** the failure temperature for the main inputs (902.9 C = 1657 F) comes from extrapolating the cubic fit beyond the data, whose maximum is 1500 F. MATLAB and the port agree.
- **`HADES_size_h2Inj`:** `density = 0.91` is flagged "UNSURE" in the source.
- **`HADES_size_P_plenum_fanno`:** the Sutherland viscosity mixes T and T_stag, and Re uses the guessed Mach number. Both kept as written.
- **Sonic nozzle:** `c_f = 0.6989`, whereas the ideal critical-flow factor for gamma = 1.4 is 0.6847.
- **`HADES_size_thermal_v1` CEA property table:** the comment says each phi uses a different CEA column. The CJ temperatures and viscosities match CEA at T1 = 298.15 K and P1 of about 2 bar for all four phi, but the densities do not, so the table cannot be regenerated from a single condition.
- **`HADES_size_convert`:** inch factor 39.37 per m (not 39.3701) kept.

## 6. Validation evidence

GNU Octave 11.3.0 ran the **original** `.m` files. The harness and its notes are in `tests/hades_sizing/reference/octave_harness/`.
Shims cover the functions Octave lacks (`readmatrix`, `contains`, `splitlines`, `optimoptions`, `fit`/`fittype`, `pdepe`). Three
Octave-only patches were applied to /tmp copies, never to the originals:
- `inputParser` PartialMatching off, working around an Octave bug where `'P0'` matched `'P0Units'`;
- string `+` concatenation in `h2Inj` replaced;
- script-local functions moved to separate files.

`FCEA2.exe` was replaced by a shim that replays the committed CEA2 `.out` files, and only if the deck MATLAB just wrote is byte-identical
to the committed deck. It was, for both decks.

| Comparison | Result / tolerance |
|---|---|
| `HADES_size_main` vs Octave run of `HADES_size_main_v1.m` | All scalars including thrust rel 1e-12 (`P_plenum` 1e-7); CEA struct fields identical (k_eq within 1e-4). Full stdout identical line by line, including Thrust 863.770 N. Chamber pressure 142.7451485 -> `%f` = `142.745149`, as in the committed `cea_run.inp`. |
| Regenerated decks | Byte-identical to committed `cea_det.inp` / `cea_run.inp`. |
| Individual functions (throat, pdrop_fanno x2, HoopStressTemps x2, P_plenum_fanno, fill height, chamberPressure, h2Inj, airInj, geometry, thrust, P0, max_flux, convert x14, sonic sweep, cell_size2 fit) | rel 1e-12 / 1e-10 / 1e-7 (algebra / fzero / fsolve). |
| Cooldown v1 and v3 vs Octave `ode45` | Cooling times within 1.5 s (output spacing about 3 s; Octave's ode45 is not MATLAB's); final wall T within 1e-4 K. |
| Thermal v1 vs Octave (`pdepe` shim + `ode15s`) | h, T_rec, Re, Nu rel 1e-12; time-to-limit [0.04155, 0.14186, 0.98290] s vs [0.04155, 0.14184, 0.98268] s, rel 5e-4; T_inner within 0.5 C. |
| `rdelib` helpers | Hand and table values (isentropic, Fanno, Moody), analytic ODE/PDE solutions, `ode45` failure paths, native-crash guards for CEA. |
| Packaging (`test_packaging.py`) | `uv build` wheel built from a pristine temporary copy of the source (no `build/`, `*.egg-info`, `.venv` or `__pycache__`). The wheel must contain the sizing xlsx, the feed-sizing and data-analysis CSVs, the BOS assets and the MoC LICENSE, and must contain no legacy path (`SIZING MODEL`/`FEED SIZING`/`DATA ANALYSIS`/`MODELING_size_*`, `*.m`, `*.exe`, `*.asv`). The legacy folders are deleted from the tree, and `namespaces = false` in `[tool.setuptools.packages.find]` means only folders with an `__init__.py` become packages, so an extracted copy of the originals inside the checkout would not be packaged either; the test asserts both lists on the wheel contents. The wheel is then installed `--no-deps` into a fresh venv outside the checkout; `HADES_size_main` runs and reads the bundled xlsx (thrust and Pc rel 1e-12). `pyproject.toml` declares package data (xlsx/csv/json/svg/png/md/LICENSE, plus `HADES.data_analysis/bos_print_assets/*`) for all areas. |
