# HADES feed sizing: MATLAB/Python to Python port notes

Source folder: `HADES/FEED SIZING/`. Target package: `HADES/feed_sizing/`.

## Where each original went

The originals were scripts with no function signatures, so the Python functions have descriptive names. There are no alias shims. Run the CLI from the repo root with `uv run python -m <module> [subcommand]`.

| Original script | Python function | CLI |
| --- | --- | --- |
| `SonicNozzleSizing.py` (prototype) | `sonic_nozzle_sizing.size_sonic_nozzle_with_feed_line` | `HADES.feed_sizing.sonic_nozzle_sizing feedline` |
| `SonicNozzleSizingV1.m` | `sonic_nozzle_sizing.size_single_sonic_nozzle` | `HADES.feed_sizing.sonic_nozzle_sizing v1` |
| `SonicNozzleSizingV2.m` | `sonic_nozzle_sizing.size_fixed_sonic_orifice` | `HADES.feed_sizing.sonic_nozzle_sizing v2` |
| `SonicNozzleSizingV3.m` | `sonic_nozzle_sizing.size_two_stage_choked` | `HADES.feed_sizing.sonic_nozzle_sizing v3` |
| `SonicNozzleSizingV4.m` | `sonic_nozzle_sizing.manifold_pressure_from_injector` | `HADES.feed_sizing.sonic_nozzle_sizing v4` |
| `SonicNozzleSizingV5.m` | `sonic_nozzle_sizing.size_injector_sonic_orifice` | `HADES.feed_sizing.sonic_nozzle_sizing` (default `v5`) |
| `blowdown_orifice_cd.py` | `blowdown_orifice_cd.blowdown_orifice_cd` (`Cds`, `average_Cd`), `plot_cd_vs_time` | `HADES.feed_sizing.blowdown_orifice_cd` |
| `blowdown_orifice_cd_chat.py` | `blowdown_orifice_cd.blowdown_orifice_cd` (`Cd`, `Cd_integrated`), `plot_mass_calibration` | `HADES.feed_sizing.blowdown_orifice_cd` |
| `Orifice_CD_Results.m` | dropped (see below) | none |

## File mapping

| Original | Python | Notes |
| --- | --- | --- |
| `SonicNozzleSizing.py` | `size_sonic_nozzle_with_feed_line` (CLI `feedline`) | The Feb-2026 first prototype, kept as a separately labelled workflow because it is the only one where the throat depends on the feed line: it computes the feed-line Mach number, converts static to total temperature, and sizes the throat from the isentropic area-Mach relation. Its area-Mach bug is corrected; see issue 1. Later history: it was translated to `SonicNozzleSizing.m` (1db24d9), the bug was fixed there (f8c4cd3), and the file was then rewritten into the V1 approach (6bc6d88). |
| `SonicNozzleSizingV1.m` | `size_single_sonic_nozzle` (CLI `v1`) | Single-point air nozzle: flow margin, throat static state, feed-line ID. |
| `SonicNozzleSizingV2.m` | `size_fixed_sonic_orifice` (CLI `v2`) | One fixed orifice per gas, fixed downstream pressure per gas. |
| `SonicNozzleSizingV3.m` | `size_two_stage_choked` (CLI `v3`) | Injector area sized from chamber pressure, then the sonic orifice. |
| `SonicNozzleSizingV4.m` | `manifold_pressure_from_injector` (CLI `v4`) | Gives the same manifold pressure and choke flags as V5 (V5 is a superset). Kept as its own short function because it writes `manifold_pressure_results.csv` and has no positive-flow filter. |
| `SonicNozzleSizingV5.m` | `size_injector_sonic_orifice` (CLI default `v5`) | Canonical workflow. |
| inline choked-flow equations (all of the above) | `choked_flow.py` | Shared helpers; good candidates to move into `rdelib`. |
| `blowdown_orifice_cd.py` | `blowdown_orifice_cd.py` | Merged with the chat variant. |
| `blowdown_orifice_cd_chat.py` | `blowdown_orifice_cd.py` | Same preprocessing and ideal-flow model as `blowdown_orifice_cd.py`, with a second Cd estimator added, so it is a variant. |
| `Orifice_CD_Results.m` | dropped | Header says "WASTE PROGRAM". It is unfinished: it reads the CSV, sets `start_pressure`/`end_pressure` and stops, with no output except `gas = `. That was confirmed by running it in Octave. It also mixes units (`end_pressure` is never converted to Pa). |
| `flow_requirements.csv`, `flow_requirements_manifold_input.csv`, `orifice_blowdown_input.csv`, `orifice_blowdown_input_template.csv` | copied into this folder | The bundled default inputs, resolved relative to this module. |
| `flow_requirements_sonic_orifice_results.csv`, `flow_requirements_two_stage_choked_results.csv`, `manifold_pressure_results.csv`, `flow_requirements_manifold_input_combined_results.csv` | copied into this folder; reference copies in `tests/hades_feed_sizing/reference/matlab/` | The committed MATLAB results. The Python port regenerates them byte-for-byte, into the CWD. |
| `orifice_blowdown_results.csv` | `tests/hades_feed_sizing/reference/original_python/` only | Stale. It was written by the c7402e4 version of `blowdown_orifice_cd.py` (real-gas sonic state, `np.gradient` mass flow, incompressible Cd), which was replaced in 6c5df96. No current code produces it, so it is only used as a cross-check of the pressures and tank mass. |

## Behaviour changes (no physics changed except issue 1)

- **File locations.** Bundled default inputs resolve relative to the module.
  - Generated result CSVs go to the current working directory by default, named `<input name><suffix>.csv`, or `manifold_pressure_results.csv` for V4. `--output` (or `run(..., output_file=...)`) picks another path.
  - So a run never overwrites the bundled reference CSVs and does not need a writable package folder.
  - MATLAB, run on its default bare-filename input, also wrote into `pwd`. With an input in another folder, MATLAB wrote next to the input; Python still writes to the CWD.
  - The "Results written to" line prints the absolute path.
  - `blowdown_orifice_cd` writes no files; it only prints and plots.
- **CSV parsing.** CSVs are read with `float_precision="round_trip"`, which parses numbers with correct rounding like MATLAB `readtable`. Without it, pandas' fast parser differed from MATLAB by up to 1 ulp in the input columns.
- **Blowdown input checks.** The original scripts called CoolProp one sample at a time, so an invalid state raised `ValueError`, e.g. `For now, we don't support T [nan K] below Ttriple`. The port calls CoolProp vectorized, which returns `inf` for an invalid sample instead of raising.
  - Before computing Cd, the port therefore rejects non-finite `time_ms` and a run duration of 0 or less.
  - It also rejects a supply pressure or supply temperature that is non-finite or not above zero (absolute).
  - It also rejects any `Dmass`, `Cpmass` or `Cvmass` that is non-finite or not above zero, and any `gamma` not above 1.
  - Each case raises `ValueError` naming the 0-based data row indexes (CSV line = index + 2) and the CoolProp reason. The CLI exits with status 1 and prints the message.
  - If no row at all can be evaluated, vectorized `PropsSI` raises `No outputs were able to be calculated`. The port then falls back to per-row scalar calls, so the failing rows and CoolProp's reason are still reported.
  - Non-numeric cells are coerced to NaN, so they are rejected too.
  - Upstream and downstream pressures are only converted, not used for Cd, so they are not checked (same as the originals).
- **Blowdown output.** The merged script prints the four lines of `blowdown_orifice_cd_chat.py` plus `Average Cd (run-average mdot / ideal mdot)`, a value the old `blowdown_orifice_cd.py` only showed in its plot legend. Both original plots are now separate figures: in the old scripts both used `plt.figure(1)`.
- **Speed-up.** CoolProp property calls are vectorized. Valid-data results are bit-identical to the original per-sample loops.
- **Warnings.** `warning(...)` became `warnings.warn(..., RuntimeWarning)` with the same message text. The feed-line workflow also warns if the feed-line Mach number is 1 or more.

## Issues found

1. **Bug fixed: area-Mach relation in `SonicNozzleSizing.py`.** The prototype had `2/(gamma+1)` outside the power:
   - before: `A_factor = 2/(g+1) * (1 + (g-1)/2*M^2)**(0.5*(g+1)/(g-1))`
   - after: `A_factor = (2/(g+1) * (1 + (g-1)/2*M^2))**(0.5*(g+1)/(g-1))`
   - Both feed `A_t = A_2 / A_factor * M_2`. With the correct form, `A_2/A_t` equals the isentropic `A/A*` at `M_2`.

   The corrected form is exactly what the authors wrote in the MATLAB translation (f8c4cd3). Results on the prototype's default inputs (3 lb/s, `P_m` = 20 bar, wall 0.095 in, so ID = 0.810 in):

   | | before (original `.py`) | after (port) |
   | --- | --- | --- |
   | M_2 | 0.419004 | 0.419004 |
   | A_factor | 0.924234 | 0.641829 |
   | A_t [m^2] | 1.507175e-4 | 2.170333e-4 |
   | D_t [mm] | 13.852781 | 16.623337 |
   | P_t,u [bar] | 38.232314 | 26.550218 |
   | legacy threshold `P_m/crit_ratio` [bar] | 37.858583 | 37.858583 |
   | legacy pressure-ratio check | "Flow is choked!" (passed) | not passed |

   With a 0.049 in wall, D_t goes from 14.102979 to 16.923575 mm. Other prototype quirks are kept as written:
   - `mass_flow_margin = 0.85` only scales the flow used for `M_2`.
   - `P_t,u` uses the unscaled flow with no Cd.
   - `R = 287`.

   **The old "Flow is choked!" verdict is now a legacy diagnostic, not a choking determination.** The historical calculation is kept as `legacy_pressure_check_passed = P_t,u > P_m/crit_ratio`, and its threshold as `P_t_u_req`.
   - The printout labels it "Legacy pressure-ratio check" and prints a NOTE on every run saying it is not a valid choking criterion.
   - The old field `choked` and the old verdict strings are gone.
   - No replacement choking criterion is provided.

   Why it is invalid: the model already puts the throat at exactly M = 1, with isentropic recovery to the feed-line static pressure `P_m` downstream. The legacy check compares `P_m` with the throat critical pressure `crit_ratio * P_t,u`, which ignores that recovery.

   Counterexample (`mass_flow_margin = 1`, other defaults):
   - M_2 = 0.492946.
   - `P_t,u` = 23.613659 bar, equal to the isentropic total pressure of the feed-line state, `P_m (1 + 0.2 M_2^2)^3.5`, to 1e-13.
   - `A_t` = A* from independent sonic-state continuity, `mdot sqrt(T_t) / (P_t,u K)`, to 1e-13.
   - The throat static pressure is 12.474666 bar, below `P_m` = 20 bar (pressure recovery).
   - So the nozzle is exactly sonic by construction, yet the legacy check reports "not passed" (23.61 < 37.86 bar).

   With the default margin of 0.85, `P_t,u` is 1/0.85 = 1.176471 times the isentropic total pressure of the feed-line state (26.550218 vs 22.567686 bar), because the scaled and unscaled flows are mixed.
2. **No numerical bugs in V1-V5.** The Python reproduces them bit-for-bit. V1's comment says `P_down ... = 20 bar`, but the code adds a 20 psi margin: `20e5 + 20*psi_to_Pa`, which is 21.38 bar. The code is kept as written.
3. **Inconsistent psi to Pa constants**, all kept per script:
   - 6894.76 in V1-V3
   - 6894.757293 in V5
   - 6894.75729 in the blowdown scripts
   - 6894.757293168 in the old blowdown version

   This causes relative differences of about 4e-7 in psi outputs between versions. Recommend using one constant from `rdelib` later; doing so changes the 7th significant digit.
4. **Blowdown, `average_Cd` method.** It divides one run-average mass flow (total tank mass loss / duration) by the instantaneous ideal flow, so the spread of `Cds` over time (0.534924 to 0.585122 on the committed data) is not a real Cd history.
5. **Blowdown, the two Cd estimators disagree.** On the committed input, the least-squares `Cd` is 0.640279 and `Cd_integrated` is 0.564165.
   - The measured mass loss is not proportional to the ideal choked prediction. The measured tank outflow (the current supply-based tank mass, `-d m/dt` by `np.gradient`) falls from 0.697225 to 0.199661 kg/s.
   - The current ideal choked flow (supply P, T; Cd = 1) only falls from 0.748208 to 0.684018 kg/s.
   - So the cumulative ratio drifts: 0.884 at 100 ms, 0.719 at 500 ms, 0.564 at the end.
   - For comparison, the stale `orifice_blowdown_results.csv` evaluated its ideal flow at UPSTREAM conditions: ideal-gas 0.741611 to 0.677685 kg/s, real-gas 0.739866 to 0.676820 kg/s. Its gradient mass flow, 0.697225 to 0.199661 kg/s, matches the current one.
   - The committed input looks synthetic, so neither Cd should be used for hardware until real data is run.
6. **Blowdown, real-gas gamma in an ideal-gas formula.** The ideal-gas choked formula is given a real-gas `gamma = cp/cv` of 1.584 to 1.586 at about 10 MPa. In the stale results CSV (upstream conditions, t = 0), this gives 25409.5 kg/(m^2 s), against 25349.7 for the real-gas isentropic sonic flux: 0.24 % high. Using gamma = 1.4 instead would be about 4 % low. Kept as written.

## Validation

Ran `uv run pytest tests/hades_feed_sizing`.

- **V1-V5 vs GNU Octave 11.3.0.** The original `.m` files were copied unmodified to `/tmp`. Octave has no `table`, `readtable` or `writetable`, so they were replaced with minimal shims: a classdef stand-in supporting `.Properties.VariableNames`, `data{:,i}`, `data.X = ...` and `height`. Results were written at 17 significant digits.
  - Octave also reproduces the committed MATLAB CSVs to within 4.3e-15, the 15-digit rounding.
  - The Python port is bit-identical to Octave in every number of all four result tables and the V1 scalars. The tests allow rtol = 1e-14.
  - Python's printed summaries are text-identical to Octave's command-window output for V1-V5.
- **V2-V5 vs the committed MATLAB result CSVs.** Column order, NaN pattern and the 1/0 flags are identical, numbers agree within rtol = 1e-13, and the regenerated files are byte-identical.
- **Feed-line workflow.**
  - The authors' corrected MATLAB (`git show f8c4cd3:"HADES/FEED SIZING/SonicNozzleSizing.m"`, 2 lb/s) was run in Octave. The port called with `m_dot = 2 lb/s` is bit-identical in all 16 variables (`reference/octave/SonicNozzleSizing_f8c4cd3_values.csv`).
  - Independent checks:
    - `A_2/A_t` matches textbook isentropic `A/A*` table values for gamma = 1.4 (M = 0.1 to 0.8, to 1e-4) and the closed-form `A/A*(M_2)` (1e-13).
    - A brentq inversion of `A/A*` recovers `M_2`.
    - `M_2` matches continuity, and `T_t` and `P_t,u` match their textbook forms.
  - The buggy formula, rebuilt in the test, reproduces the original prototype's printed output (`reference/original_python/SonicNozzleSizing_py_stdout.txt`), and the buggy/fixed ratios are checked against the table above.
  - The `mass_flow_margin = 1` counterexample to the legacy pressure-ratio check (issue 1) is a regression test.
- **Hand calculations.** Textbook choked-flow form for V1. Round-trip mass flow through the sized orifices (V2). The sizing point equals the design minimum pressure. The two-stage pressure ratio stays constant at crit/1.15 (V3). V4 and V5 give the same manifold pressure. V5 choked/valid counts and warnings checked.
- **Blowdown.** The originals were run in `/tmp` with the Agg backend (CoolProp 8.0.0). Every per-sample array and scalar is bit-identical except `average_Cd`, which differs by 1 ulp; the tests allow rtol = 1e-9 for CoolProp version drift. Printed lines match the chat script.
  - Against the stale `orifice_blowdown_results.csv`: pressures, tank mass and the run-average mass flow agree to 4.6e-10, which is exactly the psi to Pa constant difference.
  - Independent checks: textbook choked flow, `np.trapezoid`, and `np.linalg.lstsq` for the least-squares Cd.
  - Invalid samples are covered by tests: an interior `supply_temperature_K` of NaN, 0, -1, 40 K or blank, and a supply pressure below vacuum. Each raises `ValueError` naming row 5 (or the bad row), and the CLI exits non-zero.
  - An all-invalid batch (every supply temperature 40 K) still reports rows [0, 1, ...] (1001 rows) and the `Tmelt` reason.

### Reproducing the Octave comparison

1. Copy the `.m` files and the input CSVs into a scratch folder.
2. Put shims for `readtable` (parse with `strsplit(..., ',', 'CollapseDelimiters', false)` so blank cells become NaN), `writetable` (`%.17g`, `NaN`, logical as `%d`) and a small table class on the Octave path.
3. Run `octave --no-gui --eval "addpath('<shims>'); SonicNozzleSizingV5"`.

The reference outputs are in `tests/hades_feed_sizing/reference/octave/`.

## Candidates for `rdelib`

These all live in `choked_flow.py`: `critical_pressure_ratio`, `choked_flow_parameter`, `choked_mass_flow`, `choked_area`, `choked_upstream_pressure`, `area_to_diameter`. The unit constants (`IN_TO_M`, `LB_TO_KG`, the `PSI_TO_PA` variants, `L_TO_M3`, `P_ATM_PSI`) could move too. Two expressions are written inline to keep the original operation order: the blowdown ideal flow (the same as `choked_mass_flow(1, A, P0, T0, gamma, R)`) and the feed-line area-Mach factor.

`rdelib/gasdynamics.py` appeared while this port was being finished. Its equivalents are `critical_pressure_ratio`, `choked_flow_factor`, `choked_throat_area`, `choked_mass_flow` and `area_ratio`; `rdelib/units.py` has `PSI_TO_PA = 6894.76`, `M_PER_IN` and `LBM_TO_KG`. The package was not switched over yet. When it is, re-run `tests/hades_feed_sizing`: the Octave comparison uses rtol = 1e-14, so a different operation order is fine, but a different psi constant is not.
