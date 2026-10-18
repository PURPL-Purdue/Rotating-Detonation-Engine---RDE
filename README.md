# PURPL Rotating Detonation Engine (RDE)

Design and analysis tools for PURPL's (Purdue) rotating detonation engine work:
the HADES engine sizing model, its feed-system sizing and test-data analysis,
the DEIMOS orifice and Fanno-flow calculators, and a method-of-characteristics
(MoC) nozzle solver with the team's RDE MoC model. Everything is Python. It was
ported from the team's MATLAB, and each port was checked number-for-number
against the original run in GNU Octave. The MATLAB originals were then removed
from the working tree (they are still in git history, see
[Original MATLAB](#original-matlab-and-live-octave-checks)). Shared physics and
MATLAB-compatibility helpers live in `rdelib/`.

## Setup

Install [uv](https://docs.astral.sh/uv/), then from the repository root run:

```bash
uv sync            # creates .venv with Python 3.12 and the locked dependencies
```

NASA CEA comes from the `cea` package, which has wheels for macOS, Linux and
Windows. The Windows-only `FCEA2.exe` is no longer needed.

## Running the tools

Run the commands from the repository root. To run them from another folder,
use `uv run --project /path/to/repo python -m ...`. Output files (result CSVs,
BOS images, synthetic data) go to the current directory unless an option says
otherwise. Bundled default inputs are read from the package, so they are found
from any folder. Figures open in a window; set `MPLBACKEND=Agg` to run without
one.

| Tool | Command |
|---|---|
| HADES sizing model (full chain, prints the report) | `uv run python -m HADES.sizing_model.HADES_size_main` |
| Chamber cooldown, two-node model (v1) | `uv run python -m HADES.sizing_model.HADES_size_chambCooltime` |
| Chamber cooldown strategies (v3, four scenarios) | `uv run python -m HADES.sizing_model.HADES_size_chambCooltime_v3` |
| Hot-fire wall-temperature limit (`pdepe`) | `uv run python -m HADES.sizing_model.HADES_size_thermal_v1` |
| Air sonic-nozzle choke-diameter sweep | `uv run python -m HADES.sizing_model.HADES_size_sonicnozzlearea_air` |
| Feed sizing, canonical V5 (injector + sonic orifice) | `uv run python -m HADES.feed_sizing.sonic_nozzle_sizing` (same as `... v5`) |
| Feed sizing, other workflows | `uv run python -m HADES.feed_sizing.sonic_nozzle_sizing v1` (or `v2`, `v3`, `v4`, `feedline`) |
| Feed sizing, own inputs | `uv run python -m HADES.feed_sizing.sonic_nozzle_sizing v3 --input my.csv --output out.csv` (or `--no-write`) |
| Blowdown orifice Cd | `uv run python -m HADES.feed_sizing.blowdown_orifice_cd [--input test.csv --gas Air --tank-volume-L 50 --orifice-diameter-in 0.24 --no-show]` |
| Heat flux from thermocouples (model choice **required**) | `uv run python -m HADES.data_analysis.heat_flux_calc --equations cubic --depths-m 0.25e-3 0.575e-3` |
| Heat flux, original equations / depths | `uv run python -m HADES.data_analysis.heat_flux_calc --equations legacy --matlab-depths` (warns: known-invalid, units unconfirmed) |
| Synthetic ITP pressure data | `uv run python -m HADES.data_analysis.purpl_fake_data [--output-dir DIR \| --no-csv]` |
| ITP FFT (v2 averaged spectrum; `--raw` = April pipeline) | `uv run python -m HADES.data_analysis.purpl_fft [--raw]` |
| Mass flow from plenum pressure | `uv run python -m HADES.data_analysis.purpl_mass_flow_from_plenum_pressure` |
| Choked-nozzle Cd / temperature sweeps | `uv run python -m HADES.data_analysis.choked_nozzle_analysis [cd\|temp]` |
| CDR design-point conditions | `uv run python -m HADES.data_analysis.cdr_conditions` |
| Upstream (plenum) conditions | `uv run python -m HADES.data_analysis.upstream_conditions` |
| SCFM conversion | `uv run python -m HADES.data_analysis.scfm` |
| BOS random-dot background (SVG + 600 dpi PNG) | `uv run python -m HADES.data_analysis.bos_background_generation [--output-dir DIR] [--no-plot]` |
| DEIMOS orifice Cd | `uv run python -m DEIMOS.DEIMOS_sys_CdCalc_v1` |
| DEIMOS Fanno-flow calculator | `uv run python -m DEIMOS.DEIMOS_sys_FannoFlowCalc_v2` |
| MoC nozzle + plume (team inputs, 3 figures) | `uv run python -m MODELING.moc.MOC_2D_steady_irrotational_main` |
| RDE MoC model (injection, triple point, domain) | `uv run python -m MODELING.moc.rde.MoC_Main_Model` |
| MoC helper scripts | `uv run python -m MODELING.moc.minimum_length_nozzle` (or `MODELING.moc.expansion`, `MODELING.moc.oblique_shock`) |

Every tool is also a function with the MATLAB inputs as keyword defaults. For
example, `HADES_size_main(plot=False)` returns a results object, and
`HeatFluxCalc(path, equations="cubic", depths_m=(x1, x2))` does the heat-flux
calculation. See each package's `PORTING_NOTES.md` and the module docstrings.

## Repository layout

| Path | Contents |
|---|---|
| `HADES/sizing_model/` | Engine sizing chain (`HADES_size_*.py`), CEA wrappers, cooldown/thermal models, cell-size data (`.xlsx`) |
| `HADES/feed_sizing/` | Sonic nozzle/orifice sizing V1-V5 + feed-line prototype, blowdown Cd, bundled CSV inputs |
| `HADES/data_analysis/` | Heat flux, fake ITP/plenum data, FFT, mass flow, CDR/upstream conditions, BOS backgrounds |
| `DEIMOS/` | `DEIMOS_sys_CdCalc_v1.py`, `DEIMOS_sys_FannoFlowCalc_v2.py` |
| `MODELING/moc/` | MoC nozzle/plume solver (GPL-3.0, see below); `rde/` = team RDE MoC model |
| `rdelib/` | Shared helpers: `units`, `gasdynamics`, `cea`, `matlab_compat`, `ode` (`ode45`, `pdepe`) |
| `tests/` | One folder per area, with committed Octave/MATLAB reference data |
| `Pre-Detonator/` | Pre-detonator scripts (Python, maintained separately, not part of the package or test suite) |

## Tests

```bash
uv run pytest -q          # full suite; live Octave checks skip unless set up (below)
uv run pytest tests/hades_sizing -q
```

The tests pin the ports to reference numbers from the original MATLAB, run in
Octave 11.3 and committed under `tests/*/reference` and `tests/modeling/data`.
They also include hand calculations and textbook values. CI
(`.github/workflows/tests.yml`) runs the suite on Ubuntu, macOS and Windows
with Python 3.12.

What was ported, the bugs fixed (with before/after numbers), the solver
substitutions and the validation evidence are in each area's notes:

- `HADES/sizing_model/PORTING_NOTES.md`
- `HADES/feed_sizing/PORTING_NOTES.md`
- `HADES/data_analysis/PORTING_NOTES.md`
- `DEIMOS/PORTING_NOTES.md`
- `MODELING/moc/PORTING_NOTES.md` (plus `MODELING/moc/README.md`)

## Original MATLAB and live Octave checks

Commit `2a30c95` is `main` before the rewrite and has every original file.
Extract it outside the working tree:

```bash
mkdir -p /tmp/rde_src
git archive 2a30c95 | tar -x -C /tmp/rde_src      # or: git archive 2a30c95 "HADES/SIZING MODEL" | tar -x -C /tmp/rde_src
```

To rerun the originals in GNU Octave next to the Python, point the tests at
that tree:

```bash
RDE_MATLAB_SOURCE_ROOT=/tmp/rde_src uv run pytest -q
```

`RDE_MATLAB_SOURCE_ROOT` is the root of the extracted tree, the folder that
holds `DEIMOS/`, `HADES/` and `MODELING/`. With it set and GNU Octave on
`PATH`, the live checks for DEIMOS, HADES data analysis and MoC run. The
heat-flux check also needs the Octave Forge `signal` package (`sgolayfilt`):
`octave --eval "pkg install -forge signal"`, or the OS package (e.g.
`apt install octave-signal`). If anything is missing, the check skips and the
skip reason says what to install or gives this recovery command. The HADES
sizing and feed-sizing Octave references came from a manual harness (no live
test). To redo them, see
`tests/hades_sizing/reference/octave_harness/HARNESS_NOTES.txt` and the
"Reproducing the Octave comparison" section of
`HADES/feed_sizing/PORTING_NOTES.md`. To regenerate the MoC fixtures, see
`MODELING/moc/PORTING_NOTES.md`.

## Licensing

- `MODELING/moc/` is a Python translation of *Nozzle_MOC_supersonic* by
  Xavier Dechamps and is **GPL-3.0** (`MODELING/moc/LICENSE`, modification
  notice in `MODELING/moc/README.md`).
- The rest of the repository has **no license**.
- The `rde` package (`pyproject.toml`) bundles both. The team must decide on
  licensing **before distributing** it, for example by licensing the whole
  repository under GPL-compatible terms or by splitting `MODELING/moc` out.

## Open questions for the team

These are physics and engineering decisions the ports flagged. Each port kept
the original behaviour, unless marked "fixed" in its PORTING_NOTES.md. Details
and numbers are in the notes named in brackets.

**HADES sizing model** [`HADES/sizing_model/PORTING_NOTES.md` §4-5]

| Question | Current behaviour |
|---|---|
| Thrust station | The legacy default uses the CEA **throat** column (Isp 927.8 m/s): **863.8 N**. The opt-in exit-plane thrust (`station="exit"`, `exit_thrust_from_cea`) gives **843.7 N** at sea level. It relies on CEA's 10:1 `supar` fallback, because no nozzle geometry is defined. Choose the station and define the nozzle. |
| Fill-height CJ inputs | Hand-entered (1931.5 m/s, a2 1075.6 m/s, P2 32.33). They are not the CEA results (1929.7 m/s, 1074.2 m/s, 28.927 bar). |
| Max-flux unit label | Says kg/(s·mm²), but the value is lbm/(s·mm²). |
| Annulus gap | Fixed 9.525 mm. The 1.4 × cell-width correlation (8.917 mm) is only printed. |
| Hoop-stress failure temperature | 902.9 °C comes from extrapolating the cubic fit beyond the data (max 1500 °F). |
| Other inputs to confirm | `h2Inj` density 0.91 is flagged "UNSURE". Sonic-nozzle c_f is 0.6989 (ideal for γ = 1.4: 0.6847). The `P_plenum_fanno` Sutherland viscosity mixes T and T_stag, and its Re uses the guessed Mach. The `thermal_v1` CEA property table cannot be regenerated from one condition. |

**HADES feed sizing** [`HADES/feed_sizing/PORTING_NOTES.md` issues 1-6]

| Question | Current behaviour |
|---|---|
| Feed-line prototype choking check | The legacy pressure-ratio check (`P_t,u > P_m/crit_ratio`) is **not a valid choking criterion**. It is printed as a labelled legacy diagnostic, and no replacement criterion exists. Prototype quirks are kept: the 0.85 margin only scales the Mach flow, `P_t,u` has no Cd, and R = 287. |
| Blowdown Cd method | The two estimators disagree on the committed (synthetic-looking) data: least-squares **0.640** against integrated **0.564**. `average_Cd` divides a run-average flow by an instantaneous ideal flow. Use neither until real data are run. |
| Real-gas γ in an ideal-gas choked formula | γ = cp/cv ≈ 1.585 at about 10 MPa gives choked flux about **0.24 % high** compared with real-gas isentropic flow. γ = 1.4 would be about 4 % low. |
| V1 downstream pressure | The comment says 20 bar, but the code adds 20 psi (21.38 bar). |

**HADES data analysis** [`HADES/data_analysis/PORTING_NOTES.md`, open questions 1-6]

| Question | Current behaviour |
|---|---|
| Heat flux Eq. 3/4 | **Invalid.** Their units are inconsistent (K·m³/s² against K/m), the transient sign is wrong, and they fail an exact heat-equation solution (−8288 against +7800 W/m²). There is no default: `--equations cubic` (derived replacement) or `legacy` must be chosen. Check the source paper. |
| Sensor depths | The MATLAB has 0.25 / 0.575 with no units, **probably mm**. Read as metres, the cubic model gives a mean q0 of 1.4e8 W/m². Read as mm, it gives 3.41 MW/m². Confirm the units. |
| Fake plenum data | `sin(t·f)` treats f as rad/s, `f3` is never used, and the oscillation is `abs`-rectified. |
| Fake ITP data | The comment says PSI 30-70, but the values span 30-200/218 and the axis label says kPa. |
| Mass flow from plenum | Uses static plenum pressure as stagnation pressure. |
| CDR air Reynolds number | Puts the whole 0.907 kg/s through one 2 mm hydraulic diameter (Re 3.2e7). It may need per-passage flow. |

**DEIMOS** [`DEIMOS/PORTING_NOTES.md`, behaviour notes]

| Question | Current behaviour |
|---|---|
| Methane orifice past critical ratio | P2/P01 = 0.5135 < 0.5439, so the subcritical orifice equation is out of range. It gives **32.25 bar**, where a choked formula gives **32.18 bar**. The port warns and sets `orifice_choked`. |
| Re from the Mach guess | Re uses `M_G`, not the solved Mach. Using the solved Mach changes P01 by at most 0.006 %. |
| Inputs to confirm | Viscosity uses T_stag while everything else uses static T. `Rmid` is used as a radius (the comment says diameter). `L` is assumed to be in mm. |

**MoC / RDE model** [`MODELING/moc/PORTING_NOTES.md`, "Behaviour preserved on purpose"]

| Question | Current behaviour |
|---|---|
| Triple point | **No solution for the default inputs.** The flow-angle residual is at least 0.1485 rad, so the printed β, slip angle, P and M3 are a minimum-residual point, not a solution. |
| `InjV.Cp = 315` | Equal to `InjV.R = 315`. About 1200 J/(kg·K) is expected for γ = 1.357. |
| Sharp-corner `rhod = 1e-6·yt` | Causes 17 "weak shock" patch deletions and about **4 % mass loss** along the characteristics (0.977 → 0.937 of 1-D choked). |
| `out.Vwave(invalid) = NaN` | Writes to a field that does not exist. `V4` was probably intended. |
| Far-plume complex nodes | Occur in axisymmetric cases. They are not physical, and the port warns. |

**Constants that disagree across scripts** (each port keeps its script's own value)

| Quantity | Values |
|---|---|
| psi → Pa | 6894.76 (feed V1-V3, CEA wrappers), 6894.757293 (feed V5), 6894.75729 (blowdown) |
| psi per bar | 14.5038 (`rdelib.units`, sizing), 14.504 (DEIMOS, `P_plenum_fanno` printout) |
| lb per kg | 2.205 (sizing), 2.20462 (data analysis mass flow) |
| R_H2 [J/(kg·K)] | 4157, 4124.2, 4126, 4124 (feed), CoolProp |
| H2/air stoichiometry | FAR 0.0292, O/F 34.3, AFR 34.5 |
| Air / H2 injector area [mm²] | 304.77 / 304.52 / 316.04; 25.00 / 25.13 |

## Naming convention

ALL USERS MUST FOLLOW THE NAMING CONVENTION BELOW. NOTE: THIS CONVENTION IS
SUBJECT TO CHANGE.

`PROJECT_subteam_category_description_vN.py`, for example
`HADES_size_thermal_heat_flux_v1.py`.

How it applies to the Python code:

- **Modules** follow `PROJECT_subteam_category_description_vN` where the
  original file did: `HADES/sizing_model/HADES_size_*.py` and
  `DEIMOS/DEIMOS_sys_*_v1/_v2.py`. Superseded versions were merged into one
  canonical module, and the PORTING_NOTES.md file maps say which. For example,
  `HADES_size_main_v1.m` became `HADES_size_main.py`.
- **Folders** are lowercase, importable packages without spaces
  (`sizing_model`, `feed_sizing`, `data_analysis`, `moc`). Folder names like
  `SIZING MODEL` cannot be imported.
- **Descriptive names** were used where the originals did not follow the
  convention: `HADES/feed_sizing` (`sonic_nozzle_sizing.py` for
  `SonicNozzleSizingV1-V5.m`, `blowdown_orifice_cd.py`, `choked_flow.py`) and
  `HADES/data_analysis` (`heat_flux_calc.py`, `purpl_fft.py`, ...).
  `MODELING/moc` keeps the upstream MATLAB function names, one module per
  function.
- New code should use the convention, in a lowercase package folder.

**For general bookkeeping, DO NOT NAME VERSIONS as FINAL until the code is used in CDR.**

### HADES

`HADES_subteam_category_description_versionnumber`

| Subteam | Code |
| ------------------------- | ------- |
| RDE Sizing Model | `size` |
| CAD | `cad` |
| CFD | `cfd` |
| FEA | `fea` |
| Integration | `int` |
| Testing / Instrumentation | `test` |
| Pre-Detonator | `pd` |
| Logistics / Systems | `sys` |
| Manufacturing | `mfg` |
| Media | `media` |

Examples: `HADES_size_thermal_heat_flux_v1.py`, `HADES_size_thermal_heat_flux_v2.py`,
`HADES_cfd_solver_axial_v3.py`, `HADES_fea_struct_wall_stress_v1.py`,
`HADES_pd_flow_ignition_delay_v1.py`, `HADES_test_val_hotfire_tc_v2.py`.

### DEVIL

`DEVIL_subteam_category_description_versionnumber`

| Subteam | Code |
| ------------------- | ------- |
| Sizing | `size` |
| CAD | `cad` |
| CFD | `cfd` |
| FEA | `fea` |
| Integration | `int` |
| Testing | `test` |
| Pre-Detonator | `pd` |
| Systems / Logistics | `sys` |
| Manufacturing | `mfg` |
| Media | `media` |

Examples: `DEVIL_size_thermal_heat_flux_v1.py`, `DEVIL_size_thermal_heat_flux_v2.py`,
`DEVIL_cfd_solver_axial_v3.py`, `DEVIL_fea_struct_wall_stress_v1.py`,
`DEVIL_pd_flow_ignition_delay_v1.py`, `DEVIL_test_val_hotfire_tc_v2.py`,
`DEVIL_int_perf_mass_budget_v1.py`, `DEVIL_mfg_util_tolerance_stackup_v1.py`.

### Category tags

Use consistent, appropriate `<category>` tags:

| Category | Code |
| ------------ | --------- |
| Geometry | `geom` |
| Thermal | `thermal` |
| Structural | `struct` |
| Flow | `flow` |
| Performance | `perf` |
| Optimization | `opt` |
| Solver | `solver` |
| Validation | `val` |
| Plotting | `plot` |
| Utilities | `util` |
| Run Script | `run` |
