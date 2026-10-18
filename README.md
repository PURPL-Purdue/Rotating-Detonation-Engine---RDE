# PURPL RDE

PURPL's rotating detonation engine code, all in Python now: HADES sizing, feed sizing and test-data analysis, DEIMOS calcs, and the MoC nozzle/RDE model. It was ported from our old MATLAB and checked against the originals run in Octave. The MATLAB still lives in git history at `2a30c95`.

## Setup

Install [uv](https://docs.astral.sh/uv/), then:

```bash
uv sync
```

No more `FCEA2.exe`: CEA now comes from NASA's `cea` package, which works on Mac, Linux and Windows.

## Running stuff

Run these from the repo root. Outputs land in whatever folder you're in. Set `MPLBACKEND=Agg` if you don't want plot windows.

```bash
uv run python -m HADES.sizing_model.HADES_size_main          # full sizing chain
uv run python -m HADES.feed_sizing.sonic_nozzle_sizing       # feed sizing (v5; also v1-v4, feedline)
uv run python -m HADES.feed_sizing.blowdown_orifice_cd       # blowdown Cd
uv run python -m HADES.data_analysis.heat_flux_calc --equations cubic --depths-m 0.25e-3 0.575e-3
uv run python -m HADES.data_analysis.purpl_fft               # ITP FFT
uv run python -m DEIMOS.DEIMOS_sys_FannoFlowCalc_v2
uv run python -m MODELING.moc.MOC_2D_steady_irrotational_main
uv run python -m MODELING.moc.rde.MoC_Main_Model
```

Every module runs this way (`python -m <package>.<module>`). The ones with options (feed sizing, data analysis) show them with `--help`.

## Layout

| Folder | What's in it |
|---|---|
| `HADES/sizing_model/` | engine sizing chain, CEA wrappers, cooldown/thermal |
| `HADES/feed_sizing/` | sonic nozzle/orifice sizing, blowdown Cd |
| `HADES/data_analysis/` | heat flux, FFT, fake data, mass flow, BOS backgrounds |
| `DEIMOS/` | Cd and Fanno-flow calcs |
| `MODELING/moc/` | MoC nozzle solver (GPL-3.0) + RDE model in `rde/` |
| `rdelib/` | shared helpers: units, gas dynamics, CEA, ODE |
| `Pre-Detonator/` | pre-det scripts (separate, untouched) |

Each package has a `PORTING_NOTES.md` covering what moved where and what bugs got fixed.

## Tests

```bash
uv run pytest -q
```

CI runs the tests on Linux, Mac and Windows. To also compare against the original MATLAB (needs Octave plus its `signal` package):

```bash
mkdir -p /tmp/rde_src && git archive 2a30c95 | tar -x -C /tmp/rde_src
RDE_MATLAB_SOURCE_ROOT=/tmp/rde_src uv run pytest -q
```

## Stuff the team needs to decide

The ports kept the original behavior. These came up along the way (details in each `PORTING_NOTES.md`):

- **Thrust:** the default is still the MATLAB's throat Isp → 863.8 N. The exit-plane option (`station="exit"`) gives 843.7 N, but we have no real nozzle geometry. Pick one.
- **Heat flux:** the original Eq. 3/4 are dimensionally wrong, so you have to choose `--equations cubic` or `legacy`. The sensor depths are probably mm, not m.
- **Feed-line choke check:** the prototype's "is it choked" check isn't a valid criterion. It's labeled legacy now.
- **Blowdown Cd:** the two methods disagree (0.640 vs 0.564), and the test data looks synthetic.
- **DEIMOS methane orifice:** it's past critical, so the subcritical formula is out of range (32.25 vs 32.18 bar).
- **MoC:** the triple point has no solution for the defaults. Also `Cp = 315`, which equals R (should be ~1200), and the sharp-corner geometry loses ~4% of the mass flow.
- **Constants:** psi→Pa, R_H2, stoich ratio etc. don't match across scripts.
- **License:** `MODELING/moc` is GPL-3.0 and the rest has no license. Sort this out before sharing the package outside the team.

## Naming convention

**Everyone follows this (subject to change):** `PROJECT_subteam_category_description_vN`, e.g. `HADES_size_thermal_heat_flux_v1.py`. Folders are lowercase with no spaces, so Python can import them. Don't call anything FINAL until it's used in CDR.

| Subteam | Code | | Category | Code |
|---|---|---|---|---|
| Sizing | `size` | | Geometry | `geom` |
| CAD | `cad` | | Thermal | `thermal` |
| CFD | `cfd` | | Structural | `struct` |
| FEA | `fea` | | Flow | `flow` |
| Integration | `int` | | Performance | `perf` |
| Testing / Instrumentation | `test` | | Optimization | `opt` |
| Pre-Detonator | `pd` | | Solver | `solver` |
| Systems / Logistics | `sys` | | Validation | `val` |
| Manufacturing | `mfg` | | Plotting | `plot` |
| Media | `media` | | Utilities / Run script | `util` / `run` |

Projects: `HADES_...` and `DEVIL_...`, e.g. `HADES_cfd_solver_axial_v3.py`, `DEVIL_fea_struct_wall_stress_v1.py`.
