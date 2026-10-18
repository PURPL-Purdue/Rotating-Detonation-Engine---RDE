# Porting notes: MODELING (MoC nozzle + RDE MoC model)

## Canonical sources

* **Solver**: `MODELING_size_Nozzle_MOC_supersonic-master/src2/`. Compared with
  the unmodified upstream copy (`MODELING_size_MoC_unchanged_backup/src/`),
  `src2` differs only in two places:
  * Line breaks in `internal_point.m`. There is no numerical change.
  * The team's patch-colour fix in `postprocess.m`. Characteristic lines are
    now coloured by region, and patches past `patches_xlim` are handled with
    NaN/zero padding.
* The four `MOC_2D_*` files in `RDE_Toolbox_2/src/` are byte-identical to the
  team's `src2/` copies. Compared with upstream `src/` they are identical only
  after newline/whitespace normalisation: `extend_expansion`, `get_geometry`
  and `wall` use different line endings. They are not separate versions.
  Separately, `src2/IVLINE` differs from upstream in whitespace only.
* **Driver**: the team's `MOC_2D_steady_irrotational_main.m` (planar,
  `rhod = 1e-6 yt`, `xplume = 10 yt`, `NIexpansion = 5`,
  `patches_data = 0`). The upstream inputs are axisymmetric, with
  `rhou = 2 yt`, `rhod = 0.5 yt`, `xplume = 20 yt`, `NIexpansion = 10` and
  `patches_data = 1`. They are covered by the regression tests and the
  README example.

## Bugs fixed (behaviour differs from MATLAB)

1. **`get_Laval_theory` uses the wrong area ratio for planar nozzles.** MATLAB
   always used `area = Y.^2/Y(1)^2`, which is axisymmetric. The team's driver
   is planar (`geom.delta = 0`), so the "1D" comparison curves in figures 2
   and 3 were wrong. The port now uses `A/A* = Y/Y(1)` when `delta = 0`; the
   axisymmetric expression is unchanged. Example: team case, last wall node
   before the lip (x = 9.39 yt). The MoC wall solution there is M = 2.502,
   p/P0 = 0.0541.

   |  | Before (MATLAB) | After (port) |
   |---|---|---|
   | 1-D M | 3.4262 | 2.5214 |
   | 1-D p/P0 | 0.00948 | 0.0522 |

   Only the plotted reference curve changes. The MoC solution does not.
   Test: `test_laval_theory_planar_fix_vs_original`.
2. **`minimum_length_nozzle.m` uses angles as slopes.** The straight
   characteristic and wall segments used the mean angle in radians
   (`angle/rad2deg`) as their slope instead of `tan(angle)`. The exit area
   ratio is therefore wrong. For M = 2.4 and gamma = 1.4 (isentropic
   A/A* = 2.4031):

   | Rays | Before (MATLAB) | After (port) |
   |---|---|---|
   | 7 | 2.7601 (+14.9 %) | 2.4045 (+0.06 %) |
   | 60 | 2.7872 | 2.4048 |

   The error does not shrink as rays are added, so it is not a
   discretisation effect. `minimum_length_nozzle(use_tangent_slopes=False)`
   reproduces the MATLAB numbers exactly; Octave validation uses that mode.
3. **`MoC_Injection_Velocity.m` does not parse.** Its last `fprintf` argument
   is `out.);`, so `MoC_Main_Model.m` could not run in MATLAB. The header
   names a "Mach #" column, and `InjV.gamma1` (the unreacted gas gamma) is
   read but never used. The port therefore prints
   `M4 = V4/sqrt(gamma1*R*T4)` in the 5th column (0.201 for the default
   inputs) and stores it as `out.M4`. The "infill angle phi" header column
   has no defined value in the original and is still not printed. The other
   outputs are unchanged and match Octave run on a copy with that one
   argument replaced.
4. **`oblique_shock.m` axis labels were swapped.** x is the deflection angle
   theta and y is the shock angle beta. This is cosmetic.
5. **The team driver's `addpath('./src/')` points to a folder that does not
   exist in the team copy**, which uses `src2/`. Python imports from the
   package, so this does not apply.

## Behaviour preserved on purpose (report to the team)

* **Complex nodes in the far plume.** In the axisymmetric cases the
  free-pressure point eventually has no real solution: `sqrt` of a negative
  number in `MOC_2D_steady_irrotational_free_pressure`. MATLAB silently
  carries complex x, y, u, v from there on. The upstream default case is
  affected only in the discarded column past `xe + xplume`; with `te = 0` or
  `PRatio = 0.8`, the last 1-2 columns are affected; with `NI = 21`, the last
  39 columns are. The port reproduces this with MATLAB semantics: values with
  a zero imaginary part are narrowed back to real, `<` and `>` compare real
  parts, `max` compares moduli, and LAPACK's complex pivoting is used. It
  also emits a `RuntimeWarning` naming the first complex column. Those nodes
  are not physical. GNU Octave orders complex numbers by modulus, so Octave
  and MATLAB results differ in that region. `_matlab.complex_ordering("octave")`
  exists only so that tests can reproduce the Octave reference there.
* **Team geometry with `PRatio = 0.8` (weak lip shock).** The original stops
  with "The maximum of iterations for the predictor-corrector algorithm has
  been reached". The port raises `MocConvergenceError`, as in MATLAB.
* **Team geometry gives many "weak shock" patch deletions.** The sharp corner
  (`rhod = 1e-6 yt`) produces 12 + 5 "weak shock" messages in the nozzle.
  Mass flow along the characteristics drops from 0.977 to 0.937 of the 1-D
  choked value. The upstream axisymmetric case conserves mass to within
  0.6 %.
* **`MoC_Resolve_Triple_Point` has no solution for the default inputs.** On
  the pressure-matched curve the flow-angle residual is at least 0.1485 rad
  (8.5 deg) and never changes sign. MATLAB's/Octave's `fsolve` and SciPy's
  `fsolve` all report non-convergence: MATLAB prints
  `warning('Triple point solver did not converge.')` and the port raises a
  `RuntimeWarning`. All three stop at the point of minimum residual
  (beta = 70.34 deg, M3 = 2.228), and Python and Octave agree there to
  4e-6 relative. The printed beta, slip angle, P and M3 are therefore not a
  triple-point solution.
* **`InjV.Cp = 315` equals `InjV.R = 315` in `MoC_Main_Model.m`.** For
  gamma = 1.357 one would expect Cp ≈ 1200 J/kg-K. The input is kept.
* **`out.Vwave(invalid) = NaN` on a field that does not exist.** This
  produces `[0, NaN]` for the default inputs. MATLAB semantics are kept,
  although `out.V4(invalid) = NaN` was probably intended.
* **Upstream quirks kept.**
  * The predictor-corrector convergence test is `abs(max([dx, dy])) < eps`,
    a signed maximum, not the maximum absolute value.
  * The region colouring in the post-processing advances one region per patch.
  * The `intersection` flag persists between marching steps. If the
    internal-point loop is empty on the first step, the port raises
    `NameError`, as MATLAB errors on the undefined variable.
* **Safety guard.** The open-ended `while 1` marching loops in
  `extend_expansion` and `extend_plume` stop after 100 000 columns with a
  `RuntimeError`. MATLAB would hang.
* **Memory.** `_matlab.GrowArray` mimics MATLAB's zero-padded growth with
  amortised doubling per axis: only the axis whose capacity is exceeded
  grows. One row grown to 8 193 columns allocates 16 x 16 384 doubles
  (2 MiB).

## Numerical-method mapping

| MATLAB | Python |
|---|---|
| `fzero(f,[a b])` | `_matlab.fzero`, a port of MATLAB's bracketed Brent-Dekker algorithm and its stopping rule |
| `fminbnd` | `scipy.optimize.fminbound` (same Brent golden-section/parabolic method), `xtol = 1e-4` (MATLAB default `TolX`) |
| `fsolve` (trust-region-dogleg, tol 1e-8) | `scipy.optimize.fsolve` (MINPACK hybrd, Powell dogleg), `xtol = 1e-8` |
| `A\b` (2x2) | `_matlab.mldivide2` (LAPACK getrf/getrs order; singular gives Inf/NaN, no exception) |
| `patch(...,'FaceColor','interp')` | `tripcolor(..., shading='gouraud')` with quads split into two triangles |
| `patch(...,'EdgeColor','flat','FaceColor','none')` | `PolyCollection(facecolors='none', edgecolors=<region colour>)` |
| `pbaspect([xMax yMax 1])` | `Axes.set_box_aspect(yMax/xMax)` |

## Validation

GNU Octave 11.3.0 runs the original, unmodified `.m` files through
`tests/modeling/octave/generate_references.m`. Input overrides are applied
only to the input section of a temporary copy of the driver. The results are
stored in `tests/modeling/data/*.mat`.

### Nozzle and plume

There are eight cases:

* team default;
* team with `te = 5`;
* team with `PRatio = 0.4` (strong shock: early stop);
* team with `PRatio = 0.8` (original errors);
* upstream default;
* upstream with `te = 0`;
* upstream with `PRatio = 0.8` (weak shock);
* upstream with `NI = 21`.

Each case has 1 574 to 5 038 nodes. For every node (x, y, u, v), the maximum
relative difference is at most 3.7e-13. Mach differs by at most 1.2e-12 and
p/P0 by at most 1e-13. `indI`, `LENG_INDI`, all patch index arrays and
`LOC_PATCHES` are identical. These figures hold over the full arrays under
Octave complex ordering; under the default MATLAB ordering they hold up to
the first complex column.

Test tolerances:

| Quantity | Tolerance |
|---|---|
| Nodes | 1e-10 relative |
| Mach | 1e-9 absolute |
| p/P0 | 1e-10 relative |
| 1-D area-Mach curve | 1e-9 relative (ill-conditioned near M = 1; fzero variants differ by 6e-12) |

The console trace (shock and expansion messages, lip Mach 2.5359 -> 2.9057,
deflection 11.3206 deg) is identical. The figures for the upstream inputs
reproduce the upstream `Images/Pressure.jpg` and `Images/Mach_number.jpg`.

### Scripts and RDE functions

| Item | Agreement with Octave |
|---|---|
| `minimum_length_nozzle` (original mode) | ≤ 1e-13 |
| `expansion` | delta to 2e-13 |
| `oblique_shock` | weak/strong beta to 6e-13 deg; fminbnd beta_max to 1e-6 deg |
| `MoC_Injection_Velocity` | exact |
| `MoC_Calculate_Domain_Size` | exact |
| `MoC_Resolve_Triple_Point` | 4e-6 (non-converged system, see above) |

### Theory (`tests/modeling/test_moc_theory.py`)

* Prandtl-Meyer function against tables, and inversion to 1e-12.
* 1-D area-Mach relation, planar and axisymmetric.
* Jet-boundary pressure equals ambient to 1e-12.
* Exit-lip Mach within 2 % of the 1-D value for the exit area ratio
  (observed 0.8 % planar, 1.0 % axisymmetric).
* Mass conservation along characteristics: within 1 % in the upstream case
  (observed 0.6 %), within 7 % in the team case.
* Fixed minimum-length-nozzle area ratio within 0.5 % of isentropic for
  M = 1.5 to 3 (observed ≤ 0.21 %).
* Theta-beta-M residuals, and theta_max for M = 5 is 41.12 deg.
* Injection-velocity root satisfies continuity, momentum and energy to 1e-12.

### Precision notes

The degree trigonometry is MATLAB-accurate (exact argument reduction).
Octave's `sind` reduces with `mod(x-180,360)-180`, and that is one source of
the remaining ulp-level differences.

## Regenerating the Octave references (after the MATLAB originals are deleted)

The `.mat` fixtures in `tests/modeling/data/` are committed, so the tests do
not need Octave or the MATLAB sources. To regenerate the fixtures, or to run
the live Octave test (`test_live_octave_team_case`), recover the originals
from git history. Commit `2a30c95` was `main` HEAD before the deletion.

```bash
mkdir -p /tmp/moc_src
git archive 2a30c95 MODELING | tar -x -C /tmp/moc_src
export RDE_MATLAB_SOURCE_ROOT=/tmp/moc_src/MODELING     # folder holding MODELING_size_*
(cd "$RDE_MATLAB_SOURCE_ROOT" && shasum -a 256 -c checksums.txt)   # optional: save the list below as checksums.txt
octave --no-gui --norc --eval "addpath([pwd '/tests/modeling/octave']); generate_references(pwd, 'tests/modeling/data')"
uv run pytest tests/modeling -q          # the live Octave test now runs instead of skipping
```

`generate_references.m` and the live test read `RDE_MATLAB_SOURCE_ROOT` and
fall back to `<repo>/MODELING`. If the sources are missing, the generator
stops with these instructions and the live test's skip reason repeats them.
Verified on 2026-10-18: regenerating from the `git archive` tree reproduced
all 12 fixture files exactly.

SHA-256 of the key original files at `2a30c95`, with paths relative to
`$RDE_MATLAB_SOURCE_ROOT`:

```text
84625f1f9a60cbb574a2d597dc0ea8a3b3761b6c7d817f9f1419f69abf0dce88  MODELING_size_Nozzle_MOC_supersonic-master/MOC_2D_steady_irrotational_main.m
b0ec7c52aec1d28340e20518a97f556534834895499eb4c60e709cf878ed822c  MODELING_size_Nozzle_MOC_supersonic-master/src2/FinalWaveExpansion.m
0d6e96f50c7e3103037c1c816324ad85054b7f5c444045d0b16216f63ac72b48  MODELING_size_Nozzle_MOC_supersonic-master/src2/MOC_2D_steady_irrotational_IVLINE.m
d57c25ae56838607e0569f9ae68b5a0724eb35286b98cf218c55b52a24a0562e  MODELING_size_Nozzle_MOC_supersonic-master/src2/MOC_2D_steady_irrotational_extend_expansion.m
c611c7e45518d78e01a2874992fc4cddd6cd1aa078271d037baf8dc671c0a770  MODELING_size_Nozzle_MOC_supersonic-master/src2/MOC_2D_steady_irrotational_extend_initial_expansion.m
bdcb70ad10e8b679eb466fc65497cd207f6354b1f6c0ab94b0bba7c7685cbb69  MODELING_size_Nozzle_MOC_supersonic-master/src2/MOC_2D_steady_irrotational_extend_initial_line.m
4c6a14ce17dfb41e4d453a9ad695cad497c4df03f0fcfa6273ee1ffcc67fc6c4  MODELING_size_Nozzle_MOC_supersonic-master/src2/MOC_2D_steady_irrotational_extend_plume.m
4cc77b72d090dcc85fb5c1bed473987910d21c361b9847ce75e5338177c5cf9b  MODELING_size_Nozzle_MOC_supersonic-master/src2/MOC_2D_steady_irrotational_free_pressure.m
22990bdf666d40eee85881f3dcd5e8a18b4ea064d532984dc13fe93dc18933a9  MODELING_size_Nozzle_MOC_supersonic-master/src2/MOC_2D_steady_irrotational_get_geometry.m
cde9fa1c53612b6425457c5513b4108facc7702ddc49b7807ff000eed29f5906  MODELING_size_Nozzle_MOC_supersonic-master/src2/MOC_2D_steady_irrotational_get_thermo.m
7eef68cd11024b56106c4037a48b7b815ac063d3c3c62211b8108e17a23de0e4  MODELING_size_Nozzle_MOC_supersonic-master/src2/MOC_2D_steady_irrotational_get_thermo2.m
cb225dcd3e258cebe0f37abe9069511e8115700c4cc9dc85ad62f55d0484b885  MODELING_size_Nozzle_MOC_supersonic-master/src2/MOC_2D_steady_irrotational_internal_point.m
30316a9962292284688ff279d822b514853b2b5cf79c87e4d3fb124f11bc1c34  MODELING_size_Nozzle_MOC_supersonic-master/src2/MOC_2D_steady_irrotational_postprocess.m
e36089f7dc059932608618dfdd35654cb547b4f1a735ca7a6ce0d4e77edca09f  MODELING_size_Nozzle_MOC_supersonic-master/src2/MOC_2D_steady_irrotational_solve_free_pressure.m
4bb08bdd42abc03e8eaf5a692967b68c6bff4d1c2dc784d6925b940dd52c8ca6  MODELING_size_Nozzle_MOC_supersonic-master/src2/MOC_2D_steady_irrotational_wall.m
745ccef86f5c14c59d57ae22ff21635ccf82f4d1576227f7008886c69610d298  MODELING_size_Nozzle_MOC_supersonic-master/src2/MOC_2D_steady_irrotational_wall_inverse.m
094fa9e5eed60cb7b3e10a5e30543fa01b878c39b70322b3e551cf9dcb8fb1d4  MODELING_size_Nozzle_MOC_supersonic-master/src2/expansion.m
b80c729fd8dbed64cce2823267620c84ee868b12907190db0d1f48a3811cff71  MODELING_size_Nozzle_MOC_supersonic-master/src2/get_Laval_theory.m
573f261dd4cce76804530fd6afdf9033fd68dd4c09b3a3d50c1a6340cde27d61  MODELING_size_Nozzle_MOC_supersonic-master/src2/get_Mach_from_nu.m
e12e28760c460800809f2cf9dab0f6f8e4b7049081ede856e4c7e4109f4d2266  MODELING_size_Nozzle_MOC_supersonic-master/src2/get_prandtl_meyer_function.m
e1519e2c579734e516d4c3beec7a5f69e1efc6d4433134fbcb993d09fa71c9c2  MODELING_size_Nozzle_MOC_supersonic-master/src2/get_speed_sound.m
3362f5e195b00d010120fc1c6f423764820baa4c2eb8317208185b9cac70bced  MODELING_size_Nozzle_MOC_supersonic-master/src2/minimum_length_nozzle.m
23c1c70805aebd449699df6354b39fa2b773b7b840521d7e5817539feb58f4c7  MODELING_size_Nozzle_MOC_supersonic-master/src2/oblique_shock.m
bd7f6ab233e46c11a0a2bb6b26f9ffd130b52341f5e70b3e323f6047d989b816  MODELING_size_Nozzle_MOC_supersonic-master/RDE_Toolbox_2/MoC_Main_Model.m
1f3b4cd4ebb08dc9dfecc747f953a23e90db7eb0f06b1e0adff03eb3094fd84b  MODELING_size_Nozzle_MOC_supersonic-master/RDE_Toolbox_2/src/MoC_Calculate_Domain_Size.m
c80f6dd9169af22ce8a79a8de866e6bf9d2649bc3cbf8057ecb5cacee917cdfc  MODELING_size_Nozzle_MOC_supersonic-master/RDE_Toolbox_2/src/MoC_Injection_Velocity.m
5e968b4edfc0a9e630536eb5a254ddd8e01b48e9d9d50a14eb3770afec3b5fc1  MODELING_size_Nozzle_MOC_supersonic-master/RDE_Toolbox_2/src/MoC_Resolve_Triple_Point.m
213c301fa505d2ed61455d2cd6efacb9e9e5ed5ec4a197049d1c8f143ee6a14d  MODELING_size_Nozzle_MOC_supersonic-master/RDE_Toolbox_2/post_processing/Post_Process_Main.m
69146b84c3d55842e454ea506990333822dd097261890d006fece4f3af547aab  MODELING_size_Nozzle_MOC_supersonic-master/RDE_Toolbox_2/post_processing/Setup_Graph_Domain.m
7770d6e449015d081a2234a0f27eefb36ec9ffedee2fa673d80ff7db74e4e12d  MODELING_size_MoC_unchanged_backup/MOC_2D_steady_irrotational_main.m
```

Aggregate for the 21 upstream solver files:
`(cd MODELING_size_MoC_unchanged_backup/src && shasum -a 256 *.m | shasum -a 256)`
= `4166a1e07be8881cb96c4819fb4d082dc0e7fb705b3a41383ebe3ea291ee6b6c`.
