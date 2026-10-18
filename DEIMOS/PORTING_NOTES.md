# DEIMOS MATLAB -> Python porting notes

## File mapping

| Original | Python | Notes |
|---|---|---|
| `DEIMOS_sys_CdCalc_v1.m` | `DEIMOS_sys_CdCalc_v1.py` | `DEIMOS_sys_CdCalc_v1(...)` -> `CdCalcResult`; `main()` prints the MATLAB output |
| `DEIMOS_sys_FannoFlowCalc_v2.m` | `DEIMOS_sys_FannoFlowCalc_v2.py` | canonical Fanno module; `DEIMOS_sys_FannoFlowCalc_v2(...)` -> `FannoFlowResult`; `main(diagnostic=True, show=True)` prints + plots |
| `DEIMOS_sys_FannoFlowCalc_v1.m` | superseded by v2 | no unique logic (see below) |
| `test` | not ported | 1-byte placeholder (a single newline) |

The MATLAB files were scripts. Every user input is now a keyword argument with
the MATLAB default, in the order the MATLAB declared them.

## Validation

The original `.m` files were run unmodified in GNU Octave 11.3.0. Octave has no
`optimoptions`, so a temporary shim mapped it to `optimset`.

- CdCalc: every intermediate and output agrees to rel 1e-12, and the console
  output is identical.
- FannoFlowCalc v2, with Octave fsolve/fzero at 1e-15 tolerances: all solved
  quantities agree to rel 1e-10 (observed ~1e-15). Console output is identical,
  including the residual line.
- FannoFlowCalc v2, with fsolve at its default TolX/TolFun = 1e-6 (as in
  MATLAB): solved quantities agree to rel 2e-6. The FLOW CALCULATION SUMMARY
  block is identical. The residual line is not: Octave prints
  `-0.000000006411511`, Python prints `-0.000000006411619`, because fsolve
  stops about 3e-7 short of the converged root. Full stdout equality holds only
  against the tight-tolerance Octave run.
- Fanno relations match the gamma = 1.4 textbook table (NACA 1135 / Anderson
  Table A.4). Colebrook matches an independent fixed-point solution, the
  Haaland and Swamee-Jain fits, and the Moody chart. The orifice equation
  matches the Mach-number form.

The tests are in `tests/deimos/`. `test_deimos_octave_live.py` reruns the
originals in Octave and skips itself when Octave or the `.m` files are absent.

## Solver substitutions (same roots)

- **Colebrook** (`fsolve` from `f_G`): bracketed `brentq` in x = 1/sqrt(f). The
  Colebrook root is unique, so the start point `f_G` is no longer needed. It is
  still accepted and length-checked.
- **Fanno inlet Mach** (`fsolve` from `M_G`): `brentq` on the branch that `M_G`
  selects, (0, 1) or (1, inf). `M_G == 1` raises, as the MATLAB note warns.
  Asking for a supersonic root with f L/D above the limit raises.
- **P01** (`fzero` on [1.01 P2, 1.01 P2 + 1e8]): `brentq` on the same window,
  with the same error message when there is no sign change.

## MATLAB bugs

v2 and CdCalc have no coding bugs that change the numbers. v1 had two real
bugs, and v2 had already fixed both. Results below are from Octave runs of the
originals at default inputs.

1. **Reynolds number in v1 used D in mm**, missing the `/1000`, so Re was 1000x
   too large: v1 Re = 3.40e8 / 2.18e8, v2 Re = 3.41e5 / 2.19e5.
2. **v1 orifice equation was missing the 2/(gamma-1) factor**:
   `sqrt(gamma/(R T) ...)` instead of `sqrt(2 gamma/(R T (gamma-1)) ...)`.
   - P01 in v1: 63.60 bar (O2), 173.16 bar (CH4).
   - P01 in v2: 22.32 bar (O2), 32.25 bar (CH4).

Other v1 differences:
- L = [9.5 9.5] mm instead of [9.48944 9.56818].
- The P01 window was [1.01 P2, 1e7 * 1.01 P2], and a missing root gave NaN
  instead of an error.
- The Fanno expression was written twice; the second definition is
  algebraically identical to v2's.

## Behaviour notes and open questions (numbers unchanged)

- **The methane orifice is past its critical ratio.** At default inputs,
  P2/P01 = 0.5135 for CH4, below its critical ratio of 0.5439, so the
  subcritical orifice equation is out of range. O2 is fine: 0.617 against
  0.528.
  - A choked-orifice formula would give CH4 P01 = 32.18 bar instead of
    32.25 bar (-0.2%).
  - The port keeps the MATLAB formula. It sets `orifice_choked` in the result
    and emits a `UserWarning`. This warning is the only output the MATLAB did
    not produce.
- **Re uses the Mach guess `M_G`, not the solved `M_TubeIN`.** This is
  documented in the MATLAB input comment, so the port keeps it.
  - Using the solved Mach changes P2 by about 0.04% and P01 by at most 0.006%:
    P01 = 22.3156 / 32.2515 bar against 22.3172 / 32.2510 bar.
  - The result does depend weakly on the guess: M_G = 0.3 gives
    P01(O2) = 22.3203 bar.
- **Viscosity and temperature.** Viscosity is evaluated at T_stag, while the
  Re density and both orifice/Fanno temperatures use the static
  T = 283.15 K. This is kept as is.
- **`Rmid` is a radius.** The MATLAB comment says "Average annulus diameter",
  but the code uses it as a radius (gap = A_t / (2 pi Rmid)). The port keeps
  the code's meaning and corrects the comment.
- **Units of `L`.** `L` has no unit in the MATLAB. It is only used as L/D with D
  in mm, so it is taken to be mm.
- **Residual printout.** The MATLAB residual covered only the last propellant
  (methane). The port prints the same line, but `P01_residual` holds both
  values.
- **Unused inputs.** `Imp_angle` was declared but unused in the MATLAB, and the
  port keeps it that way. `f_G` was used in the MATLAB as the start point of
  the Colebrook `fsolve`. It became unnecessary only because the port uses a
  bracketed solver.
- **Input checks.** These follow the MATLAB. Each case was checked by running
  a modified temporary copy of the original `.m` in Octave.
  - **Scalars accepted for `ope`, `T_stagG`, `mu_ref`, `T_ref`, `S_ref`.**
    These appear only in element-wise expressions, so MATLAB broadcasts a
    scalar to both propellants, and the port does the same. Scalar
    `T_stagG = 300` or `T_ref = 291` reproduces the default result exactly.
  - **Exactly two entries for `inj_OD`, `M_G`, `P_G`, `f_G`, `L`, `D`, `gamma`,
    `R`.** These are indexed `X(i)`, so a scalar or 1-entry array fails in
    Octave (`X(2): out of bound 1`). The port raises `ValueError`. Before this
    check, `M_G=[0.5]` silently skipped methane.
  - **More than two entries are rejected for every input.** Octave also errors
    (`nonconformant arguments`) for all inputs except `L` and `f_G`, which it
    silently truncates to the first two entries. Rejecting `L`/`f_G` is a
    deliberate divergence.
  - **`tol >= 1` or `max_iter <= 0` (or NaN) raises `ValueError`.** No
    iteration runs, and the MATLAB then fails with `'P01' undefined`.
  - **A fractional `max_iter > 0` is valid, as in MATLAB.** `0.5` runs one
    iteration, giving P01 = [3378932.9051010655, 3801528.1574528813] Pa,
    identical to Octave. `1.5` runs two.
- **CdCalc.** The incompressible mdot used to build Re_D1 leaves out the
  velocity-of-approach factor 1/sqrt(1 - beta^4). That is 1.0008 for
  beta = 0.2, so the effect is negligible. Kept as in the MATLAB.

## Candidates for `rdelib`

These helpers in `DEIMOS_sys_FannoFlowCalc_v2.py` match the functions in
`rdelib.gasdynamics` to round-off:
`sutherland_viscosity`, `colebrook_darcy`, `fanno_fLD`, `fanno_P_Pstar`,
`fanno_mach_from_fLD`, `orifice_mdot`, `critical_pressure_ratio`.

Not yet in `rdelib`: `fanno_T_Tstar` and `fanno_fLD_max`.
