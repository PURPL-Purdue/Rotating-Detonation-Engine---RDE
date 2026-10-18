# MODELING.moc: method-of-characteristics nozzle solver and RDE MoC model

> **Modification notice (GPL-3.0 §5(a)).** These files were translated from
> MATLAB to Python on 2026-10-18 by the PURPL RDE team. They derive from
> *Nozzle_MOC_supersonic* by Xavier Dechamps, including the team's own
> modifications in `MODELING_size_Nozzle_MOC_supersonic-master/src2`.
>
> Summary of modifications relative to the upstream MATLAB:
>
> * Translated function by function into Python, with MATLAB structs mapped
>   to dataclasses and figures drawn with matplotlib.
> * Kept the team's post-processing change that colours the characteristic
>   lines by region.
> * Fixed the 1-D area ratio for planar nozzles in `get_Laval_theory`.
> * Fixed `minimum_length_nozzle` to use `tan(angle)` as the segment slope.
> * Corrected swapped axis labels in `oblique_shock`.
> * Added a warning when complex-valued plume nodes appear.
> * Added a 100 000-column guard on the open-ended marching loops.
>
> Details and before/after numbers are in [`PORTING_NOTES.md`](PORTING_NOTES.md).

**Python translation of GPL-3.0 code.** The nozzle solver in this package is a
line-by-line Python port of [*Nozzle_MOC_supersonic*](https://github.com/xavierdechamps/Nozzle_MOC_supersonic)
by Xavier Dechamps (MATLAB/Octave, GPL-3.0), including the modifications made by
the PURPL RDE team in `MODELING/MODELING_size_Nozzle_MOC_supersonic-master`.
It solves the steady, 2-D (planar or axisymmetric), irrotational supersonic flow
in a de Laval nozzle of known shape and in its plume, after Zucrow & Hoffman,
*Gas Dynamics, Vol. II: Multidimensional Flow*, Wiley, 1977. The upstream
licence is kept in [`LICENSE`](LICENSE) (GPL-3.0). This derivative work is
distributed under the same licence.

`rde/` is the Python port of the team's own RDE MoC model (`RDE_Toolbox_2`).

## Running

```bash
uv run python -m MODELING.moc.MOC_2D_steady_irrotational_main   # nozzle + plume, 3 figures
uv run python -m MODELING.moc.rde.MoC_Main_Model                 # RDE model (injection, triple point, domain)
uv run python -m MODELING.moc.minimum_length_nozzle
uv run python -m MODELING.moc.expansion
uv run python -m MODELING.moc.oblique_shock
uv run pytest tests/modeling
```

As a library:

```python
from MODELING.moc.MOC_2D_steady_irrotational_main import default_inputs, run
geom, params, plots = default_inputs()          # team inputs from the MATLAB driver
sol = run(geom, params, plots=None)             # plots=None -> no figures
Mach, p, rho, T, a = sol.thermo()
```

`sol.X, sol.Y, sol.U, sol.V` are the characteristic-net nodes with the MATLAB
layout: column `I` is one right-running (C-) characteristic, row 1 is the
wall or jet boundary, and row `LENG_INDI[I-1]` is the axis. Arrays are
zero-padded like MATLAB's growing arrays. The patch index arrays
(`patch_i_tri`, `patch_j_tri`, `patch_i_quad`, `patch_j_quad`) keep MATLAB
**1-based** values.

## Conventions

* One module per MATLAB function, with the same name, argument order and
  outputs (a tuple when there are several outputs). MATLAB structs map to
  dataclasses with the same field names (`structs.Geom`, `structs.Params`,
  `structs.Plots`, `structs.Chocked`, `rde/structs.py`). Scripts get a
  `main()`.
* `_matlab.py` reproduces the MATLAB semantics the solver depends on:
  * the `a:d:b` colon operator;
  * degree trigonometry;
  * MATLAB's bracketed `fzero` algorithm;
  * `A\b` for 2x2 systems;
  * NaN-ignoring `max`;
  * growing zero-padded arrays;
  * complex-number propagation, where relational operators compare real parts as in MATLAB.
* See [`PORTING_NOTES.md`](PORTING_NOTES.md) for bugs fixed, deviations and validation.

## Original file -> Python module

| Original (MATLAB) | Python |
|---|---|
| `MOC_2D_steady_irrotational_main.m` (team) | `MOC_2D_steady_irrotational_main.py` (`default_inputs`, `run`, `main`) |
| `src2/MOC_2D_steady_irrotational_IVLINE.m` | `MOC_2D_steady_irrotational_IVLINE.py` |
| `src2/MOC_2D_steady_irrotational_extend_initial_line.m` | `MOC_2D_steady_irrotational_extend_initial_line.py` |
| `src2/MOC_2D_steady_irrotational_extend_initial_expansion.m` | `MOC_2D_steady_irrotational_extend_initial_expansion.py` |
| `src2/MOC_2D_steady_irrotational_extend_expansion.m` | `MOC_2D_steady_irrotational_extend_expansion.py` |
| `src2/MOC_2D_steady_irrotational_extend_plume.m` | `MOC_2D_steady_irrotational_extend_plume.py` |
| `src2/MOC_2D_steady_irrotational_internal_point.m` | `MOC_2D_steady_irrotational_internal_point.py` |
| `src2/MOC_2D_steady_irrotational_wall.m` | `MOC_2D_steady_irrotational_wall.py` |
| `src2/MOC_2D_steady_irrotational_wall_inverse.m` | `MOC_2D_steady_irrotational_wall_inverse.py` |
| `src2/MOC_2D_steady_irrotational_free_pressure.m` | `MOC_2D_steady_irrotational_free_pressure.py` |
| `src2/MOC_2D_steady_irrotational_get_geometry.m` | `MOC_2D_steady_irrotational_get_geometry.py` |
| `src2/MOC_2D_steady_irrotational_get_thermo.m`, `..._get_thermo2.m` | `MOC_2D_steady_irrotational_get_thermo.py`, `..._get_thermo2.py` |
| `src2/MOC_2D_steady_irrotational_postprocess.m` (team colour fix) | `MOC_2D_steady_irrotational_postprocess.py` |
| `src2/get_speed_sound.m`, `get_prandtl_meyer_function.m`, `get_Mach_from_nu.m`, `get_Laval_theory.m`, `FinalWaveExpansion.m` | same names `.py` |
| `src2/minimum_length_nozzle.m`, `expansion.m`, `oblique_shock.m` (scripts) | same names `.py` (function + `main()`) |
| `RDE_Toolbox_2/MoC_Main_Model.m` | `rde/MoC_Main_Model.py` |
| `RDE_Toolbox_2/src/MoC_Injection_Velocity.m`, `MoC_Resolve_Triple_Point.m`, `MoC_Calculate_Domain_Size.m` | `rde/` same names |
| `RDE_Toolbox_2/post_processing/Post_Process_Main.m`, `Setup_Graph_Domain.m` | `rde/post_processing/` same names |
