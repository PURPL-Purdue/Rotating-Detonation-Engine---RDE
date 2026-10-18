"""DEIMOS systems analysis tools (Python port of the original MATLAB scripts).

Modules
-------
DEIMOS_sys_CdCalc_v1
    Compressible subcritical orifice discharge-coefficient / mass-flow calculator.
DEIMOS_sys_FannoFlowCalc_v2
    Injector feed-line Fanno-flow calculator ("the big FFC"). Canonical version;
    supersedes DEIMOS_sys_FannoFlowCalc_v1 (see PORTING_NOTES.md).

Run either module as a script to reproduce the MATLAB console output, e.g.
``uv run python -m DEIMOS.DEIMOS_sys_FannoFlowCalc_v2``.
"""
