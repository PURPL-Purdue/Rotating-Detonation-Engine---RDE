"""HADES RDE sizing model (Python port of ``HADES/SIZING MODEL``).

Each MATLAB function ``HADES_size_foo`` lives in module
``HADES.sizing_model.HADES_size_foo`` and exposes a function with the same
name, the same positional arguments and the same outputs (a tuple for
multiple outputs, a dataclass for MATLAB structs, keyword arguments for
MATLAB name/value pairs). Scripts expose a function plus ``main()``.

Run the full sizing chain with::

    uv run python -m HADES.sizing_model.HADES_size_main

See ``PORTING_NOTES.md`` for the MATLAB-to-Python mapping, the validation
evidence and the MATLAB bugs that were fixed.
"""
