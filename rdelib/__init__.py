"""Shared helpers (unit conversions, gas dynamics, thermochemistry) used across RDE projects.

Modules (import them explicitly; nothing heavy is imported here):

``rdelib.units``
    ``convert`` (port of HADES_size_convert) and the conversion constants
    used throughout the MATLAB sources.
``rdelib.gasdynamics``
    Perfect-gas isentropic, choked/sonic-orifice, Fanno, friction-factor
    (Colebrook, Tkachenko) and Sutherland relations, SI units.
``rdelib.cea``
    Cross-platform NASA CEA (pip package ``cea``): Chapman-Jouguet
    detonation, rocket performance and hp/tp equilibrium. Replaces the
    Windows-only FCEA2.exe workflow.
``rdelib.matlab_compat``
    MATLAB semantics needed by faithful ports: ``colon`` (a:step:b),
    ``num2str``, ``fzero`` (MATLAB's bracket search).
``rdelib.ode``
    ``ode45`` (MATLAB's algorithm, incl. Refine=4 output) and ``pdepe``
    (Skeel-Berzins discretisation, scalar PDE).
"""
