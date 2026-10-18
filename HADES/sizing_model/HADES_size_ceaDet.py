"""Port of ``HADES_size_ceaDet.m`` (NASA CEA Chapman-Jouguet detonation wrapper).

The MATLAB wrote a CEA2 input deck, ran the Windows-only ``FCEA2.exe`` and
parsed the printed ``.out`` file. This port calls NASA's cross-platform CEA
(``cea`` package) in-process through :func:`rdelib.cea.cj_detonation`, with
the identical problem definition:

* the inputs go through the same unit conversions and the same ``%f``
  (6-decimal) formatting the deck used, then CEA2's psia -> bar factor;
* by default every returned number is rounded to the precision CEA2 printed
  in ``cea_det.out`` (which is what the MATLAB parser returned), so results
  are bit-for-bit what the MATLAB produced. Pass ``cea2_rounding=False`` for
  full precision. The full-precision solution is always in ``.raw``.

Usage (both forms work)::

    d = HADES_size_ceaDet(ox="Air", fuel="H2", phi=0.9, P0=26, P0Units="psia",
                          T0=283.15, T0Units="K")
    d = HADES_size_ceaDet('ox', 'Air', 'fuel', 'H2', 'phi', 0.9, 'P0', 26, ...)
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from rdelib import cea as _cea
from rdelib.units import CEA2_PSIA_TO_BAR, convert_pressure, convert_temperature

from ._nvargs import deck_float, parse_name_value

__all__ = ["HADES_size_ceaDet", "CeaDetResults", "cea2_input_deck"]

_NAN = math.nan


@dataclass
class CeaDetResults:
    """Fields of the MATLAB ``results`` struct (same names, same units)."""

    P_ratio: float = _NAN  # P/P1
    T_ratio: float = _NAN  # T/T1
    M_ratio: float = _NAN  # M/M1
    rho_ratio: float = _NAN  # RHO/RHO1
    detMach: float = _NAN  # detonation Mach number
    cjVel: float = _NAN  # CJ velocity [m/s]
    P_burned_bar: float = _NAN  # burned (CJ) pressure [bar]
    T_cj: float = _NAN  # burned gas temperature [K]
    Cp_eq: float = _NAN  # [kJ/(kg K)], equilibrium
    k_eq: float = _NAN  # conductivity [mW/(cm K)], equilibrium
    Pr_eq: float = _NAN
    Cp_frozen: float = _NAN
    k_frozen: float = _NAN
    Pr_frozen: float = _NAN
    Mu: float = _NAN  # viscosity [millipoise]
    Son_speed_unburned: float = _NAN  # [m/s]
    Son_speed_burned: float = _NAN  # [m/s]
    P0: float = _NAN
    T0: float = _NAN
    P0Units: str = "psia"
    T0Units: str = "K"
    raw: _cea.DetonationResult | None = None  # full-precision CEA solution


_DEFAULTS = dict(ox="", fuel="", P0=None, P0Units="psia", T0=None, T0Units="K",
                 of=None, pctFuel=None, phi=None, r=None, ceaExe=None, cea2_rounding=True)


def _prepare(args, kwargs, func="HADES_size_ceaDet"):
    o = parse_name_value(args, kwargs, _DEFAULTS, func)
    if not o["ox"] or not o["fuel"]:
        raise ValueError("Both oxidizer and fuel must be specified.")
    if o["P0"] is None:
        raise ValueError("Initial pressure P0 must be specified.")
    if o["T0"] is None:
        raise ValueError("Initial temperature T0 must be specified.")
    mix = [k for k in ("of", "pctFuel", "phi", "r") if o[k] is not None]
    if len(mix) != 1:
        raise ValueError("Specify exactly one mixture type.")
    o["_mix_type"], o["_mix_val"] = mix[0], o[mix[0]]
    o["_P0_psia"] = convert_pressure(o["P0"], o["P0Units"], "psia")
    o["_T0_K"] = convert_temperature(o["T0"], o["T0Units"], "K")
    return o


def cea2_input_deck(*args, **kwargs) -> str:
    """The ``cea_det.inp`` text the MATLAB wrote for these inputs (for reference/CEA2 users)."""
    o = _prepare(args, kwargs, "cea2_input_deck")
    mt, val = o["_mix_type"], o["_mix_val"]
    mix_line = {
        "of": "    o/f=%f\n" % val,
        "pctFuel": "    o/f=%f\n" % ((100 - val) / val),
        "phi": "    phi=%f\n" % val,
        "r": "    r,eq.ratio=%f\n" % val,
    }[mt]
    return (
        "problem\n    det\n"
        + "    p,psia=%f\n" % o["_P0_psia"]
        + "    t,k=%f\n" % o["_T0_K"]
        + mix_line
        + "reac\n"
        + "    oxid %s wt=100\n" % o["ox"]
        + "    fuel %s wt=100\n" % o["fuel"]
        + "output transport\noutput short\nonly\nend\n"
    )


def HADES_size_ceaDet(*args, **kwargs) -> CeaDetResults:
    """CJ detonation properties of an oxidizer/fuel mixture.

    Name/value options (MATLAB names): ``ox``, ``fuel`` (CEA species names,
    e.g. 'Air', 'O2', 'H2'), ``P0`` + ``P0Units`` ('psia', 'atm', 'bar',
    'mmHg'), ``T0`` + ``T0Units`` ('K', 'C', 'F', 'R'), and exactly one of
    ``of`` (O/F weight ratio), ``pctFuel`` (weight % fuel), ``phi``
    (weight-based equivalence ratio), ``r`` (valence-based equivalence ratio).
    ``ceaExe`` is accepted for compatibility and ignored; ``cea2_rounding``
    (default True) rounds outputs to CEA2's printed precision.
    """
    o = _prepare(args, kwargs)
    mix_type, val = o["_mix_type"], o["_mix_val"]

    # Same values the deck carried ("p,psia=%f", "t,k=%f", "phi=%f", ...)
    P1_bar = deck_float(o["_P0_psia"]) * CEA2_PSIA_TO_BAR
    T1_K = deck_float(o["_T0_K"])
    mix_kw = {
        "of": {"of": deck_float(val)} if mix_type == "of" else None,
        "pctFuel": {"of": deck_float((100 - val) / val)} if mix_type == "pctFuel" else None,
        "phi": {"phi": deck_float(val)} if mix_type == "phi" else None,
        "r": {"r": deck_float(val)} if mix_type == "r" else None,
    }[mix_type]

    d = _cea.cj_detonation(o["fuel"], o["ox"], P1_bar, T1_K, transport=True, **mix_kw)

    rnd = (lambda x, k: _cea.cea2_round(x, k)) if o["cea2_rounding"] else (lambda x, k: x)
    # Decimals = the formats CEA2 used in cea_det.out ('var' = magnitude dependent)
    return CeaDetResults(
        P_ratio=rnd(d.P_P1, 3),
        T_ratio=rnd(d.T_T1, 3),
        M_ratio=rnd(d.M_M1, 4),
        rho_ratio=rnd(d.rho_rho1, 4),
        detMach=rnd(d.det_mach, 4),
        cjVel=rnd(d.det_vel, 1),
        P_burned_bar=rnd(d.P_bar, "var"),
        T_cj=rnd(d.T_K, 2),
        Cp_eq=rnd(d.cp_eq, 4),
        k_eq=rnd(d.k_eq, "var"),
        Pr_eq=rnd(d.Pr_eq, 4),
        Cp_frozen=rnd(d.cp_frozen, 4),
        k_frozen=rnd(d.k_frozen, "var"),
        Pr_frozen=rnd(d.Pr_frozen, 4),
        Mu=rnd(d.visc_mP, "var"),
        Son_speed_unburned=rnd(d.son_vel1, 1),
        Son_speed_burned=rnd(d.son_vel, 1),
        P0=o["P0"],
        T0=o["T0"],
        P0Units=o["P0Units"],
        T0Units=o["T0Units"],
        raw=d,
    )
