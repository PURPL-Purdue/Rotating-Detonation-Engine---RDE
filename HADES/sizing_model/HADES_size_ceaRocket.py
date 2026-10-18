"""Port of ``HADES_size_ceaRocket.m`` (NASA CEA rocket-performance wrapper).

Runs the same ``rocket equilibrium`` problem the MATLAB wrote to its CEA2
deck, in-process through :func:`rdelib.cea.rocket` (cross-platform NASA
``cea`` package), with the same unit handling and ``%f`` deck formatting.
Reactants carry no temperature in the deck, so (as in CEA2) the chamber
enthalpy is zero.

STATION AMBIGUITY (see PORTING_NOTES.md): the MATLAB parser took the *first*
number on each PERFORMANCE PARAMETERS line, i.e. the **throat** column. The
legacy fields ``AeAt``, ``cf``, ``isp`` and ``ivac`` keep that MATLAB meaning
(throat; ``AeAt`` = 1) so existing callers and ``HADES_size_main`` are
unchanged. Both stations are also exposed under unambiguous names:
``*_throat`` (identical to the legacy fields) and ``*_exit`` (the column at
the requested ``areaRatio`` / ``Pe``), plus ``pe_exit_bar`` and ``pc_bar``.
Exit-based thrust (with the pressure term) is an explicit opt-in in
:func:`HADES.sizing_model.HADES_size_thrust.HADES_size_thrust`.

Numbers are rounded to CEA2's printed precision by default
(``cea2_rounding=True``); the full solution is in ``.raw``.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from rdelib import cea as _cea
from rdelib.units import CEA2_PSIA_TO_BAR, convert_pressure

from ._nvargs import deck_float, parse_name_value

__all__ = ["HADES_size_ceaRocket", "CeaRocketResults", "cea2_input_deck"]

_NAN = math.nan


@dataclass
class CeaRocketResults:
    """MATLAB ``results`` struct fields (legacy = THROAT values) + both stations.

    Legacy fields (MATLAB meaning, i.e. the CEA **throat** column):
    ``AeAt`` (= 1), ``cstar``, ``cf``, ``isp``, ``ivac``, ``Pc``, ``Pe``,
    ``PcUnits``, ``PeUnits``.
    """

    AeAt: float = _NAN  # legacy: throat Ae/At (= 1)
    cstar: float = _NAN  # [m/s] (same at every station)
    cf: float = _NAN  # legacy: throat thrust coefficient
    isp: float = _NAN  # legacy: throat Isp [m/s] (= throat sonic velocity)
    ivac: float = _NAN  # legacy: throat vacuum Isp [m/s]
    Pc: float = _NAN  # as given, in PcUnits
    Pe: float | None = None  # as given, in PeUnits
    PcUnits: str = "psia"
    PeUnits: str = "psia"
    # unambiguous station-tagged values
    AeAt_throat: float = _NAN
    cf_throat: float = _NAN
    isp_throat: float = _NAN  # [m/s]
    ivac_throat: float = _NAN  # [m/s]
    AeAt_exit: float = _NAN  # at the requested areaRatio / Pe
    cf_exit: float = _NAN
    isp_exit: float = _NAN  # [m/s] exit velocity (CEA Isp: ambient = exit pressure)
    ivac_exit: float = _NAN  # [m/s] = isp_exit + pe*Ae/mdot
    pe_exit_bar: float = _NAN  # exit static pressure [bar]
    pc_bar: float = _NAN  # chamber pressure used by CEA [bar]
    raw: _cea.RocketResult | None = None  # full-precision CEA solution


_DEFAULTS = dict(ox="", fuel="", Pc=None, PcUnits="psia", Pe=None, PeUnits="psia",
                 areaRatio=None, of=None, pctFuel=None, phi=None, r=None, ceaExe=None,
                 cea2_rounding=True)


def _prepare(args, kwargs, func="HADES_size_ceaRocket"):
    o = parse_name_value(args, kwargs, _DEFAULTS, func)
    if not o["ox"] or not o["fuel"]:
        raise ValueError("Both 'ox' and 'fuel' must be specified.")
    if o["Pc"] is None:
        raise ValueError("Chamber pressure 'Pc' must be specified.")
    if o["Pe"] is None and o["areaRatio"] is None:
        o["areaRatio"] = 10
    if o["Pe"] is not None and o["areaRatio"] is not None:
        raise ValueError("Specify either 'Pe' or 'areaRatio', not both.")
    mix = [k for k in ("of", "pctFuel", "phi", "r") if o[k] is not None]
    if not mix:
        raise ValueError("Must specify one mixture type: of, pctFuel, phi, or r.")
    if len(mix) > 1:
        raise ValueError("Specify only one mixture type.")
    o["_mix_type"], o["_mix_val"] = mix[0], o[mix[0]]
    o["_Pc_psia"] = convert_pressure(o["Pc"], o["PcUnits"], "psia")
    o["_Pe_psia"] = convert_pressure(o["Pe"], o["PeUnits"], "psia") if o["Pe"] is not None else None
    return o


def cea2_input_deck(*args, **kwargs) -> str:
    """The ``cea_run.inp`` text the MATLAB wrote for these inputs."""
    o = _prepare(args, kwargs, "cea2_input_deck")
    mt, val = o["_mix_type"], o["_mix_val"]
    if o["_Pe_psia"] is not None:
        exit_line = "    pi/p=%f\n" % (o["_Pc_psia"] / o["_Pe_psia"])
    else:
        exit_line = "    supar=%f\n" % o["areaRatio"]
    mix_line = {
        "of": "    o/f=%f\n" % val,
        "pctFuel": "    o/f=%f\n" % ((100 - val) / val),
        "phi": "    phi=%f\n" % val,
        "r": "    r,eq.ratio=%f\n" % val,
    }[mt]
    return (
        "problem\n    rocket equilibrium\n"
        + "    p,psia=%f\n" % o["_Pc_psia"]
        + exit_line
        + mix_line
        + "reac\n"
        + "    oxid %s wt=100\n" % o["ox"]
        + "    fuel %s wt=100\n" % o["fuel"]
        + "output short\nonly\nend\n"
    )


def HADES_size_ceaRocket(*args, **kwargs) -> CeaRocketResults:
    """Equilibrium rocket performance (infinite-area combustor).

    Name/value options (MATLAB names): ``ox``, ``fuel``, ``Pc`` + ``PcUnits``,
    optional exit condition ``Pe`` + ``PeUnits`` **or** ``areaRatio``
    (default ``areaRatio=10`` if neither is given), and exactly one of
    ``of``, ``pctFuel``, ``phi``, ``r``. ``ceaExe`` is accepted and ignored;
    ``cea2_rounding`` (default True) rounds outputs to CEA2's printed precision.
    """
    o = _prepare(args, kwargs)
    mix_type, val = o["_mix_type"], o["_mix_val"]
    Pc_psia, Pe_psia = o["_Pc_psia"], o["_Pe_psia"]

    Pc_bar = deck_float(Pc_psia) * CEA2_PSIA_TO_BAR  # "p,psia=%f"
    exit_kw = (
        {"pi_p": deck_float(Pc_psia / Pe_psia)}  # "pi/p=%f"
        if Pe_psia is not None
        else {"supar": deck_float(o["areaRatio"])}  # "supar=%f"
    )
    mix_kw = {
        "of": {"of": deck_float(val)} if mix_type == "of" else None,
        "pctFuel": {"of": deck_float((100 - val) / val)} if mix_type == "pctFuel" else None,
        "phi": {"phi": deck_float(val)} if mix_type == "phi" else None,
        "r": {"r": deck_float(val)} if mix_type == "r" else None,
    }[mix_type]

    rr = _cea.rocket(o["fuel"], o["ox"], Pc_bar, hc=0.0, **exit_kw, **mix_kw)

    rnd = (lambda x, k: _cea.cea2_round(float(x), k)) if o["cea2_rounding"] else (lambda x, k: float(x))
    ex, th = -1, 1  # exit (last) and throat stations
    throat = dict(AeAt=rnd(rr.ae_at[th], "var"), cf=rnd(rr.cf[th], 4),
                  isp=rnd(rr.isp[th], 1), ivac=rnd(rr.ivac[th], 1))
    return CeaRocketResults(
        **throat,  # legacy fields = MATLAB parser result (throat column)
        cstar=rnd(rr.cstar[th], 1),
        Pc=o["Pc"],
        Pe=o["Pe"],
        PcUnits=o["PcUnits"],
        PeUnits=o["PeUnits"],
        AeAt_throat=throat["AeAt"],
        cf_throat=throat["cf"],
        isp_throat=throat["isp"],
        ivac_throat=throat["ivac"],
        AeAt_exit=rnd(rr.ae_at[ex], "var"),
        cf_exit=rnd(rr.cf[ex], 4),
        isp_exit=rnd(rr.isp[ex], 1),
        ivac_exit=rnd(rr.ivac[ex], 1),
        pe_exit_bar=rnd(rr.P_bar[ex], "var"),
        pc_bar=float(Pc_bar),
        raw=rr,
    )
