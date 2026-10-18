"""Port of ``HADES_size_thrust.m`` (+ an opt-in exit-plane thrust option)."""

from __future__ import annotations

__all__ = ["HADES_size_thrust", "exit_thrust_from_cea"]

_LBM_TO_KG = 0.453592
_N_TO_LBF = 0.2248


def HADES_size_thrust(isp, total_mdot_lb, *, station: str = "throat", exit_pressure=None,
                      ambient_pressure=None, exit_area=None):
    """Thrust [N].

    ``station="throat"`` (default, the MATLAB formula): ``F = isp * mdot``.
    ``HADES_size_main`` passes the CEA *throat* Isp here, exactly as the
    MATLAB did, so this is the legacy sonic-exit momentum estimate.

    ``station="exit"`` (explicit opt-in): NASA rocket thrust equation

        F = mdot * Ve + (pe - p0) * Ae

    with ``isp`` = exit velocity Ve [m/s] (CEA's exit-station "Isp", which
    assumes ambient = exit pressure), ``exit_pressure`` pe [Pa],
    ``ambient_pressure`` p0 [Pa] and ``exit_area`` Ae [m^2], all required.

    ``total_mdot_lb`` is in lbm/s (converted with 0.453592 kg/lbm).
    """
    total_mdot = total_mdot_lb * _LBM_TO_KG  # lb/s -> kg/s
    if station == "throat":
        if any(v is not None for v in (exit_pressure, ambient_pressure, exit_area)):
            raise ValueError("exit_pressure/ambient_pressure/exit_area are only used with station='exit'")
        thrust = isp * total_mdot
        thrust_lbf = thrust * _N_TO_LBF
        print("Thrust: %.3f N (%.3f lbf) " % (thrust, thrust_lbf))
        return thrust
    if station == "exit":
        if exit_pressure is None or ambient_pressure is None or exit_area is None:
            raise ValueError("station='exit' requires exit_pressure [Pa], ambient_pressure [Pa] and exit_area [m^2]")
        if exit_area < 0 or exit_pressure < 0 or ambient_pressure < 0:
            raise ValueError("exit_area and pressures must be non-negative")
        momentum = total_mdot * isp
        pressure = (exit_pressure - ambient_pressure) * exit_area
        thrust = momentum + pressure
        print("Thrust: %.3f N (%.3f lbf) " % (thrust, thrust * _N_TO_LBF))
        print("  exit plane: momentum %.3f N + pressure %.3f N (pe = %.4f bar, p0 = %.4f bar, Ae = %.4e m^2)"
              % (momentum, pressure, exit_pressure / 1e5, ambient_pressure / 1e5, exit_area))
        return thrust
    raise ValueError("station must be 'throat' or 'exit'")


def exit_thrust_from_cea(ceaRock_results, total_mdot_lb, ambient_pressure):
    """Exit-plane thrust [N] from a ``HADES_size_ceaRocket`` result (opt-in).

    Uses the full-precision CEA solution (``.raw``): throat area
    ``At = mdot c* / Pc``, ``Ae = (Ae/At)_exit * At``, exit pressure and exit
    velocity from the exit station, then :func:`HADES_size_thrust` with
    ``station="exit"``. ``ambient_pressure`` in Pa.
    """
    raw = ceaRock_results.raw
    if raw is None:
        raise ValueError("ceaRock_results.raw (full CEA solution) is required")
    mdot = total_mdot_lb * _LBM_TO_KG
    pc_pa = raw.Pc_bar * 1e5
    At = mdot * raw.cstar[-1] / pc_pa
    Ae = raw.ae_at[-1] * At
    return HADES_size_thrust(raw.isp[-1], total_mdot_lb, station="exit", exit_pressure=raw.P_bar[-1] * 1e5,
                             ambient_pressure=ambient_pressure, exit_area=Ae)
