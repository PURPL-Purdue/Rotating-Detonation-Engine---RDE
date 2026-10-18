"""rdelib.units / gasdynamics / matlab_compat / ode: Octave + hand-calculation checks."""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from rdelib import gasdynamics as gd
from rdelib.matlab_compat import colon, fzero, num2str
from rdelib.ode import ode45, pdepe
from rdelib.units import convert

# Cases evaluated with the ORIGINAL HADES_size_convert.m under Octave (octave_units.json)
CONVERT_CASES = [(1, "atm", "bar"), (1, "bar", "atm"), (100, "psia", "bar"), (2, "bar", "mmHg"),
                 (10, "m/s", "ft/s"), (3, "ft/s", "m/s"), (300, "K", "C"), (25, "degC", "F"),
                 (500, "R", "K"), (77, "F", "K"), (0.25, "in", "m"), (2, "ft", "in"),
                 (1, "lbms", "kgs"), (3, "kgs", "lbms")]


def test_convert_matches_octave_except_documented_atm_fix(octave_units):
    for (val, f, t), oct_val in zip(CONVERT_CASES, octave_units["convert"]):
        py = convert(val, f, t)
        if "atm" in (f, t):
            # MATLAB table had atm = 1.01325 per bar (inverted). Fixed: 1 atm = 1.01325 bar.
            assert py == pytest.approx(1.0 / oct_val, rel=1e-12), (f, t)
        else:
            assert py == pytest.approx(oct_val, rel=1e-13), (f, t)
    assert convert(1, "atm", "bar") == pytest.approx(1.01325, rel=1e-12)
    assert convert(1, "atm", "psia") == pytest.approx(14.696, rel=1e-4)


def test_convert_errors_and_arrays():
    with pytest.raises(ValueError):
        convert(1, "bar", "m/s")
    with pytest.raises(ValueError):
        convert(1, "K", "bar")  # MATLAB: 'Target unit not recognized.'
    np.testing.assert_allclose(convert(np.array([1.0, 2.0]), "bar", "psia"), [14.5038, 29.0076])


def test_isentropic_hand_values():
    # NACA 1135 table, gamma = 1.4, M = 2
    assert gd.P0_over_P(2.0, 1.4) == pytest.approx(7.8244, rel=1e-4)
    assert gd.T0_over_T(2.0, 1.4) == pytest.approx(1.8, rel=1e-12)
    assert gd.rho0_over_rho(2.0, 1.4) == pytest.approx(4.3469, rel=1e-4)
    assert gd.area_ratio(2.0, 1.4) == pytest.approx(1.6875, rel=1e-4)
    assert gd.critical_pressure_ratio(1.4) == pytest.approx(0.528282, rel=1e-5)


def test_choked_flow_consistency():
    A = gd.choked_throat_area(1.0, 2e6, 300.0, 1.4, 287.0)
    assert gd.choked_mass_flow(A, 2e6, 300.0, 1.4, 287.0) == pytest.approx(1.0, rel=1e-12)
    # hand calc: mdot = A P0 sqrt(g/R) (2/(g+1))^((g+1)/(2(g-1))) / sqrt(T0)
    hand = 1.0 * math.sqrt(300.0) / (2e6 * math.sqrt(1.4 / 287.0) * (2 / 2.4) ** 3.0)
    assert A == pytest.approx(hand, rel=1e-12)


def test_fanno_table_values():
    # Fanno table gamma = 1.4, M = 0.5: 4fL*/D = 1.06906, p/p* = 2.1381
    assert gd.fanno_fLD(0.5, 1.4) == pytest.approx(1.06906, rel=1e-5)
    assert gd.fanno_p_over_pstar(0.5, 1.4) == pytest.approx(2.1381, rel=1e-4)
    assert gd.fanno_mach_from_fLD(1.06906, 1.4) == pytest.approx(0.5, rel=1e-5)
    # supersonic branch: M = 2 -> 4fL*/D = 0.30499
    assert gd.fanno_mach_from_fLD(gd.fanno_fLD(2.0, 1.4), 1.4, "supersonic") == pytest.approx(2.0, rel=1e-10)
    assert gd.fanno_fLD(2.0, 1.4) == pytest.approx(0.30499, rel=1e-4)


def test_friction_factors():
    Re, rr = 1e5, 1e-4
    f = gd.colebrook(Re, rr)
    assert 1 / math.sqrt(f) == pytest.approx(-2 * math.log10(rr / 3.7 + 2.51 / (Re * math.sqrt(f))), rel=1e-8)
    # Moody chart value ~0.0185 for Re = 1e5, eps/D = 1e-4
    assert f == pytest.approx(0.0185, abs=5e-4)
    # explicit Tkachenko form agrees with Colebrook to < 0.5 %
    assert gd.tkachenko_friction_factor(Re, rr) == pytest.approx(f, rel=5e-3)


def test_sutherland_and_orifice():
    assert gd.sutherland_viscosity(273.0, 84e-7, 273.0, 72.0) == pytest.approx(84e-7)
    A = gd.incompressible_orifice_area(0.1, 0.7, 1000.0, 1e5)
    assert A == pytest.approx(0.1 / (0.7 * math.sqrt(2e8)))


def test_colon_matches_matlab_semantics():
    np.testing.assert_allclose(colon(0, 0.1, 0.3), [0, 0.1, 0.2, 0.3])  # inclusive despite FP
    assert colon(0, 0.001, 0.1499).size == 150
    np.testing.assert_array_equal(colon(20, 1, 60), np.arange(20, 61))
    assert colon(1, 1, 0).size == 0
    x = colon(0, 0.1, 1.0)
    assert x.size == 11 and x[-1] == 1.0


def test_num2str_matches_matlab():
    assert num2str(25.407516934888024) == "25.4075"  # printed by Octave's num2str in the h2Inj line
    assert num2str(32.349917) == "32.3499"
    assert num2str(3.14159265) == "3.1416"
    assert num2str(0.333333333) == "0.33333"
    assert num2str(123456.789) == "123456.789"
    assert num2str(42.0) == "42"


def test_fzero_bracket_search_and_root():
    f = lambda x: (x - 1.0) * (x - 3.0) * (x + 2.0)  # noqa: E731
    # MATLAB expands around x0 = 2.6: first sign change is near 3
    assert fzero(f, 2.6) == pytest.approx(3.0, abs=1e-12)
    assert fzero(f, [0.0, 2.0]) == pytest.approx(1.0, abs=1e-12)
    assert fzero(lambda x: x**2 - 2, 1.0) == pytest.approx(math.sqrt(2), abs=1e-14)


def test_ode45_against_analytic_and_solve_ivp():
    sol = ode45(lambda t, y: -2.0 * y, [0, 3], [1.0])
    np.testing.assert_allclose(sol.y[:, 0], np.exp(-2.0 * sol.t), rtol=2e-3, atol=1e-6)
    # Refine = 4: every 4th output is a solver step; 3 interpolated points in between
    assert (sol.t.size - 1) % 4 == 0 and sol.t[0] == 0 and sol.t[-1] == 3
    # harmonic oscillator with tight tolerances
    s2 = ode45(lambda t, y: [y[1], -y[0]], [0, 10], [1.0, 0.0], rtol=1e-8, atol=1e-10)
    np.testing.assert_allclose(s2.y[:, 0], np.cos(s2.t), atol=1e-6)
    # requested output points
    s3 = ode45(lambda t, y: -y, np.linspace(0, 1, 11), [1.0], rtol=1e-6, atol=1e-9)
    np.testing.assert_allclose(s3.t, np.linspace(0, 1, 11))
    np.testing.assert_allclose(s3.y[:, 0], np.exp(-s3.t), rtol=1e-5)


def _robin_cyl_setup():
    ri, ro, k, rho, cp = 0.05, 0.0627, 15.0, 7930.0, 500.0
    hi, ho, Tg, Ta = 2000.0, 15.0, 1400.0, 25.0
    pde = lambda x, t, u, ux: (rho * cp, k * ux, 0.0)  # noqa: E731
    ic = lambda x: 25.0  # noqa: E731
    bc = lambda xl, ul, xr, ur, t: (hi * (Tg - ul), 1.0, ho * (ur - Ta), 1.0)  # noqa: E731
    return ri, ro, k, hi, ho, Tg, Ta, pde, ic, bc


def test_pdepe_steady_state_matches_analytic():
    ri, ro, k, hi, ho, Tg, Ta, pde, ic, bc = _robin_cyl_setup()
    r = np.linspace(ri, ro, 40)
    t = np.array([0.0, 2e3, 2e4])
    T = pdepe(1, pde, ic, bc, r, t)
    # analytic steady state: thermal resistances in series (per unit length)
    Rtot = 1 / (hi * 2 * math.pi * ri) + math.log(ro / ri) / (2 * math.pi * k) + 1 / (ho * 2 * math.pi * ro)
    q = (Tg - Ta) / Rtot
    Ti = Tg - q / (hi * 2 * math.pi * ri)
    Tr = Ti - q * np.log(r / ri) / (2 * math.pi * k)
    np.testing.assert_allclose(T[-1], Tr, atol=0.05)


def test_pdepe_transient_converges_with_mesh():
    ri, ro, k, hi, ho, Tg, Ta, pde, ic, bc = _robin_cyl_setup()
    t = np.linspace(0, 2, 41)
    coarse = pdepe(1, pde, ic, bc, np.linspace(ri, ro, 120), t)
    fine = pdepe(1, pde, ic, bc, np.linspace(ri, ro, 481), t)
    # inner-wall temperature history agrees to < 0.3 C (of ~1000 C rise)
    assert np.max(np.abs(coarse[:, 0] - fine[:, 0])) < 0.3


def test_pdepe_m0_against_semi_infinite_solid_analytic():
    """Slab (m = 0) heated by convection: short-time surface temperature equals the
    Carslaw & Jaeger semi-infinite-solid solution
    T_s/T_g = 1 - exp(b^2) erfc(b), b = h sqrt(alpha t)/k."""
    from scipy.special import erfc

    L, k, rc, h, Tg = 0.05, 15.0, 7930.0 * 500.0, 5000.0, 1000.0
    pde = lambda x, t, u, ux: (rc, k * ux, 0.0)  # noqa: E731
    bc = lambda xl, ul, xr, ur, t: (h * (Tg - ul), 1.0, 0.0, 1.0)  # noqa: E731
    x = np.linspace(0, L, 801)
    t = np.array([0.0, 0.1, 0.5, 1.0, 2.0])
    T = pdepe(0, pde, lambda x: 0.0, bc, x, t)
    b = h * np.sqrt(k / rc * t[1:]) / k
    analytic = Tg * (1 - np.exp(b**2) * erfc(b))
    np.testing.assert_allclose(T[1:, 0], analytic, rtol=5e-4)


@pytest.mark.parametrize("label,call,exc", [
    ("stiff, negative t (was a hang)", lambda: ode45(lambda t, y: -1e20 * y, [-1, 0], [1.0]), RuntimeError),
    ("max_step = 0 (was a hang)", lambda: ode45(lambda t, y: -y, [0, 1], [1.0], max_step=0), ValueError),
    ("NaN derivative (was silent NaN output)",
     lambda: ode45(lambda t, y: np.full_like(y, np.nan), [0, 1], [1.0]), FloatingPointError),
    ("derivative shape", lambda: ode45(lambda t, y: [1.0, 2.0], [0, 1], [1.0]), ValueError),
    ("equal tspan ends", lambda: ode45(lambda t, y: -y, [1, 1], [1.0]), ValueError),
    ("non-monotonic tspan", lambda: ode45(lambda t, y: -y, [0, 1, 0.5, 2], [1.0]), ValueError),
    ("non-finite tspan", lambda: ode45(lambda t, y: -y, [0, np.inf], [1.0]), ValueError),
    ("rtol = 0", lambda: ode45(lambda t, y: -y, [0, 1], [1.0], rtol=0), ValueError),
    ("atol < 0", lambda: ode45(lambda t, y: -y, [0, 1], [1.0], atol=-1e-6), ValueError),
    ("initial_step = 0", lambda: ode45(lambda t, y: -y, [0, 1], [1.0], initial_step=0), ValueError),
    ("NaN y0", lambda: ode45(lambda t, y: -y, [0, 1], [np.nan]), ValueError),
    ("blow-up to inf", lambda: ode45(lambda t, y: y**2, [0, 2], [1.0]), (RuntimeError, FloatingPointError)),
    ("overflowed state with finite derivative (was a silent [inf])",
     lambda: ode45(lambda t, y: np.full_like(y, 1e308), [0, 10], [1e308], max_step=1), FloatingPointError),
])
def test_ode45_failure_paths_raise(label, call, exc):
    with pytest.raises(exc):
        call()


def test_ode45_negative_and_backward_time():
    s = ode45(lambda t, y: -y, [-2, -1], [1.0])
    assert s.t[-1] == -1 and s.y[-1, 0] == pytest.approx(math.exp(-1), rel=1e-3)
    b = ode45(lambda t, y: -y, [1, 0], [1.0])
    assert b.t[-1] == 0 and b.y[-1, 0] == pytest.approx(math.e, rel=1e-3)
