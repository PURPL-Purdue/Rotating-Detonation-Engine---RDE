"""Regression + validation tests for DEIMOS/DEIMOS_sys_CdCalc_v1.py.

Reference numbers were produced by running the original, unmodified
DEIMOS_sys_CdCalc_v1.m in GNU Octave 11.3.0 (no shims needed) and dumping the
workspace with %.17g.
"""

import math
import warnings

import pytest
from scipy.optimize import brentq

from DEIMOS.DEIMOS_sys_CdCalc_v1 import DEIMOS_sys_CdCalc_v1, corner_tap_Cd, format_results, main

# Octave 11.3.0 run of the original .m (default inputs)
OCTAVE_REF = {
    "A1": 9.940195505498954e-05,
    "A2": 3.9760782021995816e-06,
    "beta": 0.19999999999999998,
    "D1_in": 0.44291338582677164,
    "PR_critical": 1.8929291587378541,
    "PR_actual": 1.4615384615384615,
    "P0": 1900636.4362041964,
    "rho0": 25.837074528786886,
    "iter": 3,
    "Cd_inc": 0.61352887141875034,
    "Re_D1": 65406.242081313016,
    "mdot": 0.013580928523425545,
    "f": 0.30160345445450742,
    "r": 0.68398141550745972,
    "Kn": 0.64685211007739984,
    "term": 0.47040005633369986,
    "Cd_comp": 0.68310323937343564,
    "mdot_comp": 0.012311683882188076,
}

OCTAVE_STDOUT = (
    "\n--- Results ---\nCd (incompressible) = 0.613529\nCd (compressible)   = 0.683103\n"
    "Mass flow rate      = 0.012312 kg/s\nRe_D1               = 65406\n"
)


@pytest.fixture(scope="module")
def res():
    return DEIMOS_sys_CdCalc_v1()


@pytest.mark.parametrize("name", sorted(OCTAVE_REF))
def test_matches_octave(res, name):
    # Pure closed-form arithmetic + fixed-point loop: agreement to round-off.
    assert getattr(res, name) == pytest.approx(OCTAVE_REF[name], rel=1e-12, abs=0)


def test_printed_output_matches_octave(res, capsys):
    assert format_results(res) == OCTAVE_STDOUT
    main()
    assert capsys.readouterr().out == OCTAVE_STDOUT


def test_hand_calcs(res):
    assert res.A2 == pytest.approx(math.pi / 4 * 0.00225**2, rel=1e-15)
    assert res.beta == pytest.approx(0.2, rel=1e-15)
    assert res.D1_in == pytest.approx(0.01125 / 0.0254, rel=1e-15)
    assert res.PR_critical == pytest.approx(1.2**3.5, rel=1e-15)  # 1.89293 for gamma = 1.4
    assert res.PR_actual == pytest.approx(19 / 13, rel=1e-15)
    assert res.rho0 == pytest.approx(res.P0 / (259.8 * 283.15), rel=1e-15)


def test_P0_satisfies_isentropic_continuity(res):
    """Independent check: rho V A is equal at the pipe (P1, A1) and the bore (P2, A2) for stagnation P0."""
    g = 1.4

    def mass_flux(P, P0):
        return P ** (1 / g) * math.sqrt(1 - (P / P0) ** ((g - 1) / g))

    lhs = res.A1 * mass_flux(19e5, res.P0)
    rhs = res.A2 * mass_flux(13e5, res.P0)
    assert lhs == pytest.approx(rhs, rel=1e-9)
    assert 19e5 < res.P0 < 19e5 * 1.001  # small approach velocity for beta = 0.2


def test_Cd_inc_is_fixed_point_of_correlation(res):
    """Solve Cd = C(Re(Cd)) with brentq independently of the MATLAB loop."""
    A2, rho1, dP = res.A2, 19e5 / (259.8 * 283.15), 19e5 - 13e5

    def g(Cd):
        Re = 4 * (Cd * A2 * math.sqrt(2 * rho1 * dP)) / (math.pi * 0.01125 * 2.35e-5)
        return corner_tap_Cd(0.2, 0.01125 / 0.0254, Re) - Cd

    Cd_star = brentq(g, 0.3, 1.0, xtol=1e-15)
    assert res.Cd_inc == pytest.approx(Cd_star, rel=1e-6)  # loop tolerance is 1e-6 relative


def test_incompressible_limit_recovers_Cd_inc():
    """As P2 -> P1 the compressible Cd must reduce to the incompressible Cd (f/Kn construction)."""
    r = DEIMOS_sys_CdCalc_v1(P2=19e5 * (1 - 1e-5))
    assert r.Cd_comp == pytest.approx(r.Cd_inc, rel=1e-4)


def test_mdot_comp_equals_isentropic_orifice_flow_times_Cd_comp(res):
    g, r = 1.4, res.r
    ideal = res.A2 * math.sqrt(2 * g / (g - 1) * res.P0 * res.rho0 * (r ** (2 / g) - r ** ((g + 1) / g)))
    assert res.mdot_comp == pytest.approx(res.Cd_comp * ideal, rel=1e-12)


def test_choked_warning():
    with pytest.warns(UserWarning, match="CHOKED"):
        DEIMOS_sys_CdCalc_v1(P1=30e5, P2=13e5)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        DEIMOS_sys_CdCalc_v1()  # default case is subcritical: no warning


@pytest.mark.parametrize(("maxiter", "expected_iter"), [(100.0, 3), (1.5, 1), (2.5, 2)])
def test_non_integer_maxiter_floors_like_matlab(maxiter, expected_iter):
    # Octave: `for iter = 1:maxiter` with these limits stops after 3, 1 and 2 iterations.
    result = DEIMOS_sys_CdCalc_v1(maxiter=maxiter)
    assert result.iter == expected_iter
    if expected_iter == 3:
        assert result.Cd_inc == pytest.approx(OCTAVE_REF["Cd_inc"], rel=1e-12)


@pytest.mark.parametrize("maxiter", [0, 0.5, float("nan")])
def test_rejects_maxiter_below_one(maxiter):
    with pytest.raises(ValueError):
        DEIMOS_sys_CdCalc_v1(maxiter=maxiter)
