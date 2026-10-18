"""Regression + validation tests for DEIMOS/DEIMOS_sys_FannoFlowCalc_v2.py.

Reference numbers come from running the original, unmodified
DEIMOS_sys_FannoFlowCalc_v2.m in GNU Octave 11.3.0. Octave lacks
``optimoptions``, so a /tmp shim mapped it to ``optimset``:
  * OCTAVE_TIGHT:   optimset('Display','off','TolX',1e-15,'TolFun',1e-15)  -> converged roots
  * OCTAVE_DEFAULT: optimset('Display','off') (fsolve TolX = TolFun = 1e-6, like MATLAB defaults)
"""

import math
import warnings

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

from DEIMOS.DEIMOS_sys_FannoFlowCalc_v2 import (  # noqa: E402
    DEIMOS_sys_FannoFlowCalc_v2,
    P01_diagnostic,
    colebrook_darcy,
    colebrook_residual,
    critical_pressure_ratio,
    fanno_fLD,
    fanno_fLD_max,
    fanno_mach_from_fLD,
    fanno_P_Pstar,
    fanno_T_Tstar,
    format_P01_residual,
    format_summary,
    main,
    orifice_mdot,
    plot_P01_diagnostic,
    solve_P01,
)

OCTAVE_COMMON = {
    "m_dotT": 0.77750997364022789,
    "m_dot": [0.58679998010583223, 0.19070999353439552],
    "A_inj": [0.00017671458676442585, 5.6548667764616275e-05],
    "A_t": 1404.3385143889795,
    "A_gap": 9.9992154495312597,
    "Rinner": 17.352887575234369,
    "Router": 27.352103024765633,
}
OCTAVE_TIGHT = {
    "iter": 3,
    "err": 7.2842344571200627e-07,
    "mu": [2.2063093228778096e-05, 1.1300468029867609e-05],
    "Re": [340577.92740586167, 218923.04944664691],
    "f": [0.032366479757164453, 0.04390301389593132],
    "M_TubeIN": [0.75336715789771191, 0.62987044452577401],
    "P01": [2231718.1405727649, 3225099.1294345274],
    "P2": [1377959.3548060455, 1656078.6139927288],
    "T_stag": [315.29104028453435, 300.56208284821508],
    "P01_residual_last": -6.4116189024687742e-09,
}
OCTAVE_DEFAULT = {
    "iter": 3,
    "mu": [2.2063092060560644e-05, 1.130046788630528e-05],
    "Re": [340578.06926223641, 218923.08282627017],
    "f": [0.032366479481505769, 0.043903013863406372],
    "M_TubeIN": [0.75336690931243655, 0.6298703613068336],
    "P01": [2231718.3043564064, 3225099.0950944368],
    "P2": [1377959.8558358378, 1656078.8454706646],
    "T_stag": [315.29101907366004, 300.56207824722134],
}

RULE = "=" * 52
OCTAVE_RESIDUAL_BLOCK = (
    f"\n{RULE}\n                P01 Solver Residual                \n{RULE}\n\n"
    "Residual                        : -0.000000006411619 \n\n"
)
OCTAVE_SUMMARY_BLOCK = (
    f"\n{RULE}\n              FLOW CALCULATION SUMMARY              \n{RULE}\n\n"
    "--- Mass Flow Rates ---\n"
    "Total Mass Flow Rate            : 0.778 kg/s\n"
    "Oxygen Mass Flow Rate           : 0.587 kg/s\n"
    "Methane Mass Flow Rate          : 0.191 kg/s\n\n"
    "--- Annular Throat Geometry ---\n"
    "Outer Radius                    : 27.35 mm\n"
    "Inner Radius                    : 17.35 mm\n"
    "Annulus Gap                     : 10.00 mm\n\n"
    "--- Fanno Output & Pressure Correlations ---\n"
    "Oxygen Inlet Mach               : 0.75 \n"
    "Methane Inlet Mach              : 0.63 \n"
    "P2 Oxygen Downstream Pressure   : 13.78 bar | 199.86 psi\n"
    "P2 Methane Downstream Pressure  : 16.56 bar | 240.20 psi\n"
    "P01 Oxygen Upstream Pressure    : 22.32 bar | 323.69 psi\n"
    "P01 Methane Upstream Pressure   : 32.25 bar | 467.77 psi\n\n"
    f"{RULE}\n\n"
)


@pytest.fixture(scope="module")
def res():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return DEIMOS_sys_FannoFlowCalc_v2()


# ---------------------------------------------------------------------------
# Full model vs the original MATLAB (run in Octave)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(OCTAVE_COMMON))
def test_sizing_matches_octave(res, name):
    np.testing.assert_allclose(getattr(res, name), OCTAVE_COMMON[name], rtol=1e-14, atol=0)


@pytest.mark.parametrize("name", [k for k in OCTAVE_TIGHT if k not in ("iter", "err", "P01_residual_last")])
def test_converged_solution_matches_octave_tight(res, name):
    # Same converged roots: brentq (here) vs Octave fsolve/fzero at 1e-15 tolerances.
    np.testing.assert_allclose(getattr(res, name), OCTAVE_TIGHT[name], rtol=1e-10, atol=0)


def test_iteration_count_and_residual_match_octave(res):
    assert res.iter == OCTAVE_TIGHT["iter"]
    assert res.err == pytest.approx(OCTAVE_TIGHT["err"], rel=1e-6)
    # Residual is ~1e-8 kg/s (P01 was solved with the previous-iteration P2); match to 1e-6 relative.
    assert res.P01_residual[-1] == pytest.approx(OCTAVE_TIGHT["P01_residual_last"], rel=1e-6)


@pytest.mark.parametrize("name", [k for k in OCTAVE_DEFAULT if k != "iter"])
def test_matches_octave_with_matlab_default_fsolve_tolerances(res, name):
    # fsolve with TolX = TolFun = 1e-6 stops ~3e-7 short of the converged root.
    np.testing.assert_allclose(getattr(res, name), OCTAVE_DEFAULT[name], rtol=2e-6, atol=0)
    assert res.iter == OCTAVE_DEFAULT["iter"]


def test_printed_output_matches_octave(res, capsys):
    assert format_P01_residual(res) == OCTAVE_RESIDUAL_BLOCK
    assert format_summary(res) == OCTAVE_SUMMARY_BLOCK
    with pytest.warns(UserWarning, match="Propellant 2: orifice pressure ratio"):
        main(diagnostic=True, show=False)
    assert capsys.readouterr().out == OCTAVE_RESIDUAL_BLOCK + OCTAVE_SUMMARY_BLOCK
    plt.close("all")


# ---------------------------------------------------------------------------
# Hand calculations / internal consistency of the converged state
# ---------------------------------------------------------------------------


def test_hand_calc_sizing(res):
    OF = 4 / 1.3
    assert res.OF_ratio == pytest.approx(3.0769230769, rel=1e-10)
    assert res.m_dotT == pytest.approx(1780 / (233.37 * 9.81), rel=1e-15)  # 0.77751 kg/s
    assert res.m_dot[0] == pytest.approx(0.77750997 * OF / (1 + OF), rel=1e-7)
    assert res.m_dot.sum() == pytest.approx(res.m_dotT, rel=1e-15)
    assert res.A_inj[0] == pytest.approx(36 * math.pi / 4 * 2.5e-3**2, rel=1e-14)
    assert res.A_inj[1] == pytest.approx(2 * 36 * math.pi / 4 * 1.0e-3**2, rel=1e-14)
    assert res.A_t == pytest.approx(0.77750997364 * 1806.2 / 1e6 * 1e6, rel=1e-10)  # mm^2
    assert res.A_gap == pytest.approx(1404.3385143889795 / (2 * math.pi * 22.3524953), rel=1e-14)
    assert res.Router - res.Rinner == pytest.approx(res.A_gap, rel=1e-12)


def test_converged_state_is_self_consistent(res):
    L = np.array([9.48944, 9.56818])
    D = np.array([2.5, 1.0])
    gam = np.array([1.4, 1.31])
    for i in range(2):
        # Fanno: f L / D equals f L*/D at the inlet Mach (tube choked at the exit)
        assert fanno_fLD(res.M_TubeIN[i], gam[i]) == pytest.approx(res.f[i] * L[i] / D[i], rel=1e-12)
        # Colebrook satisfied by the friction factor of the last iteration
        assert colebrook_residual(res.f[i], res.Re[i], 0.015, D[i]) == pytest.approx(0, abs=1e-12)
    np.testing.assert_allclose(res.P2, 1e6 * fanno_P_Pstar(res.M_TubeIN, gam), rtol=1e-15)
    np.testing.assert_allclose(res.T_stag, 283.15 * (1 + (gam - 1) / 2 * res.M_TubeIN**2), rtol=1e-15)
    assert res.err <= 1e-6
    # Orifice validity flags (P2/P01 = 0.617 oxygen, 0.513 methane vs 0.528 / 0.544 critical)
    assert list(res.orifice_choked) == [False, True]


# ---------------------------------------------------------------------------
# Flow-relation helpers vs textbook tables / independent formulas
# ---------------------------------------------------------------------------

# Fanno-flow table, gamma = 1.4 (NACA Report 1135 / Anderson, Modern Compressible Flow, Table A.4)
FANNO_TABLE_G14 = [
    # M,   T/T*,   P/P*,   4fL*/D
    (0.2, 1.1905, 5.4554, 14.533),
    (0.5, 1.1429, 2.1381, 1.0691),
    (0.8, 1.0638, 1.2893, 0.07229),
    (1.0, 1.0000, 1.0000, 0.0),
    (2.0, 0.66667, 0.40825, 0.30500),
    (3.0, 0.42857, 0.21822, 0.52216),
]


@pytest.mark.parametrize("M,TT,PP,fLD", FANNO_TABLE_G14)
def test_fanno_relations_match_textbook_table(M, TT, PP, fLD):
    assert fanno_T_Tstar(M, 1.4) == pytest.approx(TT, rel=1e-4)
    assert fanno_P_Pstar(M, 1.4) == pytest.approx(PP, rel=1e-4)
    assert fanno_fLD(M, 1.4) == pytest.approx(fLD, rel=1e-4, abs=1e-12)


def test_fanno_fLD_max_textbook():
    assert fanno_fLD_max(1.4) == pytest.approx(0.82153, rel=1e-4)


@pytest.mark.parametrize("M", [0.05, 0.3, 0.5, 0.75, 0.95, 0.999])
def test_fanno_inversion_subsonic(M):
    assert fanno_mach_from_fLD(fanno_fLD(M, 1.31), 1.31, 0.5) == pytest.approx(M, rel=1e-10)


@pytest.mark.parametrize("M", [1.01, 1.5, 2.0, 3.0, 8.0])
def test_fanno_inversion_supersonic(M):
    assert fanno_mach_from_fLD(fanno_fLD(M, 1.4), 1.4, 2.0) == pytest.approx(M, rel=1e-9)


def test_fanno_textbook_inversion():
    assert fanno_mach_from_fLD(1.0691, 1.4, 0.5) == pytest.approx(0.5, abs=1e-4)
    assert fanno_mach_from_fLD(0.3050, 1.4, 1.5) == pytest.approx(2.0, abs=1e-3)


def test_fanno_bad_guess_and_no_supersonic_root():
    with pytest.raises(ValueError, match="Don't guess 1"):
        fanno_mach_from_fLD(0.2, 1.4, 1.0)
    with pytest.raises(ValueError, match="No supersonic Fanno solution"):
        fanno_mach_from_fLD(0.9, 1.4, 2.0)


@pytest.mark.parametrize("Re", [5e3, 1e5, 2.2e5, 3.4e5, 1e7])
@pytest.mark.parametrize("eps_D", [0.0, 1e-4, 0.006, 0.015])
def test_colebrook_vs_fixed_point_and_explicit_fits(Re, eps_D):
    f = colebrook_darcy(Re, eps_D, 1.0)
    # Independent algorithm: fixed-point iteration x <- -2 log10(eps/3.7D + 2.51 x / Re)
    x = 8.0
    for _ in range(200):
        x = -2 * math.log10(eps_D / 3.7 + 2.51 * x / Re)
    assert f == pytest.approx(1 / x**2, rel=1e-12)
    # Explicit approximations of Colebrook: Haaland (~2%), Swamee-Jain (~1-2.5% for eps/D <= 1e-2)
    haaland = (-1.8 * math.log10((eps_D / 3.7) ** 1.11 + 6.9 / Re)) ** -2
    assert f == pytest.approx(haaland, rel=0.025)
    if eps_D <= 1e-2:  # Swamee-Jain validity range
        swamee_jain = 0.25 / math.log10(eps_D / 3.7 + 5.74 / Re**0.9) ** 2
        assert f == pytest.approx(swamee_jain, rel=0.03)


def test_colebrook_moody_value():
    # Smooth pipe, Re = 1e5: Colebrook / Moody chart f = 0.0180
    assert colebrook_darcy(1e5, 0.0, 1.0) == pytest.approx(0.0180, rel=5e-3)


@pytest.mark.parametrize("gam", [1.4, 1.31])
@pytest.mark.parametrize("ratio", [0.95, 0.8, 0.6, 0.53])
def test_orifice_mdot_vs_mach_form(gam, ratio):
    """Independent form: mdot = Cd A P0 sqrt(g/(R T0)) M (1 + (g-1)/2 M^2)^(-(g+1)/(2(g-1)))."""
    Cd, A, R, T0, P0 = 0.6, 1.7e-4, 259.8, 283.15, 2.2e6
    M = math.sqrt(2 / (gam - 1) * (ratio ** (-(gam - 1) / gam) - 1))
    ref = Cd * A * P0 * math.sqrt(gam / (R * T0)) * M * (1 + (gam - 1) / 2 * M**2) ** (-(gam + 1) / (2 * (gam - 1)))
    assert orifice_mdot(P0, ratio * P0, Cd, A, gam, R, T0) == pytest.approx(ref, rel=1e-12)


def test_critical_ratio_and_choked_flux_textbook():
    assert critical_pressure_ratio(1.4) == pytest.approx(0.528282, rel=1e-6)
    # Mass-flow parameter at M = 1 for gamma = 1.4: mdot sqrt(R T0) / (A P0) = 0.684731
    flux = orifice_mdot(1.0, critical_pressure_ratio(1.4), 1.0, 1.0, 1.4, 1.0, 1.0)
    assert flux == pytest.approx(0.684731, rel=1e-6)


def test_solve_P01_roundtrip_and_window_error():
    P2, mdot = 1.38e6, 0.5868
    P01 = solve_P01(mdot, P2, 0.6, 1.7671458676442585e-4, 1.4, 259.8, 283.15)
    assert orifice_mdot(P01, P2, 0.6, 1.7671458676442585e-4, 1.4, 259.8, 283.15) == pytest.approx(mdot, rel=1e-12)
    with pytest.raises(ValueError, match="Solution Window"):
        solve_P01(1e-9, P2, 0.6, 1.7e-4, 1.4, 259.8, 283.15)  # root below 1.01*P2


def test_diagnostic_curves_and_plot(res):
    x, y, ym = P01_diagnostic(res)
    assert x.shape == y.shape == ym.shape == (2, 3000)
    np.testing.assert_allclose(x[:, 0], res.P2 * (1 + 1e-6), rtol=1e-15)
    np.testing.assert_allclose(x[:, -1], res.P2 * (1 + 10**0.25), rtol=1e-15)
    assert np.all(np.diff(y, axis=1) > 0)  # MFR(P0) monotone -> unique P01 root
    fig = plot_P01_diagnostic(res)
    assert len(fig.axes[0].lines) == 6
    plt.close(fig)


# Classification verified by running modified copies of the original .m in Octave 11.3.0:
# a scalar broadcasts for these five (element-wise use only) ...
SCALAR_OK_INPUTS = ["ope", "T_stagG", "mu_ref", "T_ref", "S_ref"]
# ... and fails with "X(2): out of bound 1" for these (indexed X(i) in the MATLAB).
TWO_ENTRY_INPUTS = ["inj_OD", "M_G", "P_G", "f_G", "L", "D", "gamma", "R"]


@pytest.mark.parametrize("name", TWO_ENTRY_INPUTS)
@pytest.mark.parametrize("value", [[0.5], 0.5, [], [0.5, 0.5, 0.5]])
def test_rejects_malformed_indexed_inputs(name, value):
    # Never return a half-solved result (e.g. M_G=[0.5] used to skip methane silently).
    # 3 entries: Octave errors for all but L / f_G, which it silently truncates; the port
    # rejects those too (documented divergence, PORTING_NOTES.md).
    with pytest.raises(ValueError, match=f"{name} must have exactly 2 entries"):
        DEIMOS_sys_FannoFlowCalc_v2(**{name: value})


@pytest.mark.parametrize("name", SCALAR_OK_INPUTS)
@pytest.mark.parametrize("value", [[], [0.5, 0.5, 0.5]])
def test_rejects_malformed_broadcastable_inputs(name, value):
    # Octave: 3 entries -> "nonconformant arguments" for every one of these.
    with pytest.raises(ValueError, match=f"{name} must have exactly 2 entries .* or be a scalar"):
        DEIMOS_sys_FannoFlowCalc_v2(**{name: value})


SOLVED = ["mu", "Re", "f", "M_TubeIN", "P01", "P2", "T_stag"]
# Octave 11.3.0 (tight-tolerance shim) runs of the original v2 .m with one input made scalar.
OCTAVE_SCALAR_VARIANTS = {
    "ope=1": (
        {"ope": 1},
        {"P01": [2231718.1405727649, 9845730.3207862247], "A_inj": [0.00017671458676442585, 2.8274333882308137e-05]},
    ),
    "mu_ref=2.07e-5": (
        {"mu_ref": 2.07e-5},
        {
            "iter": 4,
            "mu": [2.2063095569355176e-05, 2.1263767098467527e-05],
            "Re": [340577.64319052146, 116444.97917177621],
            "f": [0.032366479964357528, 0.044088970465782459],
            "M_TubeIN": [0.75336715729158277, 0.62935892174492869],
            "P01": [2231717.8124260376, 3224889.7529265615],
            "P2": [1377959.3560277135, 1657502.5720615743],
            "T_stag": [315.29104023281553, 300.53381335081383],
        },
    ),
    "S_ref=127": (
        {"S_ref": 127},
        {
            "iter": 3,
            "mu": [2.2063093228778096e-05, 1.1288379692749311e-05],
            "Re": [340577.92740586167, 219157.25588609857],
            "f": [0.032366479757164453, 0.043902787012432803],
            "M_TubeIN": [0.75336715789771191, 0.62987106975456575],
            "P01": [2231718.1405727649, 3225099.3880893691],
            "P2": [1377959.3548060455, 1656076.8748872988],
            "T_stag": [315.29104028453435, 300.56211741577107],
        },
    ),
}


@pytest.mark.parametrize("label", sorted(OCTAVE_SCALAR_VARIANTS))
def test_scalar_broadcast_matches_octave(label):
    kwargs, ref = OCTAVE_SCALAR_VARIANTS[label]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = DEIMOS_sys_FannoFlowCalc_v2(**kwargs)
    for name, val in ref.items():
        if name == "iter":
            assert r.iter == val
        else:
            np.testing.assert_allclose(getattr(r, name), val, rtol=1e-10, err_msg=f"{label}: {name}")


@pytest.mark.parametrize("kwargs", [{"T_stagG": 300}, {"T_ref": 291}, {"T_stagG": [300]}, {"T_ref": 291.0}])
def test_scalar_equal_to_default_pair_reproduces_default(res, kwargs):
    # Octave: scalar T_stagG = 300 / T_ref = 291 reproduce the default run exactly.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = DEIMOS_sys_FannoFlowCalc_v2(**kwargs)
    assert r.iter == res.iter
    for name in SOLVED:
        np.testing.assert_array_equal(getattr(r, name), getattr(res, name), err_msg=name)


@pytest.mark.parametrize("kwargs", [{"max_iter": 0}, {"max_iter": -3}, {"max_iter": float("nan")}])
def test_rejects_max_iter_with_no_iteration(kwargs):
    # MATLAB/Octave: loop body never runs, then the output fails on undefined P01.
    with pytest.raises(ValueError, match="max_iter must be > 0"):
        DEIMOS_sys_FannoFlowCalc_v2(**kwargs)


def test_fractional_max_iter_matches_octave():
    # Octave: max_iter = 0.5 runs exactly one iteration (iter < 0.5 is false afterwards);
    # iter == max_iter is never true, so no non-convergence warning. P01 uses P2 = P_G = 3e6 Pa.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        r = DEIMOS_sys_FannoFlowCalc_v2(max_iter=0.5)
    assert r.iter == 1
    np.testing.assert_allclose(r.P01, [3378932.9051010655, 3801528.1574528813], rtol=1e-12)
    np.testing.assert_allclose(r.P2, [1377136.2400215755, 1655346.0472583785], rtol=1e-10)
    np.testing.assert_allclose(r.M_TubeIN, [0.75377575973268018, 0.63013391827275589], rtol=1e-10)
    assert not any("did not fully converge" in str(w.message) for w in caught)
    # max_iter = 1.5 -> two iterations in Octave
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r2 = DEIMOS_sys_FannoFlowCalc_v2(max_iter=1.5)
    assert r2.iter == 2
    np.testing.assert_allclose(r2.P01, [2231449.0873086764, 3225208.3977900683], rtol=1e-10)


@pytest.mark.parametrize("tol", [1, 1.0, 5.0, float("nan")])
def test_rejects_tol_with_no_iteration(tol):
    with pytest.raises(ValueError, match="tol must be < 1"):
        DEIMOS_sys_FannoFlowCalc_v2(tol=tol)


def test_default_tolerance_octave_differs_only_in_residual_line(res):
    # Octave with MATLAB-default fsolve tolerances printed residual -0.000000006411511; the summary
    # block is identical, only the residual's last digits differ (fsolve stops ~3e-7 short of the root).
    assert format_summary(res) == OCTAVE_SUMMARY_BLOCK
    assert res.P01_residual[-1] == pytest.approx(-6.411510600212722e-09, rel=5e-5)
    assert f"{res.P01_residual[-1]:18.15f}" != "-0.000000006411511"


def test_non_convergence_warning():
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        DEIMOS_sys_FannoFlowCalc_v2(max_iter=1)
    assert any("Iteration did not fully converge" in str(w.message) for w in caught)
