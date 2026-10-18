"""HeatFluxCalc.m -> HADES.data_analysis.heat_flux_calc.

Legacy reference numbers come from running the ORIGINAL HeatFluxCalc.m in
GNU Octave 11.3 with signal 1.4.8 on sample_thermocouple_data.csv. The
"cubic" model is checked against exact solutions of the heat equation.
"""

import numpy as np
import pytest

from HADES.data_analysis import heat_flux_calc as hf
from HADES.data_analysis.heat_flux_calc import (
    LegacyEquationWarning,
    SensorDepthUnitsWarning,
)

from ._octave_harness import requires_originals, run_octave

IDX = [0, 5, 10, 250, 1000, 1990, 1999]
MATLAB_DEPTHS = hf.MATLAB_SENSOR_DEPTHS_M
MM_DEPTHS = (0.25e-3, 0.575e-3)

# Octave run of the original script (MATLAB derivative T_s/dt, see PORTING_NOTES.md)
OCT_T1_s = [300.45284605406954, 300.81252641680345, 301.5006982487341, 412.2796016968209,
            492.8827796443178, 499.79818621414194, 499.5071839917674]  # fmt: skip
OCT_T2_s = [299.17735399132386, 300.0247344217132, 300.6678301233223, 398.8941495892718,
            492.8618008344253, 499.6521543441854, 500.0243298383548]  # fmt: skip
OCT_T0 = [306.0033192476629, 305.9995816461415, 306.7323533258506, 428.6970839869174,
          500.42127192456826, 507.5368179076366, 506.7397313349482]  # fmt: skip
OCT_q0 = [-3821.74980797266, -4418.339004739977, -4376.081104785938, 8835.620598346934,
          -8776.082927626332, -8748.245374132957, -9545.56092904021]  # fmt: skip


def _exact(t, A=300.0, B=-20.0, C=50.0, D=0.0, alpha=hf.ALPHA, k=hf.K):
    """Exact heat-equation solution T = A + Bx + Cx^2 + Dx^3 + alpha*(2C + 6Dx)*t.

    Returns (T(x, t), dT/dt(x), T0(t), q0(t)) with q0 = -k dT/dx at x = 0.
    """
    T = lambda x, tt=t: A + B * x + C * x**2 + D * x**3 + alpha * (2 * C + 6 * D * x) * tt  # noqa: E731
    Tdot = lambda x: alpha * (2 * C + 6 * D * x) + 0 * x  # noqa: E731
    return T, Tdot, T(0.0), -k * (B + 6 * alpha * D * t)


@pytest.fixture(scope="module")
def legacy():
    with pytest.warns(LegacyEquationWarning), pytest.warns(SensorDepthUnitsWarning):
        return hf.HeatFluxCalc(equations="legacy", depths_m=MATLAB_DEPTHS, plot=False)


# ------------------------------------------------------------------ API rules
def test_equations_and_depths_are_required():
    with pytest.raises(TypeError):
        hf.HeatFluxCalc(depths_m=MM_DEPTHS, plot=False)  # no default model
    with pytest.raises(TypeError):
        hf.HeatFluxCalc(equations="cubic", plot=False)  # no default depths
    with pytest.raises(ValueError):
        hf.HeatFluxCalc(equations="paper", depths_m=MM_DEPTHS, plot=False)
    with pytest.raises(ValueError):
        hf.HeatFluxCalc(equations="cubic", depths_m=(0.5e-3, 0.2e-3), plot=False)


def test_cli_requires_explicit_choices(no_show):
    for argv in ([], ["--equations", "cubic"], ["--matlab-depths"], ["--equations", "nope", "--matlab-depths"]):
        with pytest.raises(SystemExit):
            hf.main(argv)
    hf.main(["--equations", "cubic", "--depths-m", "0.25e-3", "0.575e-3"])
    with pytest.warns(LegacyEquationWarning), pytest.warns(SensorDepthUnitsWarning):
        hf.main(["--equations", "legacy", "--matlab-depths"])


def test_warnings_only_when_warranted(recwarn):
    hf.HeatFluxCalc(equations="cubic", depths_m=MM_DEPTHS, plot=False)
    assert not [w for w in recwarn if issubclass(w.category, (LegacyEquationWarning, SensorDepthUnitsWarning))]


# ------------------------------------------------------------------ cubic model
@pytest.mark.parametrize("D", [0.0, 30.0])
def test_cubic_equations_exact_for_heat_equation_solutions(D):
    T, Tdot, T0_exact, q0_exact = _exact(0.3, D=D)
    x1, x2 = MATLAB_DEPTHS
    args = (T(x1), T(x2), Tdot(x1), Tdot(x2), x1, x2)
    assert hf.surface_temperature_cubic(*args) == pytest.approx(T0_exact, rel=1e-13)
    assert hf.surface_heat_flux_cubic(*args) == pytest.approx(q0_exact, rel=1e-12)


def test_cubic_reviewer_case_and_surface_sensor():
    # T = 300 - 20x + 50x^2 + 100*alpha*t  -> q0 = +7800 W/m^2 (k = 390)
    T, Tdot, T0_exact, q0_exact = _exact(0.7, C=50.0)
    assert q0_exact == pytest.approx(7800.0)
    for x1, x2 in (MATLAB_DEPTHS, MM_DEPTHS):
        q0 = hf.surface_heat_flux_cubic(T(x1), T(x2), Tdot(x1), Tdot(x2), x1, x2)
        assert q0 == pytest.approx(7800.0, rel=1e-10)
    # sensor 1 on the surface -> T0 == T1 exactly
    assert hf.surface_temperature_cubic(410.0, 400.0, 5.0, 7.0, 0.0, 1e-3) == 410.0


@pytest.mark.parametrize("depths", [MATLAB_DEPTHS, MM_DEPTHS])
@pytest.mark.parametrize("D", [0.0, 30.0])
def test_cubic_pipeline_through_csv_recovers_exact_flux(tmp_path, depths, D):
    # The exact solutions are linear in time, so the SG smoothing/derivative is exact too.
    t = np.arange(400) * 1e-3
    T, Tdot, T0_exact, q0_exact = _exact(t, D=D)
    x1, x2 = depths
    path = tmp_path / "tc.csv"
    np.savetxt(path, np.column_stack([t, T(x1), T(x2)]), delimiter=",", header="time_s,T1_K,T2_K", comments="")
    with pytest.warns(SensorDepthUnitsWarning) if depths == MATLAB_DEPTHS else _no_warning():
        r = hf.HeatFluxCalc(path, equations="cubic", depths_m=depths, plot=False)
    np.testing.assert_allclose(r.q0, q0_exact, rtol=1e-7, atol=1e-6)
    np.testing.assert_allclose(r.T0, T0_exact, rtol=1e-8)  # 1/alpha amplifies fp rounding of dT/dt
    if D == 0.0:
        assert r.q0.mean() == pytest.approx(7800.0, rel=1e-7)


def _no_warning():
    import contextlib

    return contextlib.nullcontext()


def test_cubic_regression_with_mm_depths():
    r = hf.HeatFluxCalc(equations="cubic", depths_m=MM_DEPTHS, plot=False)
    assert r.q0.mean() == pytest.approx(3412650.9664990315, rel=1e-9)  # ~3.4 MW/m^2
    np.testing.assert_allclose(r.T0[[0, 250, 1000, 1999]], [301.4668156, 422.77234138, 492.90738266, 499.06990442], rtol=1e-9)


# ------------------------------------------------------------------ legacy model
def test_legacy_equations_fail_exact_solution():
    T, Tdot, T0_exact, q0_exact = _exact(0.3)
    x1, x2 = MATLAB_DEPTHS
    args = (T(x1), T(x2), Tdot(x1), Tdot(x2), x1, x2)
    with pytest.warns(LegacyEquationWarning):
        q0 = hf.surface_heat_flux_legacy(*args)
    with pytest.warns(LegacyEquationWarning):
        T0 = hf.surface_temperature_legacy(*args)
    assert q0_exact == pytest.approx(7800.0)
    assert q0 == pytest.approx(-8287.50019821409, rel=1e-10)  # wrong sign and size
    assert T0 == pytest.approx(292.81583016941374, rel=1e-12)


def test_legacy_steady_state_and_coefficients_hand_calculation():
    with pytest.warns(LegacyEquationWarning):
        # no transient: q0 = -k (T2-T1)/(x2-x1) = 390*10/0.325 = 12000 W/m^2
        assert hf.surface_heat_flux_legacy(400.0, 390.0, 0.0, 0.0, *MATLAB_DEPTHS) == pytest.approx(12000.0, rel=1e-12)
        assert hf.surface_temperature_legacy(400.0, 390.0, 0.0, 0.0, *MATLAB_DEPTHS) == pytest.approx(
            (0.575 * 400 - 0.25 * 390) / 0.325, rel=1e-14
        )
        # Eq. 4 coefficients with x1=1, x2=2, alpha=k=1: (4/6, 5/6); Eq. 3: (0, 3/6)
        assert hf.surface_heat_flux_legacy(0.0, 0.0, 6.0, 0.0, 1.0, 2.0, k=1.0, alpha=1.0) == pytest.approx(-4.0)
        assert hf.surface_heat_flux_legacy(0.0, 0.0, 0.0, 6.0, 1.0, 2.0, k=1.0, alpha=1.0) == pytest.approx(-5.0)
        assert hf.surface_temperature_legacy(0.0, 0.0, 6.0, 6.0, 1.0, 2.0, alpha=1.0) == pytest.approx(3.0)


def test_legacy_pinned_to_octave_with_original_derivative(legacy):
    # Feed MATLAB's "derivative" (T_s/dt) through the legacy equations: this must
    # reproduce the original script's T0 and q0.
    np.testing.assert_allclose(legacy.T1_s[IDX], OCT_T1_s, rtol=1e-13)
    np.testing.assert_allclose(legacy.T2_s[IDX], OCT_T2_s, rtol=1e-13)
    T1_dot, T2_dot = legacy.T1_s / legacy.dt, legacy.T2_s / legacy.dt
    with pytest.warns(LegacyEquationWarning):
        T0 = hf.surface_temperature_legacy(legacy.T1_s, legacy.T2_s, T1_dot, T2_dot, *MATLAB_DEPTHS)
        q0 = hf.surface_heat_flux_legacy(legacy.T1_s, legacy.T2_s, T1_dot, T2_dot, *MATLAB_DEPTHS)
    np.testing.assert_allclose(T0[IDX], OCT_T0, rtol=1e-12)
    np.testing.assert_allclose(q0[IDX], OCT_q0, rtol=1e-10, atol=1e-8)


def test_legacy_with_fixed_derivative_regression(legacy):
    # Numbers recorded in PORTING_NOTES.md (derivative fixed, legacy Eq. 3/4)
    assert hf.DATA_FILE.exists() and legacy.t.shape == (2000,)
    assert legacy.dt == pytest.approx(1e-3, rel=1e-12)
    np.testing.assert_allclose(legacy.T1_dot[[0, 250, 1000]], [39.0869, 296.903, 12.4926], rtol=1e-5)
    np.testing.assert_allclose(legacy.T0[[0, 250, 1000, 1999]], [301.436522, 422.581660, 492.899205, 499.110252], rtol=1e-8)
    np.testing.assert_allclose(legacy.q0[[0, 250, 1000, 1999]], [1528.369086, 16056.432179, 24.875189, -620.839308], rtol=1e-7)
    assert legacy.q0.mean() == pytest.approx(3266.982098, rel=1e-8)


def test_derivative_fix_is_exact_for_quadratic_signal(tmp_path):
    # SG order-2 derivative is exact for a quadratic, including the edges.
    t = np.arange(500) * 1e-3
    T1 = 300 + 80 * t + 15 * t**2
    T2 = 295 + 60 * t + 10 * t**2
    path = tmp_path / "tc.csv"
    np.savetxt(path, np.column_stack([t, T1, T2]), delimiter=",", header="time_s,T1_K,T2_K", comments="")
    r = hf.HeatFluxCalc(path, equations="cubic", depths_m=MM_DEPTHS, plot=False)
    np.testing.assert_allclose(r.T1_dot, 80 + 30 * t, rtol=0, atol=1e-8)
    np.testing.assert_allclose(r.T2_dot, 60 + 20 * t, rtol=0, atol=1e-8)
    np.testing.assert_allclose(r.T1_s, T1, rtol=1e-13)


@requires_originals("HeatFluxCalc.m", packages=("signal",))  # sgolayfilt
def test_live_octave_full_arrays(tmp_path, legacy):
    out = run_octave(tmp_path, "HeatFluxCalc;", ["t", "T1_s", "T2_s", "T1_dot", "T0", "q0"],
                     readmatrix_file=hf.DATA_FILE, packages=["signal"])  # fmt: skip
    np.testing.assert_array_equal(out["t"], legacy.t)
    np.testing.assert_allclose(out["T1_s"], legacy.T1_s, rtol=1e-13)
    np.testing.assert_allclose(out["T2_s"], legacy.T2_s, rtol=1e-13)
    T1d, T2d = legacy.T1_s / legacy.dt, legacy.T2_s / legacy.dt
    np.testing.assert_allclose(out["T1_dot"], T1d, rtol=1e-13)
    with pytest.warns(LegacyEquationWarning):
        T0 = hf.surface_temperature_legacy(legacy.T1_s, legacy.T2_s, T1d, T2d, *MATLAB_DEPTHS)
        q0 = hf.surface_heat_flux_legacy(legacy.T1_s, legacy.T2_s, T1d, T2d, *MATLAB_DEPTHS)
    np.testing.assert_allclose(out["T0"], T0, rtol=1e-12)
    np.testing.assert_allclose(out["q0"], q0, rtol=1e-9, atol=1e-8)
