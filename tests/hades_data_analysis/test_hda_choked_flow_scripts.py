"""PURPLE_MassFlow_from_PlenumPressure.m, cD_analysis.m, temp_analysis.m -> Python.

Golden numbers come from running the ORIGINAL .m files in GNU Octave 11.3
(random draws replayed through ``_octave_harness``).
"""

import math

import numpy as np
import pytest

from HADES.data_analysis import choked_nozzle_analysis as cna
from HADES.data_analysis import purpl_mass_flow_from_plenum_pressure as mf
from HADES.data_analysis.choked_flow import (
    choked_mass_flow,
    mass_flow_factor,
    mass_flow_parameter,
    stagnation_pressure_for_mass_flow,
)
from HADES.data_analysis.purpl_fake_data import PURPL_Fake_Plenum_pressure

from ._octave_harness import ReplayRNG, make_streams, octave_available, run_octave


# ---------------------------------------------------------------- choked_flow
def test_choked_flow_forms_are_equivalent():
    g, R, T0, P0, Cd, A = 1.4, 287.05, 293.15, 20e5, 0.62, 304.77e-6
    m1 = choked_mass_flow(Cd, A, g, P0 / (R * T0), P0)
    m2 = Cd * A * P0 * mass_flow_factor(g, R, T0)
    m3 = Cd * A * P0 * mass_flow_parameter(g) / math.sqrt(R * T0)
    assert m1 == pytest.approx(m2, rel=1e-14) and m1 == pytest.approx(m3, rel=1e-14)
    assert stagnation_pressure_for_mass_flow(m1, Cd, A, g, R, T0) == pytest.approx(P0, rel=1e-14)


# ------------------------------------------------------- mass flow from plenum
def test_mass_flow_hand_calculation_at_20_bar():
    P = np.full(4, 20e5)
    m_air, m_h2 = mf.plenum_mass_flow(P, P)
    # air: 0.62 * 304.77 mm^2 * 20 bar * sqrt(1.4/(287.05*293.15)) * (2/2.4)^3
    hand_air = 0.62 * 304.77e-6 * 20e5 * math.sqrt(1.4 / (287.05 * 293.15)) * (2 / 2.4) ** 3
    hand_h2 = 0.71 * (math.pi * 0.0005**2 * 32) * 20e5 * math.sqrt(1.405 / (4124.2 * 293.15)) * (2 / 2.405) ** (2.405 / 0.81)
    np.testing.assert_allclose(m_air, hand_air, rtol=1e-13)
    np.testing.assert_allclose(m_h2, hand_h2, rtol=1e-13)
    assert m_air[0] == pytest.approx(0.89202, rel=1e-4)  # kg/s  (1.967 lb/s)
    assert m_h2[0] == pytest.approx(0.022252, rel=1e-4)


def test_mass_flow_reads_rows_as_time_h2_air():
    # BUG FIX: MATLAB read row 2 as air and row 3 as H2; generator returns [time; h2; air].
    t = np.arange(5) / 1000
    pp = np.vstack([t, np.full(5, 10e5), np.full(5, 30e5)])  # H2 at 10 bar, air at 30 bar
    _, m_air, m_h2 = mf.PURPLE_MassFlow_from_PlenumPressure(pp, plot=False)
    exp_air, exp_h2 = mf.plenum_mass_flow(np.full(5, 30e5), np.full(5, 10e5))
    np.testing.assert_allclose(m_air, exp_air, rtol=1e-15)
    np.testing.assert_allclose(m_h2, exp_h2, rtol=1e-15)


def test_mass_flow_matches_octave_golden():
    pp = PURPL_Fake_Plenum_pressure(ReplayRNG(*make_streams(13)), plot=False)
    # MATLAB (swapped) row assignment reproduces the original numbers ...
    m_air_legacy, m_h2_legacy = mf.plenum_mass_flow(pp[1], pp[2])
    assert m_air_legacy.mean() == pytest.approx(0.8946632300234891, rel=1e-13)
    assert m_h2_legacy.mean() == pytest.approx(0.022282256919747348, rel=1e-13)
    assert m_air_legacy[500] == pytest.approx(0.8949237836948698, rel=1e-13)
    assert m_h2_legacy[500] == pytest.approx(0.022268961389643228, rel=1e-13)
    # ... and the fixed assignment gives the numbers recorded in PORTING_NOTES.md (lb/s)
    _, m_air, m_h2 = mf.PURPLE_MassFlow_from_PlenumPressure(pp, plot=False)
    assert m_air.mean() * mf.KG_TO_LB == pytest.approx(1.9693155524793753, rel=1e-12)
    assert m_h2.mean() * mf.KG_TO_LB == pytest.approx(0.04920066142095921, rel=1e-12)


def test_mass_flow_main_smoke(no_show):
    mf.main()


@pytest.mark.skipif(not octave_available("PURPLE_MassFlow_from_PlenumPressure.m", "PURPL_Fake_Plenum_pressure.m"),
                    reason="Octave or original .m not available")  # fmt: skip
def test_live_octave_mass_flow(tmp_path):
    streams = make_streams(29)
    out = run_octave(tmp_path, "PURPLE_MassFlow_from_PlenumPressure;",
                     ["time", "air_plenum_pressure", "h2_plenum_pressure", "mass_flow_air", "mass_flow_h2"],
                     streams=streams)  # fmt: skip
    pp = PURPL_Fake_Plenum_pressure(ReplayRNG(*streams), plot=False)
    np.testing.assert_array_equal(out["air_plenum_pressure"], pp[1])  # MATLAB's (mislabelled) air = row 2
    m_air_legacy, m_h2_legacy = mf.plenum_mass_flow(pp[1], pp[2])
    np.testing.assert_allclose(out["mass_flow_air"], m_air_legacy, rtol=1e-15)
    np.testing.assert_allclose(out["mass_flow_h2"], m_h2_legacy, rtol=1e-15)


# ------------------------------------------------------ cD / temp analysis
CD_GOLD = {
    "air": [24.798752384423537, 17.738863136233945, 13.777084658013075],
    "h2_phi1.0": [33.527228177198694, 23.982452937679867, 18.62623787622149],
    "h2_phi1.2": [40.23267381263844, 28.77894352521584, 22.351485451465795],
}
TEMP_GOLD = {
    "air": [17.397622465584174, 17.82048289978851, 18.237663825684486],
    "h2_phi1.0": [23.52110497747356, 24.09280060332218, 24.656817690829968],
    "h2_phi1.2": [28.22532597296827, 28.911360723986622, 29.588181228995964],
}
IDX = [0, 99, 199]


def test_cd_analysis_matches_octave_golden():
    s = cna.cD_analysis(plot=False)
    assert s.x[0] == 0.5 and s.x[-1] == 0.9 and len(s.x) == 200
    np.testing.assert_allclose(s.P0_air_bar[IDX], CD_GOLD["air"], rtol=1e-13)
    np.testing.assert_allclose(s.P0_h2_bar[2, IDX], CD_GOLD["h2_phi1.0"], rtol=1e-13)
    np.testing.assert_allclose(s.P0_h2_bar[4, IDX], CD_GOLD["h2_phi1.2"], rtol=1e-13)


def test_temp_analysis_matches_octave_golden():
    s = cna.temp_analysis(plot=False)
    assert s.x[0] == 273 and s.x[-1] == 300
    np.testing.assert_allclose(s.P0_air_bar[IDX], TEMP_GOLD["air"], rtol=1e-13)
    np.testing.assert_allclose(s.P0_h2_bar[2, IDX], TEMP_GOLD["h2_phi1.0"], rtol=1e-13)
    np.testing.assert_allclose(s.P0_h2_bar[4, IDX], TEMP_GOLD["h2_phi1.2"], rtol=1e-13)


def test_cd_analysis_hand_calculation():
    # Cd = 0.5, T0 = 283 K: P0 = mdot / (Cd A sqrt(gamma K/(R T0)))
    K = (2 / 2.4) ** 6
    P0 = 2 * 0.453592 / (0.5 * 0.472 * 6.4516e-4 * math.sqrt(1.4 * K / (287 * 283)))
    s = cna.cD_analysis(plot=False)
    assert s.P0_air_bar[0] == pytest.approx(P0 / 1e5, rel=1e-14)
    assert s.P0_air_bar[0] == pytest.approx(24.7988, rel=1e-5)
    # P0 ~ 1/Cd and H2 curves scale with phi
    np.testing.assert_allclose(s.P0_air_bar * s.x, s.P0_air_bar[0] * 0.5, rtol=1e-13)
    np.testing.assert_allclose(s.P0_h2_bar[:, 0] / s.P0_h2_bar[2, 0], s.phi_values, rtol=1e-13)


def test_temp_analysis_scales_with_sqrt_T():
    s = cna.temp_analysis(plot=False)
    np.testing.assert_allclose(s.P0_air_bar / np.sqrt(s.x), s.P0_air_bar[0] / np.sqrt(273), rtol=1e-13)


def test_nozzle_main_smoke(no_show):
    cna.main([])


@pytest.mark.skipif(not octave_available("cD_analysis.m", "temp_analysis.m"), reason="Octave or original .m not available")
@pytest.mark.parametrize("script", ["cD_analysis", "temp_analysis"])
def test_live_octave_nozzle(tmp_path, script):
    out = run_octave(tmp_path, f"{script};", [f"plot{k}y" for k in range(1, 7)] + ["plot1x"])
    s = getattr(cna, script)(plot=False)
    np.testing.assert_allclose(out["plot1x"], s.x, rtol=1e-14)
    np.testing.assert_allclose(out["plot1y"], s.P0_air_bar, rtol=1e-14)
    for i in range(5):
        np.testing.assert_allclose(out[f"plot{i + 2}y"], s.P0_h2_bar[i], rtol=1e-14)
