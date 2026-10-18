"""HADES sizing model: Python port vs the ORIGINAL MATLAB executed under Octave 11.3.

Reference numbers in reference/octave_*.json were produced by running the
unmodified HADES/SIZING MODEL/*.m files (harness: reference/octave_harness).
The CEA steps replayed the committed CEA2 outputs; the decks the MATLAB wrote
were byte-identical to the committed decks.
"""

from __future__ import annotations

import contextlib
import io
import math

import numpy as np
import pytest

from HADES.sizing_model import HADES_size_ceaDet as det_mod
from HADES.sizing_model import HADES_size_ceaRocket as rock_mod
from HADES.sizing_model.HADES_size_airInj import HADES_size_airInj
from HADES.sizing_model.HADES_size_annulusgap_fillheight import (
    HADES_size_annulusgap_fillheight,
    cellwidth_from_caltech_fit,
)
from HADES.sizing_model.HADES_size_chambCooltime import HADES_size_chambCooltime
from HADES.sizing_model.HADES_size_chambCooltime_v3 import HADES_size_chambCooltime_v3
from HADES.sizing_model.HADES_size_chamberPressure import HADES_size_chamberPressure
from HADES.sizing_model.HADES_size_geometry import HADES_size_geometry
from HADES.sizing_model.HADES_size_h2Inj import HADES_size_h2Inj
from HADES.sizing_model.HADES_size_HoopStressTemps import HADES_size_HoopStressTemps
from HADES.sizing_model.HADES_size_main import HADES_size_main
from HADES.sizing_model.HADES_size_max_flux_range import HADES_size_max_flux_range
from HADES.sizing_model.HADES_size_P0_v1 import HADES_size_P0_v1
from HADES.sizing_model.HADES_size_P_plenum_fanno import HADES_size_P_plenum_fanno
from HADES.sizing_model.HADES_size_pdrop_fanno import HADES_size_pdrop_fanno
from HADES.sizing_model.HADES_size_sonicnozzlearea_air import HADES_size_sonicnozzlearea_air
from HADES.sizing_model.HADES_size_thermal_v1 import HADES_size_thermal_v1
from HADES.sizing_model.HADES_size_throat import HADES_size_throat
from HADES.sizing_model.HADES_size_thrust import HADES_size_thrust, exit_thrust_from_cea

EXACT = 1e-12  # pure algebra: same operations as the MATLAB


def quiet(fn, *a, **k):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        out = fn(*a, **k)
    return out, buf.getvalue()


@pytest.fixture(scope="module")
def main_run():
    return quiet(HADES_size_main, plot=True)


# ---------------------------------------------------------------- main chain
def test_main_scalars_match_octave(main_run, octave_main):
    r, _ = main_run
    o = octave_main
    for name, attr in [("cellWidth", "cellWidth"), ("annulus_gap", "annulus_gap"),
                       ("chosen_fillHeight", "chosen_fillHeight"),
                       ("det_wave_path_length", "det_wave_path_length"),
                       ("avg_chamber_p", "avg_chamber_p"), ("Failure_temps", "Failure_temps"),
                       ("P_0", "P_0"), ("h2_area", "h2_area"), ("air_area", "air_area"),
                       ("area_annulus", "area_annulus"), ("max_flux", "max_flux"),
                       ("fuel_mdot", "fuel_mdot"), ("total_mdot", "total_mdot"),
                       ("thrust", "thrust")]:
        assert getattr(r, attr) == pytest.approx(o[name], rel=EXACT), name
    # fsolve/fzero based: Octave's fsolve vs scipy's fsolve
    assert r.P_plenum_req == pytest.approx(o["P_plenum_req"], rel=1e-7)


def test_main_chamber_pressure_reproduces_committed_deck(main_run):
    r, _ = main_run
    assert "%f" % r.avg_chamber_p == "142.745149"  # p,psia in committed cea_run.inp


def test_main_cea_det_results_match_matlab_parser(main_run, octave_main):
    d = main_run[0].ceaDet_results
    for k, v in octave_main["ceaDet"].items():
        if k == "k_eq":
            assert d.k_eq == pytest.approx(v, abs=1.0001e-4)  # CEA3 5.2102 vs CEA2 5.2103
        else:
            assert getattr(d, k) == v, k


def test_main_cea_rocket_legacy_fields_are_throat_like_matlab(main_run, octave_main):
    r, _ = main_run
    rr, o = r.ceaRock_results, octave_main["ceaRock"]
    # legacy fields keep the MATLAB parser meaning (throat column of cea_run.out)
    assert (rr.AeAt, rr.cf, rr.isp, rr.ivac, rr.cstar) == (o["AeAt"], o["cf"], o["isp"], o["ivac"], o["cstar"])
    assert (rr.AeAt, rr.cf, rr.isp, rr.ivac) == (1.0, 0.6884, 927.8, 1679.1)
    assert rr.Pc == pytest.approx(o["Pc"], rel=EXACT)
    assert (rr.AeAt_throat, rr.cf_throat, rr.isp_throat, rr.ivac_throat) == (rr.AeAt, rr.cf, rr.isp, rr.ivac)
    # exit column of cea_run.out, unambiguous names
    assert (rr.AeAt_exit, rr.cf_exit, rr.isp_exit, rr.ivac_exit, rr.pe_exit_bar) == (10.0, 1.6019, 2159.0, 2293.8, 0.09844)
    assert rr.pc_bar == pytest.approx(142.745149 * 1.01325 / 14.696006, rel=EXACT)
    # main keeps the legacy throat-based thrust
    assert r.thrust == pytest.approx(927.8 * r.total_mdot * 0.453592, rel=EXACT)


def test_exit_thrust_option_matches_nasa_thrust_equation(main_run):
    """Opt-in exit-plane thrust: F = mdot Ve + (pe - p0) Ae (NASA rocket thrust equation)."""
    r, _ = main_run
    raw = r.ceaRock_results.raw
    mdot = r.total_mdot * 0.453592  # kg/s
    # hand calculation from the cea_run.out exit column (printed values)
    Ve, pe, cstar, eps, pc = 2159.0, 0.09844e5, 1347.8, 10.0, 9.8419e5
    Ae = eps * mdot * cstar / pc
    p0 = 101325.0
    hand = mdot * Ve + (pe - p0) * Ae
    F, out = quiet(exit_thrust_from_cea, r.ceaRock_results, r.total_mdot, p0)
    assert F == pytest.approx(hand, rel=2e-4)  # printed-precision inputs
    assert F == pytest.approx(843.6, abs=0.5)  # over-expanded at sea level: below the 863.8 N legacy value
    assert out.startswith("Thrust: ")
    # CEA identity Ivac = Ve + pe Ae / mdot -> vacuum thrust = mdot * Ivac
    Fvac, _ = quiet(exit_thrust_from_cea, r.ceaRock_results, r.total_mdot, 0.0)
    assert Fvac == pytest.approx(mdot * raw.ivac[-1], rel=1e-6)
    # direct call with explicit arguments
    F2, _ = quiet(HADES_size_thrust, Ve, r.total_mdot, station="exit", exit_pressure=pe, ambient_pressure=p0,
                  exit_area=Ae)
    assert F2 == pytest.approx(hand, rel=EXACT)
    with pytest.raises(ValueError):
        HADES_size_thrust(Ve, r.total_mdot, station="exit", exit_pressure=pe)
    with pytest.raises(ValueError):
        HADES_size_thrust(Ve, r.total_mdot, exit_area=Ae)  # exit inputs need station='exit'
    with pytest.raises(ValueError):
        HADES_size_thrust(Ve, r.total_mdot, station="nozzle")


def test_main_stdout_identical_to_octave(main_run, ref_dir):
    _, out = main_run
    py = out.rstrip("\n").split("\n")
    oc = (ref_dir / "octave_main_v1_stdout.txt").read_text().rstrip("\n").split("\n")
    assert py == oc
    assert "Thrust: 863.770 N (194.176 lbf) " in py


def test_regenerated_decks_are_byte_identical_to_committed(ref_dir):
    d = det_mod.cea2_input_deck(ox="Air", fuel="H2", phi=0.9, P0=26, P0Units="psia", T0=283.15, T0Units="K")
    r = rock_mod.cea2_input_deck(ox="Air", fuel="H2", phi=0.9, Pc=142.74514850616087, PcUnits="psia")
    assert d == (ref_dir / "cea_det.inp").read_text()
    assert r == (ref_dir / "cea_run.inp").read_text()


def test_cea_wrappers_matlab_style_arguments():
    a = det_mod.HADES_size_ceaDet("ox", "Air", "fuel", "H2", "phi", 0.9, "P0", 26, "P0Units", "psia",
                                  "T0", 283.15, "T0Units", "K")
    b = det_mod.HADES_size_ceaDet(ox="Air", fuel="H2", phi=0.9, P0=26, T0=283.15)
    assert a.cjVel == b.cjVel == 1929.7
    # unit handling: same state given in bar and C
    c = det_mod.HADES_size_ceaDet(ox="Air", fuel="H2", phi=0.9, P0=26 * 6894.76 / 1e5, P0Units="bar",
                                  T0=10.0, T0Units="C")
    assert c.cjVel == a.cjVel and c.P_ratio == a.P_ratio
    with pytest.raises(ValueError):
        det_mod.HADES_size_ceaDet(ox="Air", fuel="H2", phi=0.9, of=30, P0=26, T0=283.15)
    with pytest.raises(ValueError):
        rock_mod.HADES_size_ceaRocket(ox="Air", fuel="H2", phi=0.9, Pc=140, Pe=14.7, areaRatio=8)
    with pytest.warns(UserWarning):
        det_mod.HADES_size_ceaDet(ox="Air", fuel="H2", phi=0.9, P0=26, T0=283.15, ceaExe="FCEA2.exe")
    # inputParser semantics: case-insensitive, exact match first, then unique prefix
    d = det_mod.HADES_size_ceaDet("OX", "Air", "fu", "H2", "ph", 0.9, "p0", 26, "p0u", "psia", "t0", 283.15)
    assert d.cjVel == a.cjVel
    with pytest.raises(TypeError, match="ambiguous"):
        det_mod.HADES_size_ceaDet("o", "Air", "fuel", "H2", "phi", 0.9, "P0", 26, "T0", 283.15)  # ox / of
    with pytest.raises(TypeError, match="ambiguous"):
        rock_mod.HADES_size_ceaRocket(ox="Air", fuel="H2", phi=0.9, P=140)  # Pc / PcUnits / Pe / PeUnits
    with pytest.raises(TypeError, match="unknown"):
        det_mod.HADES_size_ceaDet(ox="Air", fuel="H2", phi=0.9, P0=26, T0=283.15, oxidiser="O2")
    # Pe option -> pi/p deck line, exit pressure honoured
    pe = rock_mod.HADES_size_ceaRocket(ox="Air", fuel="H2", phi=0.9, Pc=142.745149, Pe=14.7)
    assert pe.raw.P_bar[-1] == pytest.approx(14.7 * 0.0689475729, rel=1e-5)
    assert pe.pe_exit_bar == pytest.approx(1.0135, abs=5e-5)
    assert pe.isp == 927.8  # legacy field = throat, independent of the exit condition


# ------------------------------------------------------- individual functions
def test_individual_functions_match_octave(octave_units):
    o = octave_units
    assert HADES_size_throat(0.93, 9.8e5, 2330, 1.22, 333) == pytest.approx(o["throat"][0], rel=EXACT)
    assert HADES_size_throat(1, 2e6, 300, 1.4, 287) == pytest.approx(o["throat"][1], rel=EXACT)
    for key, args in [("pdrop1", (20e5, 4e-3, 0.1, 3.2e-6, 0.01, 1.4, 1.8e-5)),
                      ("pdrop2", (40e5, 6e-3, 0.3, 1e-5, 0.02, 1.4, 8.8e-6))]:
        np.testing.assert_allclose(HADES_size_pdrop_fanno(*args), o[key], rtol=1e-10)
    for key, args in [("hoop1", (2, 2.5, 0.101325, 30.7, 1.5)), ("hoop2", (1.375, 2, 0.101325, 45, 2.0))]:
        np.testing.assert_allclose(quiet(HADES_size_HoopStressTemps, *args, plot=False)[0], o[key], rtol=1e-10)
    p2, _ = quiet(HADES_size_P_plenum_fanno, 8e5, 5e-4, math.pi * 0.25 * 0.0015**2, 8, 1.5, 290, 0.65, 1.4,
                  4124.2, 84e-7, 273, 72, 3.2e-3, 0.4)
    assert p2 == pytest.approx(o["plenum2"], rel=1e-7)
    f2, _ = quiet(HADES_size_annulusgap_fillheight, 1.1, 0.25, 2, 1980, 1100, 30, 3, 420, plot=False)
    np.testing.assert_allclose(f2, o["fill2"], rtol=EXACT)
    assert quiet(HADES_size_chamberPressure, 2, 10.0, 1950.3, 15.2, 30, plot=False)[0] == pytest.approx(o["chp2"], rel=EXACT)
    assert quiet(HADES_size_h2Inj, 0.06, 150)[0] == pytest.approx(o["h2"], rel=EXACT)
    np.testing.assert_allclose(quiet(HADES_size_airInj, 1.8, 8, 2.2, 290)[0], o["air"], rtol=EXACT)
    assert quiet(HADES_size_geometry, 2.5, 0.5, 8)[0] == pytest.approx(o["geom"], rel=EXACT)
    assert quiet(HADES_size_thrust, 1500, 2.1)[0] == pytest.approx(o["thrust"], rel=EXACT)
    assert quiet(HADES_size_P0_v1)[0] == pytest.approx(o["P0"], rel=EXACT)
    assert quiet(HADES_size_max_flux_range, 2.1, 2000)[0] == pytest.approx(o["flux"], rel=EXACT)


def test_sonic_nozzle_matches_octave_and_documented_R_fix(octave_units):
    p, d_orig, _ = HADES_size_sonicnozzlearea_air(R=8.134, plot=False)  # original typo
    np.testing.assert_array_equal(p, octave_units["sonic_p"])
    np.testing.assert_allclose(d_orig, octave_units["sonic_dia_mm"], rtol=EXACT)
    _, d_fixed, _ = HADES_size_sonicnozzlearea_air(plot=True)  # R = 8.314
    np.testing.assert_allclose(d_fixed / d_orig, (8.314 / 8.134) ** 0.25, rtol=EXACT)


def test_superseded_cell_size2_logic_preserved(octave_units):
    assert cellwidth_from_caltech_fit(0.9) == pytest.approx(octave_units["cell_size2_cellwidth"], rel=1e-10)


# ------------------------------------------------------------ ODE / PDE scripts
def test_cooltime_v1_matches_octave(octave_cool):
    r, _ = quiet(HADES_size_chambCooltime, plot=True)
    o = octave_cool["v1"]
    # first-output-below-target time: solver output grids differ (Octave's ode45
    # is not MATLAB's), spacing ~3 s -> compare within 1.5 s (2e-4 relative)
    assert r.tCool95_minutes * 60 == pytest.approx(o["v1_tCool95"], abs=1.5)
    assert r.T_wall[-1] - 273.15 == pytest.approx(o["v1_Twall_end_C"], abs=1e-4)
    # air node (tau ~ 0.9 s) is integrated at ode45's stability limit: step-sequence dependent
    assert r.T_air[-1] - 273.15 == pytest.approx(o["v1_Tair_end_C"], abs=0.01)
    assert r.T95_wall_C == pytest.approx(o["v1_T95_wall_C"], abs=0.01)
    assert r.t[-1] == o["v1_t_end"] == 10800


def test_cooltime_v3_matches_octave_and_v2_is_identical(octave_cool):
    o3, o2 = octave_cool["v3"], octave_cool["v2"]
    assert o2 == o3  # v2 and v3 give identical numbers in Octave -> v2 superseded
    out, _ = quiet(HADES_size_chambCooltime_v3, plot=True)
    for i, key in enumerate(["base", "N2", "rough", "combo"]):
        s = out[key]
        assert s.t_95 == pytest.approx(o3["t95"][i], abs=1.5), key
        assert s.t_90 == pytest.approx(o3["t90"][i], abs=1.5), key
        assert s.T_wall[-1] == pytest.approx(o3["Twall_end"][i], abs=1e-4), key


@pytest.mark.parametrize("bad", [0, -1, 5, 2.0, True, "2", None])
def test_thermal_v1_rejects_bad_phi_index(bad):
    with pytest.raises(ValueError):
        HADES_size_thermal_v1(bad, plot=False)


def test_thermal_v1_matches_octave(octave_cool):
    r, _ = quiet(HADES_size_thermal_v1, plot=True)
    o = octave_cool["thermal"]
    assert r.h_inner == pytest.approx(o["h_inner"], rel=EXACT)
    assert r.T_rec_g == pytest.approx(o["T_rec_g"], rel=EXACT)
    assert r.Re_g == pytest.approx(o["Re_g"], rel=EXACT)
    assert r.Nu_g == pytest.approx(o["Nu_g"], rel=EXACT)
    # pdepe: Python (Skeel-Berzins + scipy BDF) vs Octave shim (Skeel-Berzins + ode15s)
    np.testing.assert_allclose(r.time_to_limit, o["time_to_limit"], rtol=5e-4)
    assert np.max(np.abs(r.T_inner - np.array(o["T_inner"]))) < 0.5  # C
    assert np.max(np.abs(r.T[-1] - np.array(o["T_end_profile"]))) < 0.1  # C
