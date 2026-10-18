"""Regression tests for HADES.feed_sizing.sonic_nozzle_sizing (SonicNozzleSizingV1-V5.m).

References (tests/hades_feed_sizing/reference/):
  matlab/  result CSVs committed to the repo, written by MATLAB from the
           original V2-V5 scripts (15 significant digits).
  octave/  the original, unmodified V1-V5 .m scripts run under GNU Octave 11.3
           (table/readtable/writetable shimmed), CSVs written with 17 significant
           digits, plus each script's command-window output.
"""

import contextlib
import io
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from HADES.feed_sizing import sonic_nozzle_sizing as sns

REF = Path(__file__).parent / "reference"

CASES = {
    "v2": (sns.size_fixed_sonic_orifice, "flow_requirements_sonic_orifice_results.csv"),
    "v3": (sns.size_two_stage_choked, "flow_requirements_two_stage_choked_results.csv"),
    "v4": (sns.manifold_pressure_from_injector, "manifold_pressure_results.csv"),
    "v5": (sns.size_injector_sonic_orifice, "flow_requirements_manifold_input_combined_results.csv"),
}


def _compute(model):
    func, _ = CASES[model]
    with pytest.warns() if model == "v5" else contextlib.nullcontext():
        return func()


def _assert_table_matches(table: pd.DataFrame, ref: pd.DataFrame, rtol: float):
    assert list(table.columns) == list(ref.columns)
    assert len(table) == len(ref)
    for column in ref.columns:
        got = table[column].to_numpy(dtype=float)
        exp = ref[column].to_numpy(dtype=float)
        np.testing.assert_array_equal(np.isnan(got), np.isnan(exp), err_msg=column)
        if table[column].dtype == bool:
            np.testing.assert_array_equal(got, exp, err_msg=column)
        else:
            np.testing.assert_allclose(got, exp, rtol=rtol, atol=0, equal_nan=True, err_msg=column)


@pytest.mark.parametrize("model", sorted(CASES))
def test_matches_committed_matlab_results(model):
    """Re-run on the committed inputs; compare with the committed MATLAB CSVs (15 sig. digits)."""
    result = _compute(model)
    ref = pd.read_csv(REF / "matlab" / CASES[model][1], float_precision="round_trip")
    _assert_table_matches(result.table, ref, rtol=1e-13)


@pytest.mark.parametrize("model", sorted(CASES))
def test_matches_octave_full_precision(model):
    """Compare with the original .m run in Octave, written at 17 significant digits."""
    result = _compute(model)
    ref = pd.read_csv(REF / "octave" / CASES[model][1], float_precision="round_trip")
    _assert_table_matches(result.table, ref, rtol=1e-14)


@pytest.mark.parametrize("model", ["v1", "v2", "v3", "v4", "v5"])
def test_printed_summary_matches_octave(model, tmp_path):
    out = None if model == "v1" else tmp_path / (CASES[model][1])
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), pytest.warns() if model == "v5" else contextlib.nullcontext():
        sns.run(model, output_file=out)
    text = buf.getvalue()
    if out is not None:
        # MATLAB printed the absolute output path (V4: the bare file name).
        text = text.replace(str(out), out.name if model == "v4" else f"<OUTPUT_DIR>/{out.name}")
    expected = (REF / "octave" / f"SonicNozzleSizing{model.upper()}_stdout.txt").read_text()
    assert text == expected


@pytest.mark.parametrize("model", sorted(CASES))
def test_written_csv_matches_matlab_format(model, tmp_path):
    out = tmp_path / "out.csv"
    with pytest.warns() if model == "v5" else contextlib.nullcontext():
        sns.run(model, output_file=out, verbose=False)
    text = out.read_text()
    ref_text = (REF / "matlab" / CASES[model][1]).read_text()
    assert text.splitlines()[0] == ref_text.splitlines()[0]
    # logicals as 1/0 and missing values as NaN, like MATLAB writetable
    assert ("NaN" in text) == ("NaN" in ref_text)
    assert "True" not in text and "False" not in text
    ref = pd.read_csv(REF / "matlab" / CASES[model][1], float_precision="round_trip")
    got = pd.read_csv(out, float_precision="round_trip")
    pd.testing.assert_frame_equal(got, ref, check_dtype=False, rtol=1e-13, atol=0)


def test_v1_matches_octave_values():
    result = sns.size_single_sonic_nozzle()
    ref = pd.read_csv(REF / "octave" / "SonicNozzleSizingV1_values.csv", float_precision="round_trip")
    for name, value in zip(ref["name"], ref["value"]):
        assert getattr(result, name) == pytest.approx(value, rel=1e-14), name


def test_v1_hand_calculation():
    """Independent textbook form: mdot = Cd*A*P0*sqrt(gamma/(R*T0))*(2/(gamma+1))^((gamma+1)/(2(gamma-1)))."""
    r = sns.size_single_sonic_nozzle()
    psi = 6894.76
    P_down = 20e5 + 20 * psi  # 21.378952 bar
    gamma, R, T0 = 1.4, 287.05, 283.0
    P0 = P_down / 0.5282817877171742  # (2/2.4)^3.5
    mdot = 0.9072 * 1.15
    flux = P0 * math.sqrt(gamma / (R * T0)) * (2 / 2.4) ** 3.0
    D = math.sqrt(4 * (mdot / flux) / math.pi)
    assert r.P_tu_min == pytest.approx(P0, rel=1e-12)
    assert r.D_t_m == pytest.approx(D, rel=1e-12)
    assert r.D_t_m * 1000 == pytest.approx(11.689, abs=5e-4)  # printed "11.689 mm"
    assert r.T_throat == pytest.approx(283 / 1.2, rel=1e-15)
    assert r.feed_ID / 0.0254 == pytest.approx(0.902, rel=1e-12)


def test_v2_sizing_point_and_round_trip():
    r = sns.size_fixed_sonic_orifice()
    data = pd.read_csv(sns.FLOW_REQUIREMENTS_CSV)
    for g, gas in enumerate(sns.GAS_NAMES):
        props = sns.GAS_PROPERTIES[gas]
        m = data[gas].to_numpy(float)
        valid = ~np.isnan(m)
        # choked flow through the sized orifice at the reported P0 returns the requested flow
        flux = (
            r.P0_bar_all[valid, g] * 1e5
            * math.sqrt(props.gamma / (props.R * props.T0))
            * (2 / (props.gamma + 1)) ** ((props.gamma + 1) / (2 * (props.gamma - 1)))
        )
        np.testing.assert_allclose(r.Cd_all[g] * r.orifice_area[g] * flux, m[valid], rtol=1e-12)
        # at the minimum flow the upstream pressure is exactly the design minimum
        i_min = np.nanargmin(m)
        assert r.P0_bar_all[i_min, g] * 1e5 == pytest.approx(r.P0_design_min_all[g], rel=1e-12)
        assert r.choked_all[valid, g].all()
    assert r.P_down_all / 1e5 == pytest.approx([21.378952, 21.378952, 25.378952, 37.378952])


def test_v3_two_stage_invariants():
    r = sns.size_two_stage_choked()
    P_ch = 10e5
    for g in range(4):
        crit = r.crit_ratio_all[g]
        assert r.Pman_design_min_all[g] == pytest.approx(1.15 * P_ch / crit, rel=1e-14)
        assert r.P0_orifice_design_min_all[g] == pytest.approx(1.15 * r.Pman_design_min_all[g] / crit, rel=1e-14)
        valid = ~np.isnan(r.P0_bar_all[:, g])
        # both stages scale linearly with mdot, so the orifice pressure ratio is constant
        ratio = r.P_manifold_bar_all[valid, g] / r.P0_bar_all[valid, g]
        np.testing.assert_allclose(ratio, crit / 1.15, rtol=1e-12)


def test_v5_warnings_and_validity_counts():
    with pytest.warns(RuntimeWarning) as record:
        r = sns.size_injector_sonic_orifice()
    messages = [str(w.message) for w in record]
    assert any(m.startswith("H2 injector is NOT choked at 1 test point(s).") for m in messages)
    assert any(m.startswith("CH4 injector is NOT choked at 5 test point(s).") for m in messages)
    assert r.model_valid_all.sum(axis=0).tolist() == [4, 5, 11, 16]
    assert r.valid_all.sum(axis=0).tolist() == [5, 5, 16, 16]
    # sonic orifice sized so that P0 = margin * Pman / crit at the minimum flow
    np.testing.assert_allclose(
        r.P0_design_minflow_all, 1.15 * r.Pman_minflow_all / r.critical_ratio_all, rtol=1e-14
    )


def test_v4_matches_v5_manifold_pressure():
    v4 = sns.manifold_pressure_from_injector()
    with pytest.warns(RuntimeWarning):
        v5 = sns.size_injector_sonic_orifice()
    np.testing.assert_array_equal(v4.Pman_bar_all, v5.Pman_bar_all)
    np.testing.assert_array_equal(v4.choked_all, v5.injector_choked_all)


def test_v4_keeps_matlab_semantics_for_zero_flow():
    data = pd.DataFrame({"H2": [0.0, 0.02], "Air": [0.9, np.nan], "CH4": [0.1, 0.1], "GOx": [0.5, 0.5]})
    r = sns.manifold_pressure_from_injector(data)
    assert r.Pman_bar_all[0, 0] == 0.0  # V4 has no positive-flow filter
    assert not r.choked_all[0, 0]
    assert np.isnan(r.Pman_bar_all[1, 1]) and not r.choked_all[1, 1]


def test_missing_column_raises(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text("H2,Air,CH4\n0.1,0.2,0.3\n")
    with pytest.raises(ValueError, match='CSV must contain a column named "GOx"'):
        sns.read_flow_requirements(bad)


def test_case_insensitive_columns_and_empty_gas_warning():
    data = pd.DataFrame({"h2": [0.02, 0.03], "AIR": [np.nan, np.nan], "ch4": [0.1, 0.2], "gox": [0.5, 0.4]})
    with pytest.warns(RuntimeWarning, match="No positive flow rates found for Air."):
        r = sns.size_fixed_sonic_orifice(data)
    assert np.isnan(r.orifice_diameter_m[1])
    assert "nan" not in r.summary() and "NaN" in r.summary()


def test_overrides_and_defaults():
    base = sns.size_fixed_sonic_orifice()
    bigger_cd = sns.size_fixed_sonic_orifice(Cd={"Air": 1.0})
    assert bigger_cd.orifice_area[1] == pytest.approx(base.orifice_area[1] * 0.85, rel=1e-14)
    assert bigger_cd.orifice_area[0] == base.orifice_area[0]
    with pytest.raises(ValueError, match="Unknown gas"):
        sns.size_fixed_sonic_orifice(Cd={"N2": 0.9})


def test_inputs_module_relative_outputs_default_to_cwd(tmp_path, monkeypatch):
    module_dir = Path(sns.__file__).resolve().parent
    assert sns.FLOW_REQUIREMENTS_CSV.parent == module_dir
    assert sns.MANIFOLD_INPUT_CSV.parent == module_dir
    monkeypatch.chdir(tmp_path)
    expected = {
        "v2": "flow_requirements_sonic_orifice_results.csv",
        "v3": "flow_requirements_two_stage_choked_results.csv",
        "v4": "manifold_pressure_results.csv",
        "v5": "flow_requirements_manifold_input_combined_results.csv",
    }
    for model, name in expected.items():
        default_input = sns.MODELS[model][1]
        assert sns.default_output_path(model, default_input) == tmp_path / name
    # an input elsewhere still writes into the CWD
    assert sns.default_output_path("v2", "/some/where/my_flows.csv") == tmp_path / "my_flows_sonic_orifice_results.csv"


def test_default_run_writes_cwd_and_leaves_bundled_results(tmp_path, monkeypatch):
    bundled = {
        name: (sns.MODULE_DIR / name).read_bytes() for _, name in CASES.values()
    }
    monkeypatch.chdir(tmp_path)
    for model, (_, name) in CASES.items():
        with pytest.warns() if model == "v5" else contextlib.nullcontext():
            sns.run(model, verbose=False)
        assert (tmp_path / name).read_bytes() == (REF / "matlab" / name).read_bytes()
    for name, content in bundled.items():
        assert (sns.MODULE_DIR / name).read_bytes() == content


def test_cli_no_write(capsys, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    sns.main(["v3", "--no-write"])
    assert "TWO-STAGE CHOKED FLOW SIZING SUMMARY" in capsys.readouterr().out
    assert list(tmp_path.iterdir()) == []
    out = tmp_path / "sub" / "x.csv"
    out.parent.mkdir()
    sns.main(["v2", "--output", str(out)])
    assert out.exists()
    assert f"Results written to:\n{out}\n" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# Feed-line prototype (SonicNozzleSizing.py), area-Mach bug corrected
# ---------------------------------------------------------------------------

#: NACA 1135 / Anderson isentropic table, gamma = 1.4: A/A* (subsonic)
AREA_RATIO_TABLE = {0.1: 5.8218, 0.2: 2.9635, 0.3: 2.0351, 0.4: 1.5901, 0.5: 1.3398, 0.6: 1.1882, 0.8: 1.0382}


def _area_ratio(M, g):
    """Textbook isentropic A/A*."""
    return (1 / M) * ((2 / (g + 1)) * (1 + (g - 1) / 2 * M**2)) ** ((g + 1) / (2 * (g - 1)))


def _prototype_buggy(m_dot, wall):
    """Equations of the ORIGINAL SonicNozzleSizing.py, including its A_factor bug."""
    in_to_m, P_m, R, g, T, margin = 0.0254, 2e6, 287, 1.4, 283, 0.85
    D = (16 / 16 - wall * 2) * in_to_m
    A_2 = math.pi * D**2 / 4
    M_2 = (m_dot * margin) / (P_m * A_2) * math.sqrt(R * T / g)
    A_factor = 2 / (g + 1) * (1 + (g - 1) / 2 * M_2**2) ** (0.5 * (g + 1) / (g - 1))
    A_t = A_2 / A_factor * M_2
    T_t = T * (1 + (g - 1) / 2 * M_2**2)
    P_t_u = m_dot / A_t * math.sqrt(T_t * R / g) * ((g + 1) / 2) ** ((g + 1) / (2 * (g - 1)))
    return {"M_2": M_2, "A_2": A_2, "A_factor": A_factor, "A_t": A_t, "D_t": 2 * math.sqrt(A_t / math.pi), "P_t_u": P_t_u}


def test_feedline_matches_authors_corrected_matlab_in_octave():
    """f8c4cd3 SonicNozzleSizing.m (the authors' corrected translation, 2 lb/s) run in Octave."""
    r = sns.size_sonic_nozzle_with_feed_line(m_dot=2 * sns.LB_TO_KG)
    ref = pd.read_csv(REF / "octave" / "SonicNozzleSizing_f8c4cd3_values.csv", float_precision="round_trip")
    for name, value in zip(ref["name"], ref["value"]):
        assert getattr(r, name) == pytest.approx(value, rel=1e-14), name


@pytest.mark.parametrize("M", sorted(AREA_RATIO_TABLE))
def test_feedline_area_ratio_matches_textbook_table(M):
    # choose the flow that gives feed-line Mach M for the default 0.810 in ID line
    base = sns.size_sonic_nozzle_with_feed_line()
    r = sns.size_sonic_nozzle_with_feed_line(m_dot=base.m_dot * M / base.M_2)
    assert r.M_2 == pytest.approx(M, rel=1e-12)
    assert r.A_2 / r.A_t == pytest.approx(AREA_RATIO_TABLE[M], rel=1e-4)
    assert r.A_2 / r.A_t == pytest.approx(_area_ratio(M, 1.4), rel=1e-13)


def test_feedline_independent_equations_and_root_round_trip():
    from scipy.optimize import brentq

    r = sns.size_sonic_nozzle_with_feed_line()
    g, R, T, P_m = r.gamma, r.R, r.T, r.P_m
    # continuity at static P_m, T: mdot*margin = P*A*M*sqrt(g/(R T))
    assert r.m_dot * r.mass_flow_margin == pytest.approx(P_m * r.A_2 * r.M_2 * math.sqrt(g / (R * T)), rel=1e-13)
    assert r.feed_diameter / 0.0254 == pytest.approx(0.810, rel=1e-13)
    # subsonic root of A/A*(M) = A_2/A_t recovers the feed-line Mach number
    M = brentq(lambda m: _area_ratio(m, g) - r.A_2 / r.A_t, 1e-6, 1.0, xtol=1e-15)
    assert M == pytest.approx(r.M_2, rel=1e-10)
    assert r.T_t == pytest.approx(T * (1 + 0.2 * r.M_2**2), rel=1e-15)
    # choked flow of m_dot (Cd = 1) through A_t at (P_t_u, T_t)
    mdot = r.A_t * r.P_t_u * math.sqrt(g / (R * r.T_t)) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)))
    assert mdot == pytest.approx(r.m_dot, rel=1e-13)
    assert r.P_t_u_req == pytest.approx(P_m / 0.5282817877171742, rel=1e-15)


def test_feedline_bug_fix_before_after_numbers():
    lb = 0.453592
    # the buggy equations reproduce the original SonicNozzleSizing.py printout exactly
    buggy = _prototype_buggy(3 * lb, 0.095)
    printed = (REF / "original_python" / "SonicNozzleSizing_py_stdout.txt").read_text()
    assert f"M_2 = {buggy['M_2']}\n" in printed
    assert f"A_t = {buggy['A_t']} m^2\n" in printed
    assert f"Throat Diameter: {buggy['D_t'] * 1000} mm\n" in printed
    assert f"Upstream Pressure: {buggy['P_t_u'] * 10 ** (-5)} bar\n" in printed
    assert "Flow is choked!" in printed

    fixed = sns.size_sonic_nozzle_with_feed_line()
    assert fixed.M_2 == buggy["M_2"]  # Mach number is unaffected by the fix
    assert buggy["A_factor"] == pytest.approx(0.924234, abs=5e-7)
    assert fixed.A_factor == pytest.approx(0.641829, abs=5e-7)
    assert buggy["D_t"] * 1e3 == pytest.approx(13.852781, abs=5e-7)
    assert fixed.D_t * 1e3 == pytest.approx(16.623337, abs=5e-7)
    assert buggy["P_t_u"] / 1e5 == pytest.approx(38.232314, abs=5e-7)
    assert fixed.P_t_u / 1e5 == pytest.approx(26.550218, abs=5e-7)
    # the difference is exactly the A_factor ratio
    assert fixed.A_t / buggy["A_t"] == pytest.approx(buggy["A_factor"] / fixed.A_factor, rel=1e-13)
    assert buggy["P_t_u"] / fixed.P_t_u == pytest.approx(fixed.A_t / buggy["A_t"], rel=1e-13)
    # the fix flips the (legacy) pressure-ratio check on the prototype's default inputs
    assert fixed.P_t_u_req / 1e5 == pytest.approx(37.858583, abs=5e-7)
    assert not fixed.legacy_pressure_check_passed
    assert "Legacy pressure-ratio check (Upstream Pressure > threshold): not passed\n" in fixed.summary()
    # throat now depends on feed-line geometry (wall 0.049 in)
    thin = sns.size_sonic_nozzle_with_feed_line(wall_thickness_in=0.049)
    assert thin.M_2 == pytest.approx(0.337890, abs=5e-7)
    assert _prototype_buggy(3 * lb, 0.049)["D_t"] * 1e3 == pytest.approx(14.102979, abs=5e-7)
    assert thin.D_t * 1e3 == pytest.approx(16.923575, abs=5e-7)


def test_feedline_summary_format_and_cli(capsys):
    r = sns.size_sonic_nozzle_with_feed_line()
    printed = (REF / "original_python" / "SonicNozzleSizing_py_stdout.txt").read_text()
    # sizing lines keep the original prototype layout (only the numbers changed)
    labels = [line.split("=")[0].split(":")[0] for line in printed.splitlines()]
    assert [line.split("=")[0].split(":")[0] for line in r.summary().splitlines()][:6] == labels[:6]
    # the old choke verdict is relabelled as a legacy diagnostic, with its limitation
    assert "Flow is choked" not in r.summary() and "Flow does not choke" not in r.summary()
    assert f"Legacy threshold P_m / critical ratio: {r.P_t_u_req * 10 ** (-5)} bar\n" in r.summary()
    assert "not a valid choking criterion" in r.summary()
    sns.main(["feedline"])
    out = capsys.readouterr().out
    assert out == r.summary()
    assert "NOTE: legacy diagnostic from SonicNozzleSizing.py, not a valid choking criterion" in out


def test_feedline_legacy_check_is_not_a_choking_criterion():
    """Counterexample: with mass_flow_margin = 1 the model is exactly sonic at the throat
    (A_t = A* from independent sonic-state continuity, P_t_u = isentropic total pressure
    of the feed-line state), yet the legacy check reports 'not passed'."""
    r = sns.size_sonic_nozzle_with_feed_line(mass_flow_margin=1.0)
    g, R = r.gamma, r.R
    P0 = r.P_m * (1 + (g - 1) / 2 * r.M_2**2) ** (g / (g - 1))  # isentropic total pressure
    A_star = r.m_dot * math.sqrt(r.T_t) / (P0 * math.sqrt(g / R) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1))))
    assert r.M_2 == pytest.approx(0.492946, abs=5e-7)
    assert r.P_t_u / 1e5 == pytest.approx(23.613659, abs=5e-7)
    assert r.P_t_u == pytest.approx(P0, rel=1e-13)
    assert r.A_t == pytest.approx(A_star, rel=1e-13)
    assert P0 * r.crit_ratio / 1e5 == pytest.approx(12.474666, abs=5e-7)  # throat static < P_m (recovery)
    assert not r.legacy_pressure_check_passed  # 23.61 bar < 37.86 bar threshold
    # default: P_t_u mixes the unscaled flow with the 0.85-scaled Mach number
    d = sns.size_sonic_nozzle_with_feed_line()
    P0_d = d.P_m * (1 + (g - 1) / 2 * d.M_2**2) ** (g / (g - 1))
    assert d.P_t_u / P0_d == pytest.approx(1 / 0.85, rel=1e-13)
    assert P0_d / 1e5 == pytest.approx(22.567686, abs=5e-7)


def test_feedline_warns_when_feed_line_mach_reaches_one():
    with pytest.warns(RuntimeWarning, match="M_2 = .* >= 1"):
        sns.size_sonic_nozzle_with_feed_line(m_dot=20.0)
