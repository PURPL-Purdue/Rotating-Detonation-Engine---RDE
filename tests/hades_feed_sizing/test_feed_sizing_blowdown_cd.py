"""Regression tests for HADES.feed_sizing.blowdown_orifice_cd.

References (tests/hades_feed_sizing/reference/original_python/):
  blowdown_original_python_*.csv / blowdown_orifice_cd_chat_stdout.txt
      captured from the ORIGINAL ``blowdown_orifice_cd.py`` and
      ``blowdown_orifice_cd_chat.py`` run on the committed
      ``orifice_blowdown_input.csv`` (CoolProp 8.0.0).
  orifice_blowdown_results.csv
      committed output of the earlier (c7402e4) version of
      ``blowdown_orifice_cd.py``; it used psi_to_Pa = 6894.757293168 instead of
      6894.75729, i.e. a 4.6e-10 relative difference in every pressure.
"""

import io
from contextlib import redirect_stdout
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from HADES.feed_sizing import blowdown_orifice_cd as bd  # noqa: E402

REF = Path(__file__).parent / "reference" / "original_python"

# CoolProp-derived values: allow for small differences between CoolProp releases.
COOLPROP_RTOL = 1e-9


@pytest.fixture(scope="module")
def result():
    return bd.blowdown_orifice_cd()


def test_per_sample_values_match_original_scripts(result):
    ref = pd.read_csv(REF / "blowdown_original_python_per_sample.csv", float_precision="round_trip")
    np.testing.assert_array_equal(result.time_ms, ref["time_ms"].to_numpy(float))
    for column in ref.columns[1:]:
        np.testing.assert_allclose(
            getattr(result, column), ref[column].to_numpy(float), rtol=COOLPROP_RTOL, atol=1e-15,
            err_msg=column,
        )


def test_scalars_match_original_scripts(result):
    ref = pd.read_csv(REF / "blowdown_original_python_scalars.csv", float_precision="round_trip")
    for name, value in zip(ref["name"], ref["value"]):
        assert getattr(result, name) == pytest.approx(value, rel=COOLPROP_RTOL), name
    # headline numbers
    assert result.Cd == pytest.approx(0.6402793, abs=1e-7)
    assert result.Cd_integrated == pytest.approx(0.5641651, abs=1e-7)
    assert result.average_Cd == pytest.approx(0.5645317, abs=1e-7)


def test_printed_output_matches_chat_script(result):
    expected = (REF / "blowdown_orifice_cd_chat_stdout.txt").read_text()
    assert result.summary().startswith(expected)


def test_against_committed_blowdown_results_csv(result):
    old = pd.read_csv(REF / "orifice_blowdown_results.csv", float_precision="round_trip")
    np.testing.assert_allclose(result.supply_pressure, old["supply_pressure_abs_pa"], rtol=1e-9)
    np.testing.assert_allclose(result.upstream_pressure, old["upstream_pressure_abs_pa"], rtol=1e-9)
    np.testing.assert_allclose(result.downstream_pressure, old["downstream_pressure_abs_pa"], rtol=1e-9)
    np.testing.assert_allclose(result.tank_masses, old["tank_mass_kg_raw"], rtol=1e-8)
    np.testing.assert_allclose(result.orifice_area, old["orifice_area_m2"], rtol=1e-15)
    # the old script's "full-run blowdown avg mdot" is the run-average m_dot here
    old_avg = (old["tank_mass_kg_used"].iloc[0] - old["tank_mass_kg_used"].iloc[-1]) / (
        old["time_s"].iloc[-1] - old["time_s"].iloc[0]
    )
    assert result.m_dot == pytest.approx(old_avg, rel=1e-8)


def test_independent_hand_checks(result):
    # ideal choked flow, textbook form, first sample
    g, R = result.gammas[0], result.R
    P0, T0 = result.supply_pressure[0], result.supply_temperature_K[0]
    mdot0 = result.orifice_area * P0 * np.sqrt(g / (R * T0)) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)))
    assert result.m_dot_ideals[0] == pytest.approx(mdot0, rel=1e-13)
    # (1500 + 14.6959) psi * 6894.75729
    assert P0 == pytest.approx(1514.6959 * 6894.75729, rel=1e-15)
    # trapezoidal integral and least-squares fit through the origin
    assert result.mass_ideal_total == pytest.approx(np.trapezoid(result.m_dot_ideals, result.time_s), rel=1e-12)
    x = result.mass_ideal_cumulative[1:, None]
    y = result.mass_actual_cumulative[1:]
    assert result.Cd == pytest.approx(np.linalg.lstsq(x, y, rcond=None)[0][0], rel=1e-12)
    assert result.Cd_integrated == pytest.approx(result.mass_actual_total / result.mass_ideal_total)
    assert result.m_dot == pytest.approx(result.mass_actual_total / 1.0, rel=1e-15)
    # dry air: R = 8.314.../0.02896... = 287.05 J/(kg K)
    assert R == pytest.approx(287.05, rel=1e-3)


@pytest.mark.parametrize("bad_temperature", [np.nan, 0.0, -1.0])
def test_invalid_interior_temperature_rejected(bad_temperature):
    """The original scripts raised ValueError from scalar CoolProp; vectorized CoolProp
    returns inf, so the port must reject the sample itself and name the row."""
    data = bd.read_blowdown_csv()
    data.loc[5, "supply_temperature_K"] = bad_temperature
    with pytest.raises(ValueError, match=r"supply_temperature_K must be finite and > 0.*\[5\]"):
        bd.blowdown_orifice_cd(data)


def test_temperature_outside_coolprop_range_rejected():
    data = bd.read_blowdown_csv()
    data.loc[5, "supply_temperature_K"] = 40.0  # below the melting line of air
    with pytest.raises(ValueError, match=r"CoolProp could not evaluate Dmass for Air at .*\[5\].*Tmelt"):
        bd.blowdown_orifice_cd(data)


def test_all_rows_outside_coolprop_range_rejected():
    """If NO row can be evaluated, vectorized PropsSI raises 'No outputs were able to be
    calculated'; the port must still report the failing rows and CoolProp's reason."""
    data = bd.read_blowdown_csv()
    data["supply_temperature_K"] = 40.0
    with pytest.raises(
        ValueError,
        match=r"CoolProp could not evaluate Dmass for Air at data row index\(es\) \[0, 1, .*\(1001 rows\)\].*Tmelt",
    ):
        bd.blowdown_orifice_cd(data)


def test_pressure_below_vacuum_rejected():
    data = bd.read_blowdown_csv()
    data.loc[[7, 9], "supply_pressure_psi"] = -20.0  # psig -> negative absolute
    with pytest.raises(ValueError, match=r"supply_pressure_psi .*\[7, 9\]"):
        bd.blowdown_orifice_cd(data)


def test_non_finite_time_and_bad_scalars_rejected():
    data = bd.read_blowdown_csv()
    bad = data.copy()
    bad.loc[3, "time_ms"] = np.nan
    with pytest.raises(ValueError, match=r"time_ms must be finite.*\[3\]"):
        bd.blowdown_orifice_cd(bad)
    with pytest.raises(ValueError, match="tank_volume_L"):
        bd.blowdown_orifice_cd(data, tank_volume_L=0.0)
    with pytest.raises(ValueError, match="orifice_diameter_in"):
        bd.blowdown_orifice_cd(data, orifice_diameter_in=-0.24)


def test_cli_rejects_missing_temperature_csv(tmp_path, capsys):
    lines = (Path(bd.__file__).parent / "orifice_blowdown_input.csv").read_text().splitlines()
    fields = lines[6].split(",")  # data row index 5
    fields[4] = ""  # blank supply_temperature_K
    lines[6] = ",".join(fields)
    bad_csv = tmp_path / "missing_temperature.csv"
    bad_csv.write_text("\n".join(lines) + "\n")
    with pytest.raises(SystemExit) as exc:
        bd.main(["--input", str(bad_csv), "--no-show"])
    assert exc.value.code == 1
    captured = capsys.readouterr()
    assert "nan" not in captured.out
    assert "supply_temperature_K must be finite and > 0" in captured.err and "[5]" in captured.err
    plt.close("all")


def test_missing_columns_raise(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text("time_ms,supply_pressure_psi\n0,1500\n")
    with pytest.raises(ValueError, match="missing required columns"):
        bd.read_blowdown_csv(bad)


def test_template_has_required_columns():
    header = (Path(bd.__file__).parent / "orifice_blowdown_input_template.csv").read_text().strip()
    assert tuple(header.split(",")) == bd.REQUIRED_COLUMNS


def test_plots_and_cli(result, capsys):
    fig1 = bd.plot_cd_vs_time(result)
    fig2 = bd.plot_mass_calibration(result)
    assert fig1 is not fig2
    plt.close("all")
    bd.main(["--no-show"])
    out = capsys.readouterr().out
    assert out.startswith("Calibrated Cd: 0.6403\n")
    plt.close("all")
