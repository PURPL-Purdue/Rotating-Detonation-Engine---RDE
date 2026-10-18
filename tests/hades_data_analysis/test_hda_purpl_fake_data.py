"""Fake ITP / plenum data generators -> HADES.data_analysis.purpl_fake_data."""

import numpy as np
import pytest

from rdelib.matlab_compat import as_rng, matlab_rng, matlab_round, randi
from HADES.data_analysis.purpl_fake_data import (
    PURPL_Fake_ITP_Pressure_Data,
    PURPL_Fake_Plenum_pressure,
    PURPL_test_fake_ITP_data_v2,
    main,
)

from ._octave_harness import ReplayRNG, make_streams, requires_originals, run_octave


def test_matlab_twister_stream():
    # MATLAB: rng(0) (= rng default); rand(1,3) -> 0.8147 0.9058 0.1270 ; rng(1); rand -> 0.4170
    np.testing.assert_allclose(matlab_rng(0).random(3), [0.8147236863931789, 0.9057919370756192, 0.1269868162935061], rtol=1e-15)
    assert matlab_rng(1).random() == pytest.approx(0.417022004702574, rel=1e-14)
    assert isinstance(as_rng(None), np.random.Generator)


def test_matlab_round_and_randi():
    assert matlab_round(2.5) == 3 and matlab_round(-2.5) == -3 and matlab_round(0.49) == 0
    for rng in (np.random.default_rng(0), np.random.RandomState(0)):
        draws = randi(rng, 1, 3, 3000)
        assert set(np.unique(draws)) == {1, 2, 3}


def test_v2_structure_and_csv(tmp_path):
    csv = tmp_path / "synthetic_pressure_data.csv"
    data = PURPL_test_fake_ITP_data_v2(4, csv_path=csv, plot=False)
    assert data.shape == (100000, 3)
    np.testing.assert_allclose(np.diff(data[:, 0]), 1e-5, rtol=1e-9)
    for col in (1, 2):  # scaled exactly into [1, 5] bar
        assert data[:, col].min() == pytest.approx(1.0, abs=1e-14)
        assert data[:, col].max() == pytest.approx(5.0, abs=1e-14)
    back = np.loadtxt(csv, delimiter=",")
    np.testing.assert_allclose(back, data, rtol=1e-14)  # 15 significant digits (writematrix longG)
    # two sensors: same tones, different phases/amplitudes and noise
    assert not np.allclose(data[:, 1], data[:, 2])


def test_v2_is_reproducible_with_seed_and_random_without():
    a = PURPL_test_fake_ITP_data_v2(123, csv_path=None, plot=False)
    b = PURPL_test_fake_ITP_data_v2(123, csv_path=None, plot=False)
    c = PURPL_test_fake_ITP_data_v2(None, csv_path=None, plot=False)
    np.testing.assert_array_equal(a, b)
    assert not np.array_equal(a, c)


def test_v1_delay_and_scaling():
    data = PURPL_Fake_ITP_Pressure_Data(21, plot=False, verbose=False, num_events=0)
    p1, p2 = data[:, 1], data[:, 2]
    # sensor 2 = sensor 1 delayed by 3 samples (noise only, no spikes here)
    lags = np.arange(-6, 7)
    corr = [np.corrcoef(p1[10:-10], np.roll(p2, -lag)[10:-10])[0, 1] for lag in lags]
    assert lags[int(np.argmax(corr))] == 3
    # normalised on sensor 1 -> 30..200 (sensor 1) and 30..218 (sensor 2) +- noise
    assert 27 < p1.min() < 33 and 197 < p1.max() < 203
    assert 215 < p2.max() < 221


def test_v1_spikes_are_synchronised():
    base = PURPL_Fake_ITP_Pressure_Data(ReplayRNG(*make_streams(2)), plot=False, verbose=False, num_events=0)
    spiky = PURPL_Fake_ITP_Pressure_Data(ReplayRNG(*make_streams(2)), plot=False, verbose=False)
    d1 = spiky[:, 1] - base[:, 1]
    d2 = spiky[:, 2] - base[:, 2]
    assert d1.max() >= 1000 - 1e-9
    np.testing.assert_allclose(d2[3:], 0.95 * d1[:-3], atol=1e-9)  # delayed by 3 samples, 95 %


def test_plenum_structure():
    pp = PURPL_Fake_Plenum_pressure(8, plot=False)
    assert pp.shape == (3, 1000)  # MATLAB row orientation [time; h2; air]
    np.testing.assert_allclose(pp[0], np.arange(1000) / 1000, rtol=1e-12)
    assert pp[1:].min() >= 20e5  # abs() rectification: never below 20 bar
    assert pp[1].max() <= 20e5 * (1 + 0.003 + 0.0025 + 0.0015) + 1e-6
    assert pp[2].max() <= 20e5 * (1 + 0.002 + 0.00095 + 0.00045) + 1e-6


def test_main_writes_csv_to_cwd_or_output_dir(no_show, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    main([])
    assert (tmp_path / "synthetic_pressure_data.csv").exists()
    main(["--output-dir", str(tmp_path / "out")])
    assert (tmp_path / "out" / "synthetic_pressure_data.csv").exists()
    (tmp_path / "synthetic_pressure_data.csv").unlink()
    main(["--no-csv"])
    assert not (tmp_path / "synthetic_pressure_data.csv").exists()


@requires_originals("PURPL_Fake_Plenum_pressure.m", "PURPL_test_fake_ITP_data_v2.m")
def test_live_octave_generators(tmp_path):
    streams = make_streams(17)
    out = run_octave(tmp_path / "plenum", "pp = PURPL_Fake_Plenum_pressure();", ["pp"], streams=streams)
    pp = PURPL_Fake_Plenum_pressure(ReplayRNG(*streams), plot=False)
    np.testing.assert_array_equal(out["pp"], pp.ravel(order="F"))  # column-major dump of 3 x N

    out = run_octave(tmp_path / "v2", "d = PURPL_test_fake_ITP_data_v2();", ["d"], streams=streams)
    d = PURPL_test_fake_ITP_data_v2(ReplayRNG(*streams), csv_path=None, plot=False)
    np.testing.assert_allclose(out["d"], d.ravel(order="F"), rtol=1e-14)
