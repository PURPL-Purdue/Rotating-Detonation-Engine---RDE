"""PURPL_test_FFT_v2.m / PURPL_FFT_new.m / ITP_FFT_for_PURPL_FakeData.m -> purpl_fft.

Golden numbers come from running the ORIGINAL .m files in GNU Octave 11.3 with
``rand``/``randn``/``randi`` replaced by replayed streams
(``_octave_harness.make_streams(seed)``); the Python port consumes the same
streams through ``ReplayRNG``.
"""

import numpy as np
import pytest

from HADES.data_analysis._matlab_compat import matlab_rng
from HADES.data_analysis.purpl_fake_data import (
    PURPL_Fake_ITP_Pressure_Data,
    PURPL_test_fake_ITP_data_v2,
)
from HADES.data_analysis.purpl_fft import (
    ITP_FFT_for_PURPL_FakeData,
    PURPL_test_FFT_v2,
    dominant_frequencies,
    main,
    single_sided_amplitude_spectrum,
)

from ._octave_harness import ReplayRNG, make_streams, octave_available, run_octave


def test_single_sided_amplitude_convention():
    fs, n = 100000, 100000
    t = np.arange(n) / fs
    x = 2.0 + 3.0 * np.sin(2 * np.pi * 1000 * t) + 0.5 * np.cos(2 * np.pi * 12000 * t)
    f, P1 = single_sided_amplitude_spectrum(x, fs)
    assert len(f) == n // 2 + 1 and f[1] == pytest.approx(1.0) and f[-1] == pytest.approx(fs / 2)
    assert P1[0] == pytest.approx(2.0, rel=1e-12)  # DC not doubled
    assert P1[1000] == pytest.approx(3.0, rel=1e-10)  # tone amplitude recovered
    assert P1[12000] == pytest.approx(0.5, rel=1e-10)
    assert np.sum(P1 > 1e-9) == 3


def test_odd_length_last_bin_is_doubled():
    # 9 samples of 3*cos(2*pi*4*n/9): bin 4 is an ordinary positive frequency (no Nyquist bin).
    n = np.arange(9)
    f, P1 = single_sided_amplitude_spectrum(3 * np.cos(2 * np.pi * 4 * n / 9), 9.0)
    assert len(P1) == 5 and f[-1] == pytest.approx(4.0)
    assert P1[4] == pytest.approx(3.0, rel=1e-12)  # MATLAB scripts (2:end-1 doubling) gave 1.5
    assert np.all(P1[:4] < 1e-12)
    # odd N with interior tone and DC
    n = np.arange(101)
    f, P1 = single_sided_amplitude_spectrum(0.5 + 2 * np.sin(2 * np.pi * 50 * n / 101), 101.0)
    assert P1[0] == pytest.approx(0.5, rel=1e-12) and P1[50] == pytest.approx(2.0, rel=1e-12)


def test_matlab_doubling_on_odd_length_documented():
    # The inherited MATLAB rule (double P1(2:end-1) for every N), reproduced to pin the bug.
    n = np.arange(9)
    x = 3 * np.cos(2 * np.pi * 4 * n / 9)
    p = np.abs(np.fft.fft(x) / 9)[:5]
    p[1:-1] *= 2
    assert p[4] == pytest.approx(1.5, rel=1e-12)  # Octave on the original lines: 1.5


def test_even_length_nyquist_not_doubled_and_unchanged():
    # Even N: Nyquist bin holds the full amplitude of a cos at fs/2 and must not be doubled.
    n = np.arange(8)
    f, P1 = single_sided_amplitude_spectrum(1.0 + 3 * np.cos(np.pi * n), 8.0)
    assert f[-1] == pytest.approx(4.0)
    assert P1[-1] == pytest.approx(3.0, rel=1e-12) and P1[0] == pytest.approx(1.0, rel=1e-12)
    # identical to the MATLAB rule for any even N (all fake ITP data are even length)
    x = np.random.default_rng(0).standard_normal(1000)
    p = np.abs(np.fft.fft(x) / 1000)[:501]
    p[1:-1] *= 2
    np.testing.assert_array_equal(single_sided_amplitude_spectrum(x, 1.0)[1], p)


def test_numpy_fft_equals_matlab_dft_definition():
    # MATLAB fft: X(k) = sum_n x(n) exp(-2*pi*i*(k-1)(n-1)/N), no normalisation
    x = np.array([1.0, 2.0, -1.0, 0.5, 3.0])
    n = np.arange(5)
    X = np.array([np.sum(x * np.exp(-2j * np.pi * k * n / 5)) for k in range(5)])
    np.testing.assert_allclose(np.fft.fft(x), X, atol=1e-12)


@pytest.mark.parametrize("seed", [1, 7, 2026])
def test_v2_recovers_injected_frequencies(seed):
    data = PURPL_test_fake_ITP_data_v2(seed, csv_path=None, plot=False)
    injected = 10000 + 15000 * matlab_rng(seed).random(3)  # first three rand draws
    f, A = PURPL_test_FFT_v2(data, plot=False)
    assert A[0] == 0  # DC removed
    found = np.sort(dominant_frequencies(f, A, n=3, min_separation=20))
    np.testing.assert_allclose(found, np.sort(injected), atol=0.5 + 1e-9)  # 1 Hz bins


@pytest.mark.parametrize("seed", [3, 11])
def test_v1_raw_pipeline_recovers_frequencies_despite_delay_and_spikes(seed):
    data = PURPL_Fake_ITP_Pressure_Data(seed, plot=False, verbose=False)
    u = matlab_rng(seed).random(3)
    injected = np.array([1000 + 1000 * u[0], 2000 + 2000 * u[1], 4000 + 4000 * u[2]]) / 1000
    f_kHz, M = ITP_FFT_for_PURPL_FakeData(data, plot=False)
    assert len(f_kHz) == len(M) == 50000
    found = np.sort(dominant_frequencies(f_kHz, M, n=3, min_separation=0.02))
    np.testing.assert_allclose(found, np.sort(injected), atol=0.5e-3 + 1e-12)


def test_v1_axis_fix_for_records_not_one_second_long():
    # MATLAB used linspace(1,fs,fs+1)/1000, only valid for N == fs.  With a
    # 0.25 s record the port's k*fs/N axis still puts the tone where it belongs.
    fs, n = 100000, 25000
    t = np.arange(n) / fs
    p = 3 + np.sin(2 * np.pi * 4000 * t)
    f_kHz, M = ITP_FFT_for_PURPL_FakeData(np.column_stack([t, p, p]), plot=False)
    assert f_kHz[np.argmax(M)] == pytest.approx(4.0, abs=1e-9)
    assert f_kHz[0] == pytest.approx(fs / n / 1000)
    # MATLAB's axis would have labelled that bin 1 + 999*(fs-1)/fs Hz -> wrong by 3 kHz
    matlab_f = np.linspace(1, fs, fs + 1) / 1000
    assert abs(matlab_f[np.argmax(M)] - 4.0) > 2.9


def test_v2_matches_octave_golden():
    data = PURPL_test_fake_ITP_data_v2(ReplayRNG(*make_streams(7)), csv_path=None, plot=False)
    f, A = PURPL_test_FFT_v2(data, plot=False)
    assert data[12345, 1] == pytest.approx(1.7571135094285344, rel=1e-14)
    assert A.sum() == pytest.approx(35.85470541655545, rel=1e-12)
    top = np.argsort(A)[::-1][:3]
    np.testing.assert_array_equal(top, [16576, 11145, 21699])
    np.testing.assert_allclose(f[top], [16576, 11145, 21699], rtol=1e-12)
    np.testing.assert_allclose(A[top], [0.6202950300386583, 0.4598130340225889, 0.44446194776439923], rtol=1e-12)


def test_v1_matches_octave_golden(capsys):
    data = PURPL_Fake_ITP_Pressure_Data(ReplayRNG(*make_streams(11)), plot=False)
    out = capsys.readouterr().out
    assert out == (
        "Dominant Frequency: 1.18027 (kHz)\n"
        "Prominent Frequency 2: 2.03895 (kHz)\n"
        "Prominent Frequency 3: 5.85287 (kHz)\n"
    )
    f, M = ITP_FFT_for_PURPL_FakeData(data, plot=False)
    assert data[54321, 1] == pytest.approx(146.16569258298426, rel=1e-14)
    assert M.sum() == pytest.approx(537598730.7857089, rel=1e-12)
    top = np.argsort(M)[::-1][:3]
    np.testing.assert_array_equal(top, [2038, 5852, 1179])
    np.testing.assert_allclose(f[top], [2.039, 5.853, 1.18], rtol=1e-12)
    np.testing.assert_allclose(M[top], [1479080.6797434287, 1438669.2031335854, 1328998.1114200251], rtol=1e-12)


def test_default_pipelines_and_main_smoke(no_show, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)  # v2 generator writes synthetic_pressure_data.csv to cwd (as MATLAB did)
    f, A = PURPL_test_FFT_v2()
    assert (tmp_path / "synthetic_pressure_data.csv").exists() and len(f) == 50001
    f_kHz, M = ITP_FFT_for_PURPL_FakeData()
    assert len(f_kHz) == 50000
    main([])
    main(["--raw"])


@pytest.mark.skipif(not octave_available("PURPL_test_FFT_v2.m", "PURPL_test_fake_ITP_data_v2.m"),
                    reason="Octave or original .m not available")  # fmt: skip
def test_live_octave_v2(tmp_path):
    streams = make_streams(5)
    out = run_octave(tmp_path, "PURPL_test_FFT_v2();", ["plot1x", "plot1y", "plot2y", "plot3x", "plot3y"], streams=streams)
    data = PURPL_test_fake_ITP_data_v2(ReplayRNG(*streams), csv_path=None, plot=False)
    f, A = PURPL_test_FFT_v2(data, plot=False)
    np.testing.assert_array_equal(out["plot1x"], data[:, 0])
    np.testing.assert_allclose(out["plot1y"], data[:, 1], rtol=1e-14)
    np.testing.assert_allclose(out["plot2y"], data[:, 2], rtol=1e-14)
    np.testing.assert_allclose(out["plot3x"], f, rtol=1e-14)
    np.testing.assert_allclose(out["plot3y"], A, rtol=0, atol=1e-14)
    # the MATLAB generator's CSV (writematrix) round-trips to the same matrix
    csv = np.loadtxt(tmp_path / "synthetic_pressure_data.csv", delimiter=",")
    np.testing.assert_allclose(csv, data, rtol=1e-14)


@pytest.mark.skipif(not octave_available("ITP_FFT_for_PURPL_FakeData.m", "PURPL_Fake_ITP_Pressure_Data.m"),
                    reason="Octave or original .m not available")  # fmt: skip
def test_live_octave_v1(tmp_path):
    streams = make_streams(9)
    out = run_octave(tmp_path, "ITP_FFT_for_PURPL_FakeData();", ["plot1y", "plot2y", "plot3x", "plot3y"], streams=streams)
    data = PURPL_Fake_ITP_Pressure_Data(ReplayRNG(*streams), plot=False, verbose=False)
    f_kHz, M = ITP_FFT_for_PURPL_FakeData(data, plot=False)
    np.testing.assert_allclose(out["plot1y"], data[:, 1], rtol=1e-14)
    np.testing.assert_allclose(out["plot2y"], data[:, 2], rtol=1e-14)
    np.testing.assert_allclose(out["plot3x"], f_kHz, rtol=1e-12)  # MATLAB axis == k*fs/N for N == fs
    np.testing.assert_allclose(out["plot3y"], M, rtol=1e-12, atol=1e-6)
