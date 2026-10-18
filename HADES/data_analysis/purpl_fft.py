"""FFT / spectrum analysis of the two ITP pressure transducers.

Consolidates the MATLAB FFT scripts from ``HADES/DATA ANALYSIS``:

=====================================  ==========================================
MATLAB file                            Python function (this module)
=====================================  ==========================================
``PURPL_test_FFT_v2.m``                :func:`PURPL_test_FFT_v2` (canonical pipeline)
``PURPL_FFT_new.m``                    identical to v2 (MATLAB copy called a
                                       non-existent ``PURPL_ITP_data_fft``) -> v2
``ITP_FFT_for_PURPL_FakeData.m``       :func:`ITP_FFT_for_PURPL_FakeData`
=====================================  ==========================================

Spectrum conventions (MATLAB ``fft`` and ``numpy.fft.fft`` are the same
unnormalised DFT, so no rescaling is needed):

* :func:`single_sided_amplitude_spectrum` is MATLAB's documented one-sided
  amplitude spectrum: ``P2 = |fft(x)/N|``, keep bins ``0..floor(N/2)``, and
  double every bin except DC and, for even N, Nyquist. A pure tone
  ``A*sin(2*pi*f*t)`` that falls on a bin shows up with height ``A``. The MATLAB
  scripts also skipped doubling the last bin for odd N; that is fixed here, and
  even-N results (all the fake data) are unchanged.
* :func:`ITP_FFT_for_PURPL_FakeData` keeps the older raw ``|fft(x)|``
  (unnormalised) magnitude, DC removed, plotted against kHz.

Both pipelines average the spectra of the two sensors.  These functions accept
real test data too: pass any ``(N, 3)`` ``[t, p1, p2]`` matrix as ``data``.

Run as a script::

    uv run python -m HADES.data_analysis.purpl_fft            # v2 pipeline [--output-dir DIR | --no-csv]
    uv run python -m HADES.data_analysis.purpl_fft --raw      # April raw-|FFT| pipeline
"""

from __future__ import annotations

import argparse

import numpy as np

from .purpl_fake_data import (
    SYNTHETIC_CSV,
    PURPL_Fake_ITP_Pressure_Data,
    PURPL_test_fake_ITP_data_v2,
    _figure,
    add_csv_arguments,
    csv_path_from_args,
)


def single_sided_amplitude_spectrum(x, fs: float):
    """One-sided amplitude spectrum of ``x`` sampled at ``fs``.

    Returns ``(f, P1)`` with ``f = (0:floor(N/2)) * fs/N`` (Hz). DC is never
    doubled. The Nyquist bin is not doubled either, and it exists only when N
    is even.

    FIX: the MATLAB scripts doubled ``P1(2:end-1)`` for every N. For odd N the
    last kept bin is an ordinary positive frequency, so its amplitude came out
    halved. The fix only changes odd-length records.
    """
    x = np.asarray(x, dtype=float)
    num_samples = len(x)
    p_two_sided = np.abs(np.fft.fft(x) / num_samples)
    half_N = num_samples // 2
    p_single_sided = p_two_sided[: half_N + 1].copy()
    if num_samples % 2 == 0:
        p_single_sided[1:-1] = 2 * p_single_sided[1:-1]  # keep DC and Nyquist (MATLAB 2:end-1)
    else:
        p_single_sided[1:] = 2 * p_single_sided[1:]  # odd N: no Nyquist bin
    f = np.arange(half_N + 1) * (fs / num_samples)
    return f, p_single_sided


def dominant_frequencies(f, spectrum, n: int = 3, min_separation: float = 0.0):
    """The ``n`` strongest local peaks of ``spectrum`` (Hz), strongest first.

    Peaks closer than ``min_separation`` (Hz) to a stronger peak are skipped
    (useful to ignore leakage skirts).  Not in the MATLAB scripts, which only
    plotted the spectrum; provided to read the detonation frequency off real
    data and for regression tests.
    """
    from scipy.signal import find_peaks

    f = np.asarray(f, dtype=float)
    spectrum = np.asarray(spectrum, dtype=float)
    peaks, _ = find_peaks(spectrum)
    peaks = peaks[np.argsort(spectrum[peaks])[::-1]]
    chosen: list[float] = []
    for p in peaks:
        if all(abs(f[p] - c) >= min_separation for c in chosen):
            chosen.append(f[p])
        if len(chosen) == n:
            break
    return np.array(chosen)


def PURPL_test_FFT_v2(data=None, *, plot: bool = True, rng=None, csv_path=SYNTHETIC_CSV):
    """Averaged one-sided amplitude spectrum of the two ITP sensors (canonical pipeline).

    ``data`` is an ``(N, 3)`` ``[t, p1, p2]`` matrix; by default the v2 fake data
    are generated (as in MATLAB, which also writes ``synthetic_pressure_data.csv``
    and plots the traces in figure 1).  The DC bin is zeroed ("always super high
    at 0 frequency" -- it is the mean pressure).

    Returns ``(all_f, average_fft)``; the MATLAB function returned nothing and
    only plotted (figure 2).
    """
    if data is None:
        data = PURPL_test_fake_ITP_data_v2(rng, csv_path=csv_path, plot=plot)
    data = np.asarray(data, dtype=float)

    time = data[:, 0]
    pressure_1 = data[:, 1]
    pressure_2 = data[:, 2]

    time_step = time[1] - time[0]  # time step between samples (1/100000 s)
    sampling_freq = 1 / time_step  # sampling frequency (should be 100 kHz)

    all_f, p1_single_sided = single_sided_amplitude_spectrum(pressure_1, sampling_freq)
    _, p2_single_sided = single_sided_amplitude_spectrum(pressure_2, sampling_freq)

    # average the two spectra and delete the DC point
    average_fft = (p1_single_sided + p2_single_sided) / 2
    average_fft[0] = 0

    if plot:
        fig, ax = _figure(2)
        ax.plot(all_f, average_fft, linewidth=1.25)
        ax.grid(True)
        ax.set_xlabel("Frequencies (Hz)")
        ax.set_ylabel("Relative Strength")
        ax.set_title("FFT Plot")

    return all_f, average_fft


def ITP_FFT_for_PURPL_FakeData(data=None, *, plot: bool = True, rng=None):
    """April-2026 pipeline: averaged raw ``|fft|`` magnitude of the two ITP sensors.

    ``data`` defaults to :func:`PURPL_Fake_ITP_Pressure_Data`.  Bin 0 (DC) is
    dropped (MATLAB's "shift the FFT back one entry"), so element ``k-1`` of the
    result is DFT bin ``k`` at ``k*fs/N``, ``k = 1..floor(N/2)``.

    Returns ``(f_kHz, magnitude)``.  FIX: MATLAB built the axis as
    ``linspace(1, fs, fs+1)/1000``, which is only right for exactly one second
    of data (N == fs); the port uses ``k*fs/N``.  For the 1 s fake data both
    agree to ~1e-11 Hz (see PORTING_NOTES.md).
    """
    if data is None:
        data = PURPL_Fake_ITP_Pressure_Data(rng, plot=plot)
    data = np.asarray(data, dtype=float)

    time = data[:, 0]
    itp_pt_1 = data[:, 1]
    itp_pt_2 = data[:, 2]

    dt = time[1] - time[0]  # time step
    fs = 1 / dt  # sampling frequency
    n = len(time)  # number of samples

    # magnitude of each FFT, then average
    fft1 = np.abs(np.fft.fft(itp_pt_1))
    fft2 = np.abs(np.fft.fft(itp_pt_2))
    combined_fft = (fft1 + fft2) / 2

    # drop DC: fft_final[k-1] is bin k
    fft_final = combined_fft[1:]
    half = n // 2
    k = np.arange(1, half + 1)
    f = k * fs / n / 1000  # kHz

    if plot:
        fig, ax = _figure(1)
        ax.plot(f, fft_final[:half])
        ax.set_title("Fast Fourier Transform of Data set")
        ax.set_ylabel(r"$|a_{k}|$")
        ax.set_xlabel("Frequency (kHz)")

    return f, fft_final[:half]


def main(argv: list[str] | None = None) -> None:
    import matplotlib.pyplot as plt

    p = argparse.ArgumentParser(prog="python -m HADES.data_analysis.purpl_fft")
    p.add_argument("--raw", action="store_true", help="April raw-|FFT| pipeline (no CSV is written)")
    add_csv_arguments(p)
    args = p.parse_args(argv)
    if args.raw:
        ITP_FFT_for_PURPL_FakeData()
    else:
        PURPL_test_FFT_v2(csv_path=csv_path_from_args(args))
    plt.show()


if __name__ == "__main__":
    main()
