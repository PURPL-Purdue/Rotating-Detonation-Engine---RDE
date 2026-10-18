"""Synthetic (fake) PURPL test data generators.

Consolidates the MATLAB generators from ``HADES/DATA ANALYSIS``:

=====================================  ==========================================
MATLAB file                            Python function (this module)
=====================================  ==========================================
``PURPL_test_fake_ITP_data_v2.m``      :func:`PURPL_test_fake_ITP_data_v2` (canonical)
``PURPL_ITP_data_fft_new.m``           identical copy of v2 -> same function
``PURPL_Fake_ITP_Pressure_Data.m``     :func:`PURPL_Fake_ITP_Pressure_Data`
``PURPL_Fake_Plenum_pressure.m``       :func:`PURPL_Fake_Plenum_pressure`
=====================================  ==========================================

Every generator takes ``rng`` (see ``_matlab_compat.as_rng``): ``None`` for a
fresh random stream (MATLAB behaviour), an ``int`` seed for a reproducible
stream whose ``rand`` values equal MATLAB's ``rng(seed)``, or a numpy
Generator/RandomState.  Random numbers are drawn in exactly the MATLAB order.

Arrays keep the MATLAB orientation: the ITP generators return an ``(N, 3)``
matrix ``[t, p1, p2]`` (columns); the plenum generator returns a ``(3, N)``
matrix ``[time; h2; air]`` (rows).

The v2 generator writes ``synthetic_pressure_data.csv`` to the current
directory, as MATLAB did. Change that with ``csv_path=``, or on the command line::

    uv run python -m HADES.data_analysis.purpl_fake_data [--output-dir DIR | --no-csv]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from ._matlab_compat import as_rng, rand, randi, randn

SYNTHETIC_CSV = "synthetic_pressure_data.csv"


def _figure(num):
    """MATLAB ``figure(num)`` followed by a fresh ``plot`` (replaces old content)."""
    import matplotlib.pyplot as plt

    fig = plt.figure(num)
    fig.clf()
    return fig, fig.add_subplot()


def writematrix(data, filename) -> None:
    """MATLAB ``writematrix`` for a numeric matrix to CSV (15 significant digits, ``longG``)."""
    Path(filename).parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(filename, np.asarray(data), fmt="%.15g", delimiter=",")


# ----------------------------------------------------------------------
#  ITP pressure data, v2 (Sept 2026, canonical)
# ----------------------------------------------------------------------
def PURPL_test_fake_ITP_data_v2(
    rng=None,
    *,
    csv_path: str | Path | None = SYNTHETIC_CSV,
    plot: bool = True,
    fs: float = 100000,
    duration: float = 1,
    f_min: float = 10000,
    f_max: float = 25000,
) -> np.ndarray:
    """Two ITP pressure traces sharing 3 random tones in 10-25 kHz, scaled to 1-5 bar.

    Returns ``data = [t, pressure1, pressure2]`` (N x 3) and, like MATLAB, writes
    it to ``synthetic_pressure_data.csv`` in the current directory
    (``csv_path=None`` skips the file).
    """
    rng = as_rng(rng)

    # Parameters
    n_points = int(round(fs * duration))  # 100,000 total points
    t = np.arange(n_points) / fs  # Time vector

    # 1. Pick 3 random dominant frequencies between f_min and f_max
    freqs = f_min + (f_max - f_min) * rand(rng, 3)
    # Random phases (0 to 2*pi)
    phases1 = 2 * np.pi * rand(rng, 3)
    phases2 = 2 * np.pi * rand(rng, 3)
    # Random relative amplitudes for the 3 tones
    amp1 = 0.2 + 0.8 * rand(rng, 3)
    amp2 = 0.2 + 0.8 * rand(rng, 3)

    # 2. Combine high-frequency tones with noise
    raw_signal1 = (
        amp1[0] * np.sin(2 * np.pi * freqs[0] * t + phases1[0])
        + amp1[1] * np.sin(2 * np.pi * freqs[1] * t + phases1[1])
        + amp1[2] * np.sin(2 * np.pi * freqs[2] * t + phases1[2])
        + 0.1 * randn(rng, n_points)
    )
    raw_signal2 = (
        amp2[0] * np.sin(2 * np.pi * freqs[0] * t + phases2[0])
        + amp2[1] * np.sin(2 * np.pi * freqs[1] * t + phases2[1])
        + amp2[2] * np.sin(2 * np.pi * freqs[2] * t + phases2[2])
        + 0.1 * randn(rng, n_points)
    )

    # 3. Scale signals into the [1, 5] bar range
    def scale_to_range(x, target_min, target_max):
        return target_min + (x - x.min()) * (target_max - target_min) / (x.max() - x.min())

    pressure1 = scale_to_range(raw_signal1, 1.0, 5.0)
    pressure2 = scale_to_range(raw_signal2, 1.0, 5.0)

    # 4. Save to CSV
    data = np.column_stack([t, pressure1, pressure2])
    if csv_path is not None:
        writematrix(data, csv_path)

    # Plot pressures over time (confirm randomness)
    if plot:
        fig, ax = _figure(1)
        ax.plot(data[:, 0], data[:, 1], label="Sensor 1")
        ax.plot(data[:, 0], data[:, 2], label="Sensor 2")
        ax.grid(True)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Pressure (bar)")
        ax.set_title("Fake ITP Pressure Data")
        ax.legend(loc="best")

    return data


# ----------------------------------------------------------------------
#  ITP pressure data, v1 (April 2026): delayed sensor pair with spikes
# ----------------------------------------------------------------------
def PURPL_Fake_ITP_Pressure_Data(
    rng=None,
    *,
    plot: bool = True,
    verbose: bool = True,
    n_points: int = 100000,
    fs: float = 100000,
    delay_samples: int = 3,
    num_events: int = 150,
) -> np.ndarray:
    """High-magnitude, equal-amplitude correlated pressure data with a time delay.

    Three equal tones (1-2 kHz, 2-4 kHz, 4-8 kHz); sensor 2 is the same signal
    delayed by ``delay_samples`` samples; both get Gaussian noise and 150
    synchronized (delayed) exponential spikes.  Returns ``[t, p1, p2]`` (N x 3).
    """
    rng = as_rng(rng)
    dt = 1 / fs  # 100 kHz (10 microseconds per sample)
    t = np.arange(n_points) * dt

    # 1. Define equal-magnitude dominant frequencies
    f_dom = 1000 + (2000 - 1000) * rand(rng)
    f2 = 2000 + (4000 - 2000) * rand(rng)
    f3 = 4000 + (8000 - 4000) * rand(rng)

    # The "source" signal
    shared_signal = 300 * (np.sin(2 * np.pi * f_dom * t) + np.sin(2 * np.pi * f2 * t) + np.sin(2 * np.pi * f3 * t))

    # 2. Apply time delay to the second dataset (zero-padded at the start, same length)
    phys_signal1 = shared_signal
    phys_signal2 = np.concatenate([np.zeros(delay_samples), shared_signal[: n_points - delay_samples]])

    # 3. Normalize (based on sensor 1) and scale
    s_min = phys_signal1.min()
    s_max = phys_signal1.max()
    norm1 = (phys_signal1 - s_min) / (s_max - s_min)
    norm2 = (phys_signal2 - s_min) / (s_max - s_min)

    pressure1 = norm1 * 170 + 30 + 0.5 * randn(rng, n_points)
    pressure2 = norm2 * 188 + 30 + 0.5 * randn(rng, n_points)  # Slightly different scaling

    # 4. Add synchronized spikes (also delayed).  idx is MATLAB 1-based.
    for _ in range(num_events):
        idx = randi(rng, 1, n_points - 10)
        width = randi(rng, 1, 3)
        mag = 1000
        pulse = mag * np.exp(-3 * np.arange(width))
        L = len(pulse)

        # Sensor 1 gets it at idx
        pressure1[idx - 1 : idx - 1 + L] += pulse

        # Sensor 2 gets it at idx + delay
        idx2 = idx + delay_samples
        if idx2 + L - 1 <= n_points:
            pressure2[idx2 - 1 : idx2 - 1 + L] += 0.95 * pulse

    fake_pressure_data = np.column_stack([t, pressure1, pressure2])

    if plot:
        fig, ax = _figure(2)
        ax.plot(t, pressure1)
        ax.plot(t, pressure2)
        ax.set_xlabel("Time (t)")
        ax.set_ylabel("Pressure (kPa)")
        ax.legend(["ITP 1", "ITP 2"], loc="best")
        ax.set_title("Fake ITP Data")

    if verbose:
        print("Dominant Frequency: %.5f (kHz)" % (f_dom / 1000))
        print("Prominent Frequency 2: %.5f (kHz)" % (f2 / 1000))
        print("Prominent Frequency 3: %.5f (kHz)" % (f3 / 1000))

    return fake_pressure_data


# ----------------------------------------------------------------------
#  Plenum static pressure (April 2026)
# ----------------------------------------------------------------------
def PURPL_Fake_Plenum_pressure(rng=None, *, plot: bool = True) -> np.ndarray:
    """One second of fake H2/air plenum static pressure around 20 bar (1 kHz sampling).

    Returns ``[time; plenum_pressure_h2; plenum_pressure_air]`` as a 3 x N array
    (row 2 = HYDROGEN, row 3 = AIR), in Pa.

    Faithful quirks of the MATLAB source (kept, see PORTING_NOTES.md): the
    oscillations use ``sin(time .* f)`` (``f`` acts as rad/s, not Hz), the
    third H2 term reuses ``f2`` so ``f3`` is drawn but unused, and the summed
    oscillation is rectified with ``abs``.
    """
    rng = as_rng(rng)

    base_plenum_pressure = 20  # 20 bar expected static pressure
    base_plenum_pa = base_plenum_pressure * 100000  # base static pressure
    f_samp = 1000  # sampling frequency
    dt = 1 / f_samp  # time step
    time = np.arange(f_samp) * dt  # one second of samples

    # Random frequencies causing minor oscillations in the plenum pressures
    f1 = 300 + (600 - 300) * rand(rng)
    f2 = 600 + (1200 - 600) * rand(rng)
    f3 = 1000 + (1300 - 1000) * rand(rng)  # noqa: F841  (drawn but never used, as in MATLAB)

    # Baseline pressure for the hydrogen and air plenums
    pressure_h2 = np.ones(len(time)) * base_plenum_pa
    pressure_air = np.ones(len(time)) * base_plenum_pa

    # Oscillations scaled for each frequency
    plenum_pressure_h2 = pressure_h2 + np.abs(
        (0.003 * pressure_h2 * np.sin(time * f1))
        + (0.0025 * pressure_h2 * np.sin(time * f2))
        + (0.0015 * pressure_h2 * np.sin(time * f2))
    )
    plenum_pressure_air = pressure_air + np.abs(
        (0.002 * pressure_air * np.sin(time * f1 * 2))
        + (0.00095 * pressure_air * np.sin(time * f2 * 3))
        + (0.00045 * pressure_air * np.sin(time * f2 * 1.5))
    )

    fake_plenum_pressure = np.vstack([time, plenum_pressure_h2, plenum_pressure_air])

    if plot:
        fig, ax = _figure(1)
        ax.plot(time, plenum_pressure_h2 / 100000)
        ax.plot(time, plenum_pressure_air / 100000)
        ax.set_ylabel("Pressure (bar)")
        ax.set_xlabel("Time (s)")
        ax.set_title("Static Plenum Pressure Over Time")
        ax.legend(["Hydrogen Plenum", "Air Plenum"])
        ax.grid(True)

    return fake_plenum_pressure


def add_csv_arguments(parser: argparse.ArgumentParser) -> None:
    """CLI options for the CSV written by the v2 generator (default: ./synthetic_pressure_data.csv)."""
    g = parser.add_mutually_exclusive_group()
    g.add_argument("--output-dir", type=Path, default=Path.cwd(), help="folder for synthetic_pressure_data.csv (default: cwd)")
    g.add_argument("--no-csv", action="store_true", help="do not write synthetic_pressure_data.csv")


def csv_path_from_args(args) -> Path | None:
    return None if args.no_csv else Path(args.output_dir) / SYNTHETIC_CSV


def main(argv: list[str] | None = None) -> None:
    import matplotlib.pyplot as plt

    p = argparse.ArgumentParser(prog="python -m HADES.data_analysis.purpl_fake_data")
    add_csv_arguments(p)
    args = p.parse_args(argv)
    PURPL_test_fake_ITP_data_v2(csv_path=csv_path_from_args(args))
    plt.show()


if __name__ == "__main__":
    main()
