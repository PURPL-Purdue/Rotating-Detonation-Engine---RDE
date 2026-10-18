"""Shared fixtures for the HADES sizing-model tests."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

REF = Path(__file__).with_name("reference")


@pytest.fixture(scope="session")
def ref_dir() -> Path:
    return REF


@pytest.fixture(scope="session")
def octave_main() -> dict:
    """Original HADES_size_main_v1.m run under Octave 11.3 (see reference/octave_harness)."""
    return json.loads((REF / "octave_main.json").read_text())


@pytest.fixture(scope="session")
def octave_units() -> dict:
    return json.loads((REF / "octave_units.json").read_text())


@pytest.fixture(scope="session")
def octave_cool() -> dict:
    return json.loads((REF / "octave_cool.json").read_text())


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")
