"""Shared fixtures for the HADES data-analysis tests (non-interactive plotting)."""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


@pytest.fixture
def no_show(monkeypatch):
    """Make ``plt.show()`` a no-op so module ``main()`` functions can run."""
    monkeypatch.setattr(plt, "show", lambda *a, **k: None)
