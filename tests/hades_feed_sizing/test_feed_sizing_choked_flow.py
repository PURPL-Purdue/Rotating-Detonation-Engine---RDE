"""Hand-calculation checks of the local choked-flow helpers."""

import math

import numpy as np
import pytest

from HADES.feed_sizing import choked_flow as cf


def test_critical_pressure_ratio_air():
    # (2/2.4)^3.5 = 0.528281787717...
    assert cf.critical_pressure_ratio(1.4) == pytest.approx(0.5282817877171742, rel=1e-15)
    assert cf.critical_pressure_ratio(1.31) == pytest.approx((2 / 2.31) ** (1.31 / 0.31), rel=1e-15)


def test_choked_flow_parameter_air():
    # sqrt(1.4/287.05) * (1/1.2)^3 = 0.0404149...
    assert cf.choked_flow_parameter(1.4, 287.05) == pytest.approx(
        math.sqrt(1.4 / 287.05) / 1.2**3, rel=1e-14
    )


def test_choked_relations_round_trip():
    K = cf.choked_flow_parameter(1.4, 287.05)
    A = cf.choked_area(1.0, 0.85, 40e5, 283.0, K)
    assert cf.choked_mass_flow(0.85, A, 40e5, 283.0, 1.4, 287.05) == pytest.approx(1.0, rel=1e-14)
    assert cf.choked_upstream_pressure(np.array([1.0, 2.0]), 0.85, A, 283.0, K) == pytest.approx(
        [40e5, 80e5], rel=1e-14
    )
    # textbook: mdot = A * P0 * sqrt(gamma/(R T0)) * (2/(gamma+1))^((gamma+1)/(2(gamma-1)))
    textbook = 1e-4 * 50e5 * math.sqrt(1.4 / (287.05 * 300.0)) * (2 / 2.4) ** 3
    assert cf.choked_mass_flow(1.0, 1e-4, 50e5, 300.0, 1.4, 287.05) == pytest.approx(textbook, rel=1e-14)


def test_area_to_diameter():
    assert cf.area_to_diameter(math.pi / 4) == pytest.approx(1.0, rel=1e-15)
