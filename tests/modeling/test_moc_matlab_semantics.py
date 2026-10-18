"""MATLAB-semantics helpers used by the MoC port (MODELING/moc/_matlab.py)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from MODELING.moc._matlab import (
    GrowArray,
    GrowColumns,
    GrowVector,
    asind,
    atand,
    colon,
    complex_ordering,
    cosd,
    fzero,
    mldivide2,
    mle,
    mlt,
    mmax,
    msqrt,
    narrow,
    num2str,
    sind,
    tand,
)


def test_colon_matches_matlab_construction():
    v = colon(0, 0.1, 1)
    assert v.size == 11 and v[0] == 0.0 and v[-1] == 1.0
    assert v[3] == 3 * 0.1  # built from the start ...
    assert v[7] == 1 - 3 * 0.1  # ... and symmetrically from the end (MATLAB), not 7*0.1
    assert v[5] == 0.5
    np.testing.assert_array_equal(colon(1, 1, 15), np.arange(1, 16))
    np.testing.assert_array_equal(colon(1, 5), np.arange(1, 6))
    assert colon(1, 1, 0).size == 0
    assert colon(1.1, 0.039, 5.0)[-1] == 5.0  # oblique_shock Mach grid lands on the end point


def test_degree_trig_exact_values():
    assert sind(180.0) == 0.0 and sind(360.0) == 0.0 and sind(30.0) == pytest.approx(0.5, abs=1e-16)
    assert cosd(90.0) == 0.0 and cosd(270.0) == 0.0
    assert tand(180.0) == 0.0 and tand(90.0) == math.inf
    assert tand(45.0) == pytest.approx(1.0, abs=2e-16)
    np.testing.assert_allclose(sind(np.array([0.0, 90.0, 180.0])), [0.0, 1.0, 0.0], atol=0)


def test_asind_and_sqrt_go_complex_like_matlab():
    assert asind(1.0) == 90.0
    z = asind(2.0)  # MATLAB: asind(2) = 90 - 75.4561i
    assert z.real == pytest.approx(90.0) and z.imag == pytest.approx(-75.45612929021272)
    assert msqrt(-4.0) == 2j
    assert isinstance(msqrt(4.0), np.float64)


# Reference values printed by GNU Octave 11.3.0 (format %.16g) for the same expressions.
OCTAVE_COMPLEX_TRIG = [
    ("asind([2 .5 -3])", lambda: asind(np.array([2.0, 0.5, -3.0])),
     [89.99999999999997 - 75.45612929021688j, 30 + 0j, -89.99999999999994 + 100.9979734210525j]),
    ("asind(2)", lambda: asind(2.0), [89.99999999999997 - 75.45612929021688j]),
    ("asind(-3)", lambda: asind(-3.0), [-89.99999999999994 + 100.9979734210525j]),
    ("asind(.5+.5i)", lambda: asind(0.5 + 0.5j), [25.91364618649388 + 30.40329097482185j]),
    ("sind(1i)", lambda: sind(1j), [0.01745417862959511j]),
    ("sind([180+1i,30,-180+2i])", lambda: sind(np.array([180 + 1j, 30, -180 + 2j])),
     [-0.01745417862959511j, 0.4999999999999999 + 0j, -0.03491367424101725j]),
    ("cosd(-30+1i)", lambda: cosd(-30 + 1j), [0.8661573103447925 + 0.008727089314797556j]),
    ("cosd(90+2i)", lambda: cosd(90 + 2j), [-0.03491367424101725j]),
    ("cosd([-90+1i,.5i])", lambda: cosd(np.array([-90 + 1j, 0.5j])),
     [0.01745417862959511j, 1.000038077419119 + 5.343597527217299e-19j]),
    ("tand(30+1i)", lambda: tand(30 + 1j), [0.5771158460013777 + 0.02326633209438335j]),
    ("atand(2+1i)", lambda: atand(2 + 1j), [67.5 + 9.928602006869321j]),
]


@pytest.mark.parametrize("expr, fn, expected", OCTAVE_COMPLEX_TRIG, ids=[c[0] for c in OCTAVE_COMPLEX_TRIG])
def test_complex_degree_trig_matches_octave(expr, fn, expected):
    got = np.atleast_1d(fn())
    expected = np.asarray(expected)
    assert got.shape == expected.shape
    np.testing.assert_allclose(got.real, expected.real, rtol=1e-13, atol=1e-15)
    np.testing.assert_allclose(got.imag, expected.imag, rtol=1e-13, atol=1e-15)
    # exact zeros of the real part where it vanishes mathematically
    if expr.startswith(("sind(1i)", "sind([180", "cosd(90", "cosd([-90")):
        assert got.real[0] == 0.0


def test_scalar_and_array_paths_agree():
    x = np.array([2.0, 0.5, -3.0, 1.0, -1.0])
    np.testing.assert_array_equal(asind(x), np.array([asind(v) for v in x], dtype=complex))
    z = np.array([1j, -30 + 1j, 180 + 1j, 45.0 + 0.25j])
    np.testing.assert_array_equal(sind(z), np.array([sind(v) for v in z]))
    np.testing.assert_array_equal(cosd(z), np.array([cosd(v) for v in z]))
    assert sind(1j) != 0  # exact-zero override only for real arguments
    assert np.imag(cosd(-30 + 1j)) > 0 and np.imag(cosd(30 + 1j)) < 0  # sign of the real part kept


def test_prandtl_meyer_of_subsonic_array_is_complex_like_octave():
    from MODELING.moc.get_prandtl_meyer_function import get_prandtl_meyer_function

    nu, mu = get_prandtl_meyer_function(np.array([0.5, 2.0]), 1.4)
    # Octave: nu = [0-23.59863932566508i, 26.37976081341647+7.8e-15i], mu = [90-75.45612929021688i, 30]
    np.testing.assert_allclose(nu, [-23.59863932566508j, 26.37976081341647], rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(mu, [89.99999999999997 - 75.45612929021688j, 30.0], rtol=1e-13, atol=1e-13)


def test_mmax_ignores_nan_and_orders_complex_by_modulus():
    assert mmax(np.nan, 1.0) == 1.0
    assert np.isnan(mmax(np.nan, np.nan))
    assert mmax(-4.0, 1 + 1j) == -4.0  # MATLAB max([-4, 1+1i]) = -4
    assert mmax(3.0, 1 + 5j) == 1 + 5j
    assert mmax(-3.0, 2.0) == 2.0
    assert mmax(np.complex128(-3 + 0j), 2.0) == 2.0  # zero imaginary part is narrowed


def test_relational_operators_on_complex():
    assert mlt(1 + 5j, 2)  # MATLAB: real parts
    with complex_ordering("octave"):
        assert not mlt(1 + 5j, 2)  # Octave: moduli
        assert not mle(-3 + 1j, 0)
    assert mle(-3 + 1j, 0)
    assert narrow(np.complex128(2 + 0j)) == 2.0 and isinstance(narrow(np.complex128(2 + 0j)), np.float64)


def test_mldivide2_matches_lapack_and_singular_gives_inf():
    rng = np.random.default_rng(0)
    for _ in range(50):
        A = rng.normal(size=(2, 2))
        b = rng.normal(size=2)
        x = mldivide2(A[0, 0], A[0, 1], A[1, 0], A[1, 1], b[0], b[1])
        np.testing.assert_allclose(x, np.linalg.solve(A, b), rtol=1e-12)
    x = mldivide2(1.0, 1.0, 1.0, 1.0, 1.0, 2.0)
    assert not np.all(np.isfinite(x))
    xc = mldivide2(1 + 1j, 2.0, 3.0, 4 - 1j, 1.0, 1j)
    A = np.array([[1 + 1j, 2], [3, 4 - 1j]])
    np.testing.assert_allclose(xc, np.linalg.solve(A, [1, 1j]), rtol=1e-12)


def test_fzero_bracketed():
    assert fzero(lambda x: x**2 - 2, (0, 2)) == pytest.approx(math.sqrt(2), abs=4e-16)
    assert fzero(lambda x: math.cos(x) - x, (0, 1)) == pytest.approx(0.7390851332151607, abs=1e-15)
    assert fzero(lambda x: x - 1, (1, 3)) == 1.0
    with pytest.raises(ValueError, match="differ in sign"):
        fzero(lambda x: x**2 + 1, (0, 1))
    r = fzero(lambda x: math.exp(x) - 5, (0, 30), tolx=1e-12)
    assert r == pytest.approx(math.log(5), abs=1e-11)


def test_grow_array_matlab_growth():
    A = GrowArray(np.array([[1.0]]))
    A[3, 2] = 5.0  # grows, zero padded
    np.testing.assert_array_equal(A.to_numpy(), [[1, 0], [0, 0], [0, 5]])
    with pytest.raises(IndexError):
        A[4, 1]
    A[1, 2] = 1 + 2j  # whole array becomes complex
    assert np.iscomplexobj(A.to_numpy())
    assert isinstance(A[3, 2], np.float64)  # indexing narrows
    v = GrowVector()
    v[3] = 7
    assert v.to_numpy().tolist() == [0, 0, 7]
    c = GrowColumns(3)
    c.set(2, [1, 2, 3])
    assert c.to_numpy().tolist() == [[0, 1], [0, 2], [0, 3]]


def test_grow_array_grows_each_axis_independently():
    """Only the axis whose capacity is exceeded grows (amortised doubling per axis)."""
    A = GrowArray(np.array([[1.0]]))
    for c in range(1, 8194):  # one row, 8 193 columns
        A[1, c] = float(c)
    assert A.shape == (1, 8193)
    assert A.capacity[0] == 16 and 8193 <= A.capacity[1] <= 2 * 8193
    assert A._buf.nbytes <= 16 * 16384 * 8  # 2 MiB (was 16384 x 16384 = 2 GiB)
    B = GrowArray()
    for r in range(1, 1001):  # one column, 1 000 rows
        B[r, 1] = float(r)
    assert B.shape == (1000, 1) and B.capacity[1] == 16 and B.capacity[0] <= 2000
    np.testing.assert_array_equal(B.to_numpy()[:, 0], np.arange(1, 1001))
    # number of reallocations stays logarithmic
    C = GrowArray()
    reallocs, last = 0, None
    for c in range(1, 5001):
        C[3, c] = 1.0
        if C._buf is not last:
            reallocs += 1
            last = C._buf
    assert reallocs <= 12
    assert C.to_numpy().shape == (3, 5000) and C.to_numpy()[:2].sum() == 0


def test_num2str():
    assert num2str(3.518326326701526) == "3.5183"
    assert num2str(12.3456789) == "12.3457"
    assert num2str(11.320612) == "11.3206"
    assert num2str(np.array([2, 3, 5])) == "2  3  5"
    assert num2str(7.0) == "7"
