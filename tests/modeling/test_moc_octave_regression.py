"""Node-by-node regression of the Python MoC port against the ORIGINAL MATLAB code run in GNU Octave.

Reference data: tests/modeling/data/*.mat, produced by
``tests/modeling/octave/generate_references.m`` (Octave 11.3.0) from the
unmodified sources in MODELING/MODELING_size_* (removed from the working tree
after the port; recover them from commit 2a30c95, see
``tests/matlab_sources.py``).  Tolerances are explicit below.

Complex values: in some axisymmetric cases the original free-pressure point
has no real solution in the far plume and MATLAB/Octave carry complex numbers.
MATLAB compares only real parts in ``<``/``>``; Octave compares moduli.  The
port follows MATLAB, so the full-array comparison runs under
``complex_ordering("octave")`` while the default (MATLAB) semantics are
compared on every column before the first complex one.
"""

from __future__ import annotations

import dataclasses
import subprocess
import warnings
from pathlib import Path

import numpy as np
import pytest
from _moc_ref import DATA, NOZZLE_CASES, inputs_from_reference, load, quiet, rel_err
from matlab_sources import OCTAVE, source_root, skip_reason

from MODELING.moc._errors import MocConvergenceError
from MODELING.moc._matlab import complex_ordering
from MODELING.moc.expansion import expansion
from MODELING.moc.get_Laval_theory import get_Laval_theory
from MODELING.moc.minimum_length_nozzle import minimum_length_nozzle
from MODELING.moc.MOC_2D_steady_irrotational_main import default_inputs, run
from MODELING.moc.oblique_shock import oblique_shock
from MODELING.moc.rde.MoC_Calculate_Domain_Size import MoC_Calculate_Domain_Size
from MODELING.moc.rde.MoC_Injection_Velocity import MoC_Injection_Velocity
from MODELING.moc.rde.MoC_Main_Model import default_inputs as rde_default_inputs
from MODELING.moc.rde.MoC_Resolve_Triple_Point import MoC_Resolve_Triple_Point

NODE_RTOL = 1e-10  # max|py - octave| / max|octave| for X, Y, U, V   (observed <= 4e-13)
MACH_ATOL = 1e-9  # Mach number                                      (observed <= 1.2e-12)
PRESSURE_RTOL = 1e-10  # p / P0                                      (observed <= 1e-13)
LAVAL_RTOL = 1e-9  # 1-D area-Mach inversion (fzero)                  (observed <= 6e-12)


def _solve(name, ordering="matlab"):
    geom, params = inputs_from_reference(name)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        with complex_ordering(ordering):
            sol, printed = quiet(run, geom, params, None)
    return sol, printed


def _first_complex_column(ref):
    imag = sum(np.abs(np.imag(ref[k])) for k in "XYUV")
    cols = np.flatnonzero(np.any(imag > 0, axis=0))
    return int(cols[0]) if cols.size else ref["X"].shape[1]


def test_default_inputs_match_matlab_script():
    """default_inputs() reproduces the input section of the team's MATLAB driver."""
    geom, params, plots = default_inputs()
    ref_geom, ref_params = inputs_from_reference("moc_team")
    for f in ("delta", "yt", "rhou", "rhod", "xe", "xplume", "ta", "te", "NI", "NIexpansion"):
        assert getattr(geom, f) == getattr(ref_geom, f), f
    np.testing.assert_array_equal(geom.circdownTheta, ref_geom.circdownTheta)
    assert dataclasses.astuple(params) == dataclasses.astuple(ref_params)
    assert (plots.patches, plots.patches_data, plots.patches_xlim) == (1, 0, 40.0)
    sol, _ = quiet(run, geom, params, None)
    r = load("moc_team", structs=True)
    np.testing.assert_allclose(sol.geom.circdownX, r["case_geom"].circdownX, rtol=1e-14, atol=1e-20)
    np.testing.assert_allclose(sol.geom.circdownY, r["case_geom"].circdownY, rtol=1e-15)


@pytest.mark.parametrize("name", NOZZLE_CASES)
def test_nozzle_nodes_match_octave(name):
    """Every node (x, y, u, v), Mach, p and every patch index equals the Octave run."""
    ref = load(name)
    sol, _ = _solve(name, ordering="octave")
    assert sol.indI == int(ref["indI"].item())
    np.testing.assert_array_equal(sol.LENG_INDI, ref["LENG_INDI"].ravel().astype(int))
    assert sol.ind_patch_tri == int(ref["ind_patch_tri"].item())
    assert sol.ind_patch_quad == int(ref["ind_patch_quad"].item())
    for py_name, ref_name in (
        ("patch_i_tri", "patch_i_tri"), ("patch_j_tri", "patch_j_tri"),
        ("patch_i_quad", "patch_i_quad"), ("patch_j_quad", "patch_j_quad"), ("LOC_PATCHES", "LOC_PATCHES"),
    ):
        np.testing.assert_array_equal(getattr(sol, py_name), ref[ref_name].astype(int), err_msg=py_name)
    for k in "XYUV":
        assert rel_err(getattr(sol, k), ref[k]) < NODE_RTOL, k
    Mach, pressure, *_ = sol.thermo()
    assert np.max(np.abs(Mach - ref["Mach"])) < MACH_ATOL
    assert rel_err(pressure, ref["pressure"]) < PRESSURE_RTOL
    np.testing.assert_allclose(sol.ythroat, ref["ythroat"].ravel(), rtol=1e-15, atol=1e-16)
    np.testing.assert_allclose(sol.uvnull, ref["uvnull"].ravel(), rtol=1e-15)


@pytest.mark.parametrize("name", NOZZLE_CASES)
def test_nozzle_matlab_semantics_match_octave_before_complex_region(name):
    """Default (MATLAB) semantics: identical to Octave up to the first complex-valued column."""
    ref = load(name)
    sol, _ = _solve(name, ordering="matlab")
    c0 = _first_complex_column(ref)
    np.testing.assert_array_equal(sol.LENG_INDI[:c0], ref["LENG_INDI"].ravel()[:c0].astype(int))
    for k in "XYUV":
        a = getattr(sol, k)
        b = ref[k]
        nr = min(a.shape[0], b.shape[0])
        assert rel_err(a[:nr, :c0], b[:nr, :c0]) < NODE_RTOL, k
    if c0 == ref["X"].shape[1]:  # no complex values at all: full equality incl. shape
        assert sol.X.shape == ref["X"].shape and not np.iscomplexobj(sol.X)


def test_console_messages_match_octave_order():
    """The disp() trace (shock/expansion messages, Mach numbers) is reproduced."""
    _, printed = _solve("moc_team")
    assert "    Mach number before the expansion: 2.5359" in printed
    assert "    Mach number after the expansion: 2.9057" in printed
    assert "    Deflection angle (degrees)       11.3206" in printed
    assert printed.count("Intersection of left-running characteristics inside the diverging region") == 12
    assert printed.count("Intersection of right-running characteristics inside the diverging region") == 5
    assert printed.count("Intersection of right-running characteristics inside the plume") == 3


def test_team_weak_shock_case_fails_like_original():
    """PRatio = 0.8 with the team geometry: the original stops with the predictor-corrector error."""
    ref = load("moc_team_weakshock")
    assert "maximum of iterations" in str(ref["error_message"])
    geom, params, _ = default_inputs()
    with pytest.raises(MocConvergenceError, match="maximum of iterations"):
        quiet(run, geom, dataclasses.replace(params, PRatio=0.8), None)


def test_strong_shock_stops_plume():
    """PRatio <= 0.5: plume stage returns its inputs unchanged (no plume columns)."""
    sol, printed = _solve("moc_team_strongshock")
    assert "Stopping the computation now..." in printed
    assert sol.LOC_PATCHES[3].tolist() == sol.LOC_PATCHES[2].tolist()


@pytest.mark.parametrize("name", ["moc_upstream", "moc_upstream_te0", "moc_upstream_NI21"])
def test_laval_theory_axisymmetric_matches_octave(name):
    ref = load(name)
    sol, _ = _solve(name)
    wall = ref["wallnodesJ"].ravel().astype(int) - 1
    ch = get_Laval_theory(sol.geom, sol.params, np.real(sol.X[0, wall]), np.real(sol.Y[0, wall]))
    # near A/A* = 1 the root is ill-conditioned (A/A*-1 ~ (M-1)^2): MATLAB-fzero port vs
    # Octave's fzero differ by ~1e-11 there
    np.testing.assert_allclose(ch.mach, ref["chocked_mach"].ravel(), rtol=LAVAL_RTOL)
    np.testing.assert_allclose(ch.pressure, ref["chocked_pressure"].ravel(), rtol=LAVAL_RTOL)


def test_laval_theory_planar_fix_vs_original():
    """Planar nozzle: the original used A ~ y^2 (bug); the port uses A ~ y (see PORTING_NOTES)."""
    ref = load("moc_team")
    sol, _ = _solve("moc_team")
    wall = ref["wallnodesJ"].ravel().astype(int) - 1
    y = np.real(sol.Y[0, wall])
    ch = get_Laval_theory(sol.geom, sol.params, np.real(sol.X[0, wall]), y)
    # Original (buggy) value at the last wall node is reproduced when y^2 is used
    ch_axi = get_Laval_theory(dataclasses.replace(sol.geom, delta=1), sol.params, np.real(sol.X[0, wall]), y)
    np.testing.assert_allclose(ch_axi.mach, ref["chocked_mach"].ravel(), rtol=LAVAL_RTOL)
    # last wall node (x = 9.39 yt): MoC wall Mach = 2.502
    assert ref["chocked_mach"].ravel()[-1] == pytest.approx(3.4262, abs=1e-4)  # before (bug)
    assert ch.mach[-1] == pytest.approx(2.5214, abs=1e-4)  # after (A/A* = y/y1)
    Mach, _, _, _, _ = sol.thermo()
    assert abs(ch.mach[-1] - Mach[0, wall[-1]]) < 0.03


def test_minimum_length_nozzle_original_mode_matches_octave():
    ref = load("mln", structs=True)
    res = minimum_length_nozzle(use_tangent_slopes=False)
    pairs = [("abs_", "abs"), ("ord", "ord"), ("theta", "theta"), ("nu", "nu"), ("mu", "mu"), ("Mach", "Mach"),
             ("Kminus", "Kminus"), ("Kplus", "Kplus"), ("Wall", "Wall"), ("absWall", "absWall"),
             ("ordWall", "ordWall"), ("absp", "absp"), ("ordp", "ordp"), ("Machp", "Machp")]
    for py_name, ref_name in pairs:
        np.testing.assert_allclose(getattr(res, py_name), ref[ref_name], rtol=1e-12, atol=1e-12, err_msg=py_name)
    assert res.area_ratio == pytest.approx(float(ref["area_ratio"]), rel=1e-13)
    assert res.theta_max == pytest.approx(float(ref["theta_max"]), rel=1e-14)


def test_expansion_matches_octave():
    ref = load("expansion", structs=True)
    res = expansion()
    assert res.delta == pytest.approx(float(ref["delta"]), rel=1e-11)  # fzero TolX = 1e-12
    for f in ("nu1", "mu1", "M2", "p1static", "p2static"):
        assert getattr(res, f) == pytest.approx(float(ref[f]), rel=1e-13), f


def test_oblique_shock_matches_octave():
    ref = load("oblique_shock")
    res, _ = quiet(oblique_shock)
    assert res.n_valid == [101, 101, 101]
    np.testing.assert_allclose(res.machmax[:, 0], ref["machmax"][:, 0], rtol=1e-14)
    # fminbnd (TolX = 1e-4 rad): beta at the maximum agrees to the optimiser tolerance,
    # theta_max (flat at the optimum) much better
    np.testing.assert_allclose(res.machmax[:, 1], ref["machmax"][:, 1], atol=1e-3)
    np.testing.assert_allclose(res.machmax[:, 2], ref["machmax"][:, 2], atol=1e-8)
    np.testing.assert_allclose(res.betadev, ref["betadev"], atol=1e-9)


def test_rde_injection_velocity_matches_octave():
    ref = load("rde", structs=True)
    chamber, inj, _ = rde_default_inputs()
    out, printed = quiet(MoC_Injection_Velocity, inj, chamber)
    np.testing.assert_allclose(out.V4, ref["inj_V4"], rtol=1e-14)
    np.testing.assert_allclose(out.P4, ref["inj_P4"], rtol=1e-14)
    np.testing.assert_allclose(out.T4, ref["inj_T4"], rtol=1e-14)
    np.testing.assert_array_equal(out.Vwave, ref["inj_Vwave"])  # [0, NaN]
    assert "1         121.752    1240476.626       860.328           0.201" in printed


def test_rde_triple_point_matches_octave():
    """The system has no root for the default inputs; both fsolve's stop at the same point."""
    ref = load("rde", structs=True)
    _, _, triple = rde_default_inputs()
    with pytest.warns(RuntimeWarning, match="did not converge"):
        beta, delta, P_match, M3, M1_p = MoC_Resolve_Triple_Point(triple)
    assert M1_p == pytest.approx(float(ref["M1_p"]), rel=1e-14)
    assert beta == pytest.approx(float(ref["beta"]), rel=1e-4)
    assert delta == pytest.approx(float(ref["delta"]), rel=1e-4)
    assert P_match == pytest.approx(float(ref["P_match"]), rel=1e-4)
    assert M3 == pytest.approx(float(ref["M3"]), rel=1e-4)


def test_rde_domain_size_matches_octave():
    ref = load("rde", structs=True)
    chamber, _, _ = rde_default_inputs()
    xMax, yMax = MoC_Calculate_Domain_Size(chamber)
    assert xMax == pytest.approx(float(ref["xMax"]), rel=1e-15)
    assert yMax == float(ref["yMax"])


TEAM_DIR = "MODELING_size_Nozzle_MOC_supersonic-master"


def _modeling_sources() -> Path:
    """``<root>/MODELING`` of the extracted tree; a root that already is that folder is accepted too."""
    root = source_root()
    return root if (root / TEAM_DIR).exists() and not (root / "MODELING" / TEAM_DIR).exists() else root / "MODELING"


ORIGINAL = _modeling_sources() / TEAM_DIR
LIVE_SKIP_REASON = skip_reason([ORIGINAL / "MOC_2D_steady_irrotational_main.m", ORIGINAL / "src2"])


@pytest.mark.skipif(LIVE_SKIP_REASON is not None, reason=LIVE_SKIP_REASON or "")
def test_live_octave_team_case(tmp_path):
    """Re-run the original team driver in Octave now and compare with Python."""
    script = (
        f"addpath('{DATA.parent / 'octave'}');"
        f"run_original('{ORIGINAL / 'MOC_2D_steady_irrotational_main.m'}', '{ORIGINAL / 'src2'}',"
        f" '{tmp_path / 'live.mat'}', {{}}, {{}}, {{'X','Y','U','V','LENG_INDI'}});"
    )
    subprocess.run([OCTAVE, "--no-gui", "--norc", "--eval", script], check=True, capture_output=True, timeout=300)
    import scipy.io as sio

    ref = sio.loadmat(tmp_path / "live.mat")
    geom, params, _ = default_inputs()
    sol, _ = quiet(run, geom, params, None)
    np.testing.assert_array_equal(sol.LENG_INDI, ref["LENG_INDI"].ravel().astype(int))
    for k in "XYUV":
        assert rel_err(getattr(sol, k), ref[k]) < NODE_RTOL, k
