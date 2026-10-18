"""Post-processing figures and script entry points (non-interactive Agg backend)."""

from __future__ import annotations

import dataclasses
import warnings

import matplotlib

matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402
from _moc_ref import load, quiet  # noqa: E402

from MODELING.moc.MOC_2D_steady_irrotational_main import default_inputs, run  # noqa: E402
from MODELING.moc.MOC_2D_steady_irrotational_postprocess import axis_nodes, build_patches  # noqa: E402


@pytest.fixture(autouse=True)
def _no_show(monkeypatch):
    monkeypatch.setattr(plt, "show", lambda *a, **k: None)
    yield
    plt.close("all")


UPSTREAM = dict(delta=1, rhou=2.0, rhod=0.5, xplume=20.0, NIexpansion=10)  # unchanged upstream inputs


@pytest.mark.parametrize("inputs", ["team", "upstream"])
@pytest.mark.parametrize("patches_data", range(7))
def test_postprocess_all_patch_modes(inputs, patches_data):
    """Every plots.patches_data option renders for the team inputs and for the upstream
    inputs, whose retained plume column is complex-valued (the original completes in Octave)."""
    geom, params, plots = default_inputs()
    if inputs == "upstream":
        geom = dataclasses.replace(geom, **UPSTREAM)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        sol, printed = quiet(run, geom, params, dataclasses.replace(plots, patches_data=patches_data))
    assert np.iscomplexobj(sol.X) == (inputs == "upstream")
    assert "Postprocessing the results..." in printed
    ax1 = sol.fig1.axes[0]
    x, y = ax1.lines[0].get_data()  # wall line X(1,1:indI)/yt
    np.testing.assert_allclose(x, np.real(sol.X[0, : sol.indI]) / geom.yt)
    if patches_data == 0:
        assert len(ax1.collections) == 2  # tri + quad outlines
    else:
        assert len(sol.fig1.axes) == 2  # colour bar
        assert np.all(np.isfinite(ax1.collections[0].get_array()))
    assert [t.get_text() for t in sol.fig2.axes[0].get_legend().get_texts()] == ["Wall", "Axis", "1D"]
    assert sol.fig2.axes[0].get_yscale() == "log"
    assert len(sol.fig2.axes[0].lines) == 3 and len(sol.fig3.axes[0].lines) == 3


def test_flow_angle_patch_data_on_complex_nodes():
    """patches_data = 6 is atand(V./U) (complex-capable); patches show its real part."""
    from MODELING.moc._matlab import atand

    geom, params, plots = default_inputs()
    geom = dataclasses.replace(geom, **UPSTREAM)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        sol, _ = quiet(run, geom, params, None)
    with np.errstate(divide="ignore", invalid="ignore"):
        angle = atand(sol.V / sol.U)
    assert np.iscomplexobj(angle) and np.any(np.imag(angle) != 0)
    pd = build_patches(sol.X, sol.Y, angle, sol.ind_patch_tri, sol.ind_patch_quad, sol.patch_i_tri,
                       sol.patch_j_tri, sol.patch_i_quad, sol.patch_j_quad, sol.LOC_PATCHES, 40.0)
    used = pd.faces_quad.ravel()
    assert np.all(np.isfinite(pd.nodesC_quad[used]))
    # wall nodes inside the nozzle: flow angle equals the wall angle (te = 15 deg at the lip)
    lip = int(np.flatnonzero(np.real(sol.X[0]) == geom.xe)[0])
    assert np.real(angle[0, lip]) == pytest.approx(geom.te, abs=1e-6)


def test_postprocess_unknown_patch_mode():
    geom, params, plots = default_inputs()
    with pytest.raises(ValueError, match="Unknown type"):
        quiet(run, geom, params, dataclasses.replace(plots, patches_data=7))


def test_axis_nodes_and_wall_curves_match_matlab_selection():
    """Axis nodes = last nonzero U of each column (MATLAB find(U(:,I),1,'last')), column order."""
    ref = load("moc_team")
    geom, params, plots = default_inputs()
    sol, _ = quiet(run, geom, params, plots)
    cond = axis_nodes(sol.U)
    assert cond.sum() == sol.X.shape[1]
    ax_x = sol.fig2.axes[0].lines[1].get_xdata()
    expected = [sol.X[np.flatnonzero(sol.U[:, I])[-1], I] for I in range(sol.X.shape[1])]
    np.testing.assert_allclose(ax_x, np.asarray(expected) / geom.yt)
    wall = ref["wallnodesJ"].ravel().astype(int) - 1
    np.testing.assert_allclose(sol.fig3.axes[0].lines[0].get_xdata(), sol.X[0, wall] / geom.yt)


def test_build_patches_regions_and_xlim():
    geom, params, plots = default_inputs()
    sol, _ = quiet(run, geom, params, None)
    Mach, *_ = sol.thermo()
    pd = build_patches(sol.X, sol.Y, Mach, sol.ind_patch_tri, sol.ind_patch_quad, sol.patch_i_tri,
                       sol.patch_j_tri, sol.patch_i_quad, sol.patch_j_quad, sol.LOC_PATCHES, 5.0)
    assert pd.nodes_tri.shape == (3 * sol.ind_patch_tri, 2)
    assert np.all(pd.nodes_tri[pd.faces_tri.ravel(), 0] < 5.0)
    # region colours: black (initial line) for the first LOC_PATCHES[0,0] triangles
    np.testing.assert_array_equal(pd.nodesRGB_tri[: 3 * sol.LOC_PATCHES[0, 0]], 0.0)
    np.testing.assert_array_equal(pd.nodesRGB_tri[3 * sol.LOC_PATCHES[0, 0]], [1, 0, 0])


def test_script_mains_run():
    from MODELING.moc import MOC_2D_steady_irrotational_main, expansion, minimum_length_nozzle, oblique_shock
    from MODELING.moc.rde import MoC_Main_Model

    sol, printed = quiet(MOC_2D_steady_irrotational_main.main)
    assert sol.fig1 is not None
    res, printed = quiet(expansion.main)
    assert printed.strip() == "delta = 2.1145"
    res, printed = quiet(minimum_length_nozzle.main)
    assert printed.startswith("area_ratio = 2.40448")
    res, _ = quiet(oblique_shock.main)
    assert res.betadev.shape == (3, 101, 3)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        (out, tp, post), printed = quiet(MoC_Main_Model.main)
    assert "Beta Shock: 70.34° | Slip Angle: 42.15°" in printed
    assert printed.rstrip().endswith("Graphing Results...")
    ax = plt.figure("RDE Flow Field").axes[0]
    assert ax.get_xlim() == (0.0, post.xMax) and ax.get_ylim() == (0.0, post.yMax)
    assert ax.get_box_aspect() == pytest.approx(post.yMax / post.xMax)
