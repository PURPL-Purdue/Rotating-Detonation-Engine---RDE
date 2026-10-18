"""Existing Python scripts (CDRConditions.py, UpstreamConditions.py, scfm.py) after clean-up.

``reference/*_stdout.txt`` is the stdout of the ORIGINAL scripts (run with the
Agg backend before the clean-up); the numeric goldens were read from the
original scripts' module namespaces.
"""

from pathlib import Path

import numpy as np
import pytest

from HADES.data_analysis import cdr_conditions as cdr
from HADES.data_analysis import scfm as scfm_mod
from HADES.data_analysis import upstream_conditions as up

REF = Path(__file__).parent / "reference"


# ------------------------------------------------------------------ CDR
@pytest.fixture(scope="module")
def cdr_results():
    return cdr.cdr_conditions(plot=False)


def test_cdr_printed_table_identical_to_original(cdr_results, capsys):
    cdr.print_results(cdr_results)
    assert capsys.readouterr().out == (REF / "CDRConditions_stdout.txt").read_text(encoding="utf-8")


def test_cdr_numbers_match_original(cdr_results):
    r = cdr_results
    assert r.P0_air == pytest.approx(19.648912924939477, rel=1e-10)
    np.testing.assert_allclose(r.P0_h2, [19.016002165963272, 21.39300243670868, 23.770002707454086,
                                         26.147002978199502, 28.5240032489449], rtol=1e-10)  # fmt: skip
    assert (r.P_in_air, r.P_ex_air) == pytest.approx((13.418349434528373, 9.806026709341902), rel=1e-10)
    assert (r.Ma_air, r.f_air) == pytest.approx((0.7582575849703135, 0.018330138117351683), rel=1e-10)
    np.testing.assert_allclose(r.P_in_h2, [14.085806555658882, 15.845932420760173, 17.606057497196563,
                                           19.36618199832953, 21.12630606692374], rtol=1e-10)  # fmt: skip
    np.testing.assert_allclose(r.Ma_h2, [0.6677696241917455, 0.6678135963583105, 0.6678488255735883,
                                         0.6678776836170667, 0.667901755515522], rtol=1e-10)  # fmt: skip
    assert (r.P_c_k_air, r.P_exit_k_air) == pytest.approx((10.37365463579381, 1.9261065055564282), rel=1e-10)
    np.testing.assert_allclose(r.P_i_arr, [1.885582329154349, 1.8811196806768364, 1.8767872088662747,
                                           1.8725838105305719, 1.8685083949079682], rtol=1e-10)  # fmt: skip
    np.testing.assert_allclose(r.avg_chamber_pressures[:, 0], [10.207033823206501, 10.182876615511686,
                               10.159424080119969, 10.13667024523217, 10.11460920633699], rtol=1e-10)  # fmt: skip
    # the cycle average does not depend on the number of waves
    np.testing.assert_allclose(r.avg_chamber_pressures, r.avg_chamber_pressures[:, :1].repeat(5, axis=1), rtol=1e-12)


def test_cdr_hand_checks():
    # Kessler: P_c = (2/(g+1))^(g/(g-1)) * P_p  (critical pressure ratio 0.5283 for g = 1.4)
    P_c, _ = cdr.kessler_P0(10.0, 0.1, 0.65, 1.4)
    assert P_c == pytest.approx(5.28282, rel=1e-5)
    # choked plenum pressure equals the generic choked-flow inversion
    from HADES.data_analysis.choked_flow import stagnation_pressure_for_mass_flow

    assert cdr.P0_plenum(0.9, 287.0, 293.15, 0.62, 3e-4, 1.4) == pytest.approx(
        stagnation_pressure_for_mass_flow(0.9, 0.62, 3e-4, 1.4, 287.0, 293.15) / 1e5, rel=1e-14
    )


def test_cdr_main_smoke(no_show, capsys):
    cdr.main()
    assert "RDE Average Chamber Pressures (bar):" in capsys.readouterr().out


# ------------------------------------------------------------- Upstream
@pytest.fixture(scope="module")
def grids():
    return up.build_grids(verbose=False)


def test_upstream_printed_ranges_identical_to_original(capsys):
    up.build_grids()
    assert capsys.readouterr().out == (REF / "UpstreamConditions_stdout.txt").read_text(encoding="utf-8")


def test_upstream_grids_match_original(grids):
    Z_air, Z_H2 = grids
    assert Z_air.shape == (200, 200)
    np.testing.assert_allclose([Z_air[0, 0], Z_air[-1, -1], Z_air[57, 123], Z_air.min(), Z_air.max()],
                               [25.076955995036805, 14.404449625785821, 16.9447318894078,
                                13.93164221946489, 25.92800932641448], rtol=1e-10)  # fmt: skip
    gold = {0.8: [26.616655201320675, 15.301785237837057, 17.990084207229717],
            1.0: [33.27081900165084, 19.12723154729632, 22.48760525903715],
            1.2: [39.92498280198101, 22.95267785675559, 26.98512631084458]}  # fmt: skip
    for phi, vals in gold.items():
        Z = Z_H2[phi]
        np.testing.assert_allclose([Z[0, 0], Z[-1, -1], Z[57, 123]], vals, rtol=1e-10)


def test_upstream_vectorised_grid_equals_cellwise_original(grids):
    Z_air, Z_H2 = grids
    for i, j in [(0, 0), (13, 177), (199, 42)]:
        assert Z_air[i, j] == up.calc_Pp_air(float(up.T0_GRID[i]), float(up.CD_GRID[j]))
        assert Z_H2[1.1][i, j] == up.calc_Pp_H2(float(up.T0_GRID[i]), float(up.CD_GRID[j]), 1.1)


def test_upstream_heatmaps_smoke(grids):
    fig1, fig2 = up.plot_heatmaps(*grids)
    assert len(fig2.axes) >= 5


# ------------------------------------------------------------------ SCFM
def test_scfm_identical_to_original(capsys):
    q_air, q_h2 = scfm_mod.scfm()
    assert q_air == 1568.4539638230854
    assert q_h2 == 788.871432196022
    scfm_mod.main()
    assert capsys.readouterr().out == (REF / "scfm_stdout.txt").read_text()


def test_scfm_hand_calculation():
    # ideal-gas density at 750 psia, 293 K; 1 m^3/s = 2118.88 cfm; SCFM = ACFM * P/14.7 * 519/T_R
    rho = 5.171e6 / (287 * 293)
    acfm = 0.907185 / rho * 2118.8799727597
    assert scfm_mod.scfm()[0] == pytest.approx(acfm * (5.171e6 * 0.000145038 / 14.7) * 519 / 527.73, rel=1e-14)
    assert scfm_mod.scfm()[0] == pytest.approx(1568.45, rel=1e-5)
