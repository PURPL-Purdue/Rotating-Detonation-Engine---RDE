"""rdelib.cea (NASA CEA 3 via the `cea` package) vs the committed CEA2 outputs.

Every value printed in the committed CEA2 outputs {cea_det,cea_run,h2ophi1}.out
(copied from HADES/SIZING MODEL/CEA, now removed, into reference/)
is compared at the precision CEA2 printed it (half a unit in the last digit),
except where noted.
"""

from __future__ import annotations

import numpy as np
import pytest

from rdelib import cea as rcea
from rdelib.units import CEA2_PSIA_TO_BAR


def _close(x, printed, decimals):
    """|x - printed| <= half a unit of the last printed decimal (+ float slack)."""
    return abs(x - printed) <= 0.5 * 10.0 ** (-decimals) + 1e-12


@pytest.fixture(scope="module")
def det():
    # cea_det.inp: det, p,psia=26, t,k=283.15, phi=0.9, Air/H2, output transport
    return rcea.cj_detonation("H2", "Air", 26.0 * CEA2_PSIA_TO_BAR, 283.15, phi=0.9, transport=True)


@pytest.fixture(scope="module")
def rocket():
    # cea_run.inp: rocket equilibrium, p,psia=142.745149, supar=10, phi=0.9, no reactant T (h=0)
    return rcea.rocket("H2", "Air", 142.745149 * CEA2_PSIA_TO_BAR, supar=10.0, phi=0.9)


# value printed in cea_det.out, decimals printed
DET_EXPECTED = {
    "P1_bar": (1.7926, 4), "T1_K": (283.15, 2), "H1_kJkg": (-24.38, 2), "M1": (21.586, 3),
    "gamma1": (1.4023, 4), "son_vel1": (391.1, 1),
    "P_bar": (28.927, 3), "T_K": (2904.62, 2), "rho": (2.9526, 4), "H_kJkg": (1260.53, 2),
    "U_kJkg": (280.83, 2), "G_kJkg": (-28171.0, 1), "S_kJkgK": (10.1326, 4), "M": (24.651, 3),
    "cp_eq": (2.7547, 4), "gamma_s": (1.1779, 4), "son_vel": (1074.2, 1),
    "visc_mP": (0.91685, 5), "Pr_eq": (0.4847, 4),
    "cp_frozen": (1.7219, 4), "k_frozen": (2.2983, 4), "Pr_frozen": (0.6869, 4),
    "P_P1": (16.136, 3), "T_T1": (10.258, 3), "M_M1": (1.1420, 4), "rho_rho1": (1.7964, 4),
    "det_mach": (4.9344, 4), "det_vel": (1929.7, 1),
}


@pytest.mark.parametrize("field,expected", list(DET_EXPECTED.items()))
def test_detonation_matches_cea2_printed_value(det, field, expected):
    value, dec = expected
    assert _close(getattr(det, field), value, dec), (field, getattr(det, field), value)


def test_detonation_conductivity_eq_within_2e5(det):
    # CEA2 printed 5.2103; CEA3 gives 5.21016 (2e-5 relative). Documented gap.
    assert det.k_eq == pytest.approx(5.2103, rel=3e-5)


def test_detonation_mole_fractions(det):
    printed = {"Ar": 0.00777, "CO": 0.00007, "CO2": 0.00020, "H": 0.00259, "HO2": 0.00002,
               "H2": 0.01418, "H2O": 0.28840, "NO": 0.01013, "NO2": 0.00001, "N2": 0.64247,
               "O": 0.00182, "OH": 0.01759, "O2": 0.01474}
    for sp, x in printed.items():
        assert _close(det.mole_fractions[sp], x, 5), sp
    assert det.of_ratio == pytest.approx(38.10692, abs=5e-6)


ROCKET_EXPECTED = {  # cea_run.out columns: chamber, throat, exit
    "P_bar": ([9.8419, 5.4862, 0.09844], [4, 4, 5]),
    "T_K": ([2330.13, 2094.92, 875.51], [2, 2, 2]),
    "H_kJkg": ([0.0, -430.42, -2330.67], [5, 2, 2]),
    "U_kJkg": ([-776.51, -1127.50, -2621.74], [2, 2, 2]),
    "G_kJkg": ([-23330.6, -21405.9, -11096.8], [1, 1, 1]),
    "S_kJkgK": ([10.0126] * 3, [4] * 3),
    "M": ([24.950, 24.987, 25.010], [3] * 3),
    "cp_eq": ([1.9528, 1.7828, 1.3736], [4] * 3),
    "gamma_s": ([1.2177, 1.2349, 1.3193], [4] * 3),
    "son_vel": ([972.4, 927.8, 619.7], [1] * 3),
    "mach": ([0.0, 1.0, 3.484], [3] * 3),
}


# Documented sub-print-precision differences (CEA3 vs CEA2 convergence):
#  * throat T: CEA3 2094.91497 K vs printed 2094.92 (value sits on the rounding
#    boundary; difference >= 3e-5 K, i.e. 1.5e-8 relative)
#  * chamber H: CEA3 1.4e-4 kJ/kg vs exactly 0 (hc = 0 assigned; solver residual)
ROCKET_SLACK = {("T_K", 1): 0.01, ("H_kJkg", 0): 1e-3}


@pytest.mark.parametrize("field,expected", list(ROCKET_EXPECTED.items()))
def test_rocket_matches_cea2_printed_values(rocket, field, expected):
    values, decs = expected
    arr = getattr(rocket, field)
    for i, (v, p, d) in enumerate(zip(arr, values, decs)):
        slack = ROCKET_SLACK.get((field, i))
        if slack is not None:
            assert abs(v - p) <= slack, (field, i, v, p)
        else:
            assert _close(v, p, d), (field, i, v, p)


def test_rocket_density_and_pressure_ratio(rocket):
    # RHO printed in E-format: 1.2675 0, 7.8702-1, 3.3821-2
    np.testing.assert_allclose(rocket.rho, [1.2675, 0.78702, 0.033821], rtol=5e-5)
    assert np.all(np.abs(rocket.pinf_p - [1.0, 1.7939, 99.978]) <= [5e-5, 5e-5, 5e-4])


def test_rocket_performance_parameters(rocket):
    # PERFORMANCE PARAMETERS: throat, exit
    t, e = 1, 2
    assert _close(rocket.ae_at[t], 1.0, 4) and _close(rocket.ae_at[e], 10.0, 3)
    assert _close(rocket.cstar[t], 1347.8, 1) and _close(rocket.cstar[e], 1347.8, 1)
    assert _close(rocket.cf[t], 0.6884, 4) and _close(rocket.cf[e], 1.6019, 4)
    assert _close(rocket.ivac[t], 1679.1, 1) and _close(rocket.ivac[e], 2293.8, 1)
    assert _close(rocket.isp[t], 927.8, 1) and _close(rocket.isp[e], 2159.0, 1)


def test_rocket_mole_fractions(rocket):
    printed = {"H2O": [0.31207, 0.31524, 0.31725], "N2": [0.65335, 0.65514, 0.65698],
               "O2": [0.01526, 0.01609, 0.01762], "OH": [0.00471, 0.00219, 0.0],
               "NO": [0.00410, 0.00248, 0.0], "Ar": [0.00786, 0.00787, 0.00788],
               "H2": [0.00198, 0.00061, 0.0], "CO2": [0.00026, 0.00027, 0.00027]}
    for sp, xs in printed.items():
        for v, x in zip(rocket.mole_fractions[sp], xs):
            assert _close(v, x, 5), sp


def test_h2ophi1_hp_only_h2o_public_api():
    """h2ophi1.out: H2/O2 (298.15 K), phi = 1, hp at 1 atm, products 'only H2O'.

    CEA3 cannot attach a reactant list longer than the product list (CEA2
    reduced the components instead); rdelib.cea.equilibrium validates the
    elemental feasibility in Python and solves the equivalent problem
    (same elemental composition and reactant enthalpy).
    """
    r = rcea.equilibrium("H2", "O2", 1.01325, 298.15, problem="hp", phi=1.0, only=["H2O"])
    assert r.converged
    assert r.of_ratio == pytest.approx(7.936683, abs=5e-7)
    assert _close(r.T_K, 4930.56, 2)
    assert r.rho == pytest.approx(4.4527e-2, rel=2e-5)
    assert _close(r.U_kJkg, -2275.58, 2)
    assert _close(r.G_kJkg, -86535.6, 1)
    assert _close(r.S_kJkgK, 17.5509, 4)
    assert _close(r.M, 18.015, 3)
    assert _close(r.cp_eq, 3.3823, 4)
    assert _close(r.gamma_s, 1.1580, 4)
    assert _close(r.son_vel, 1623.3, 1)
    assert r.mole_fractions["H2O"] == pytest.approx(1.0)


_GUARD_SCRIPT = r"""
import json
from rdelib import cea as rc
cases = {
    "nonstoich_only_H2O": lambda: rc.equilibrium("H2", "O2", 1.01325, 298.15, phi=0.9, only=["H2O"]),
    "foreign_element": lambda: rc.equilibrium("H2", "O2", 1.01325, 298.15, phi=1.0, only=["H2O", "CO2"]),
    "missing_element": lambda: rc.equilibrium("H2", "O2", 1.01325, 298.15, phi=1.0, only=["H2", "H"]),
    "air_missing_N_Ar_C": lambda: rc.equilibrium("H2", "Air", 1.0, 298.15, phi=1.0, only=["H2O", "O2", "H2"]),
    "det_fewer_products": lambda: rc.cj_detonation("H2", "O2", 1.0, 298.15, phi=1.0, only=["H2O"]),
    "rocket_fewer_products": lambda: rc.rocket("H2", "O2", 10.0, supar=5.0, phi=1.0, only=["H2O"]),
    "rocket_foreign_element": lambda: rc.rocket("H2", "O2", 10.0, supar=5.0, phi=1.0, only=["H2O", "OH", "CO"]),
    "bad_species": lambda: rc.equilibrium("H2", "O2", 1.0, 298.15, phi=1.0, only=["H2Oo"]),
}
out = {}
for k, f in cases.items():
    try:
        f()
        out[k] = "NO ERROR"
    except Exception as e:
        out[k] = type(e).__name__
print("RESULT" + json.dumps(out))
"""


def test_restricted_products_raise_python_errors_not_native_crashes():
    """Inputs that crash/abort the native CEA solver must raise Python exceptions.

    Run in a subprocess so that a regression shows up as a failure instead of
    killing the test session.
    """
    import json
    import subprocess
    import sys

    p = subprocess.run([sys.executable, "-c", _GUARD_SCRIPT], capture_output=True, text=True, timeout=120)
    assert p.returncode == 0, p.stderr[-2000:]
    line = [ln for ln in p.stdout.splitlines() if ln.startswith("RESULT")][-1]
    res = json.loads(line[len("RESULT"):])
    assert res == {
        "nonstoich_only_H2O": "ValueError",
        "foreign_element": "ValueError",
        "missing_element": "ValueError",
        "air_missing_N_Ar_C": "ValueError",
        "det_fewer_products": "NotImplementedError",
        "rocket_fewer_products": "NotImplementedError",
        "rocket_foreign_element": "ValueError",
        "bad_species": "RuntimeError",
    }


def test_restricted_products_valid_cases():
    # two-product set with reactants attached (supported natively) gives the same answer
    r2 = rcea.equilibrium("H2", "O2", 1.01325, 298.15, problem="hp", phi=1.0, only=["H2O", "H2"])
    assert r2.T_K == pytest.approx(4930.564, abs=1e-2)
    # restricted CJ detonation with an element-complete product list
    full = rcea.cj_detonation("H2", "Air", 1.0, 298.15, phi=1.0, transport=False)
    small = rcea.cj_detonation("H2", "Air", 1.0, 298.15, phi=1.0, transport=False,
                               only=["H2O", "N2", "O2", "H2", "Ar", "CO2", "OH", "H", "O", "NO"])
    assert small.det_vel == pytest.approx(full.det_vel, rel=2e-3)


@pytest.mark.parametrize("name,expected", [
    ("H2O", {"H": 2, "O": 1}), ("H2O(L)", {"H": 2, "O": 1}), ("C(gr)", {"C": 1}),
    ("CH3C(CH3)2CH3", {"C": 5, "H": 12}), ("C2H2,acetylene", {"C": 2, "H": 2}),
    ("(HCOOH)2", {"C": 2, "H": 4, "O": 4}), ("HO(CO)2OH", {"C": 2, "H": 2, "O": 4}), ("Ar", {"Ar": 1}),
])
def test_species_elements(name, expected):
    assert rcea.species_elements(name) == expected


@pytest.mark.parametrize("name", ["Air", "RP-1", "H2O)", "(H2O", "Xx2"])
def test_species_elements_rejects_non_formulas(name):
    with pytest.raises(ValueError):
        rcea.species_elements(name)


def test_equilibrium_hp_conserves_enthalpy_and_is_plausible():
    r = rcea.equilibrium("H2", "O2", 1.01325, 298.15, problem="hp", phi=1.0)
    assert r.H_kJkg == pytest.approx(0.0, abs=1e-3)  # reactants at 298.15 K have h ~ 0
    assert 3050.0 < r.T_K < 3110.0  # stoichiometric H2/O2 adiabatic flame temperature, 1 atm


def test_stoichiometric_h2_air_cj_literature():
    # Textbook CJ values for stoichiometric H2-air at 1 atm, 298 K: D_CJ ~ 1968 m/s,
    # P_CJ/P1 ~ 15.6, T_CJ ~ 2950 K.
    d = rcea.cj_detonation("H2", "Air", 1.01325, 298.15, phi=1.0)
    assert d.det_vel == pytest.approx(1968.0, abs=10.0)
    assert d.P_P1 == pytest.approx(15.6, abs=0.3)
    assert d.T_K == pytest.approx(2950.0, abs=30.0)


def test_mixture_specifications_are_consistent():
    base = rcea.cj_detonation("H2", "Air", 1.8, 290.0, phi=0.9, transport=False)
    via_of = rcea.cj_detonation("H2", "Air", 1.8, 290.0, of=base.of_ratio, transport=False)
    pct = 100.0 / (1.0 + base.of_ratio)
    via_pct = rcea.cj_detonation("H2", "Air", 1.8, 290.0, pct_fuel=pct, transport=False)
    for other in (via_of, via_pct):
        assert other.det_vel == pytest.approx(base.det_vel, rel=1e-10)
    with pytest.raises(ValueError):
        rcea.cj_detonation("H2", "Air", 1.8, 290.0, phi=0.9, of=30.0)


def test_cea2_round():
    assert rcea.cea2_round(1929.7153, 1) == 1929.7
    assert rcea.cea2_round(28.926726, "var") == 28.927
    assert rcea.cea2_round(0.9168523, "var") == 0.91685
    assert rcea.cea2_round(-28170.969, "var") == -28171.0
    assert rcea.cea2_round(1260.5299, "var") == 1260.53
    assert rcea.cea2_round(5.2101620, "var") == 5.2102
