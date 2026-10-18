"""Physics / theory checks of the MoC port (independent of the Octave reference)."""

from __future__ import annotations

import dataclasses
import math
import warnings

import numpy as np
import pytest
from _moc_ref import quiet

from MODELING.moc.expansion import expansion
from MODELING.moc.get_Laval_theory import get_Laval_theory
from MODELING.moc.get_Mach_from_nu import get_Mach_from_nu
from MODELING.moc.get_prandtl_meyer_function import get_prandtl_meyer_function
from MODELING.moc.minimum_length_nozzle import minimum_length_nozzle
from MODELING.moc.MOC_2D_steady_irrotational_get_thermo2 import MOC_2D_steady_irrotational_get_thermo2
from MODELING.moc.MOC_2D_steady_irrotational_IVLINE import MOC_2D_steady_irrotational_IVLINE
from MODELING.moc.MOC_2D_steady_irrotational_main import default_inputs, run
from MODELING.moc.oblique_shock import equation_Theta_Beta_Mach, oblique_shock
from MODELING.moc.rde.MoC_Calculate_Domain_Size import MoC_Calculate_Domain_Size
from MODELING.moc.rde.MoC_Injection_Velocity import MoC_Injection_Velocity
from MODELING.moc.rde.MoC_Main_Model import default_inputs as rde_default_inputs
from MODELING.moc.rde.MoC_Resolve_Triple_Point import triple_point_residuals


def area_ratio(M, g):
    """Isentropic A/A*."""
    return ((2 / (g + 1)) * (1 + 0.5 * (g - 1) * M * M)) ** ((g + 1) / (2 * (g - 1))) / M


# --------------------------------------------------------------------------- Prandtl-Meyer
@pytest.mark.parametrize(
    "M, nu_deg, mu_deg",
    [(1.0, 0.0, 90.0), (2.0, 26.3798, 30.0), (3.0, 49.7573, 19.4712)],  # NACA 1135 tables
)
def test_prandtl_meyer_tabulated(M, nu_deg, mu_deg):
    nu, mu = get_prandtl_meyer_function(M, 1.4)
    assert nu == pytest.approx(nu_deg, abs=1e-4)
    assert mu == pytest.approx(mu_deg, abs=1e-4)
    s = math.sqrt(2.4 / 0.4)
    r = math.sqrt(M * M - 1)
    assert nu == pytest.approx(math.degrees(s * math.atan(r / s) - math.atan(r)), abs=1e-12)


@pytest.mark.parametrize("gamma", [1.2, 1.4, 1.67])
def test_mach_from_nu_inverts_prandtl_meyer(gamma):
    M = np.array([1.01, 1.5, 2.4, 3.7, 6.0])
    nu, _ = get_prandtl_meyer_function(M, gamma)
    np.testing.assert_allclose(get_Mach_from_nu(nu, gamma), M, rtol=1e-12)
    assert isinstance(get_Mach_from_nu(float(nu[2]), gamma), float)


def test_mach_from_nu_asserts():
    with pytest.raises(AssertionError):
        get_Mach_from_nu(0.0, 1.4)
    with pytest.raises(AssertionError):
        get_Mach_from_nu(126.0, 1.4)


# --------------------------------------------------------------------------- 1-D relations
@pytest.mark.parametrize("delta", [0, 1])
def test_laval_theory_area_mach(delta):
    g = 1.4
    geom, params, _ = default_inputs()
    geom = dataclasses.replace(geom, delta=delta)
    params = dataclasses.replace(params, gamma=g)
    M = np.array([1.0, 1.5, 2.0, 3.0])
    ar = area_ratio(M, g)
    Y = ar ** (1.0 / (1 + delta))  # y/y* for planar / axisymmetric
    ch = get_Laval_theory(geom, params, np.zeros_like(Y), Y)
    np.testing.assert_allclose(ch.mach, M, rtol=1e-10)
    np.testing.assert_allclose(ch.pressure, (1 + 0.2 * M**2) ** -3.5, rtol=1e-10)


def test_ivline_sonic_on_axis():
    """Initial-value line: x_vnull = -eps at y = 0 and u = a* there; the sonic line lies upstream."""
    geom, params, _ = default_inputs()
    y = np.linspace(0, geom.yt, 11)
    xs, xv, u = MOC_2D_steady_irrotational_IVLINE(geom, params, y)
    astar = math.sqrt(2 * params.gamma * params.R * params.T / (params.gamma + 1))
    assert u[0] == pytest.approx(astar, rel=1e-15)
    assert np.all(u >= astar)
    assert np.all(xs[1:] < xv[1:])


# --------------------------------------------------------------------------- full nozzle
def _solve(**geom_changes):
    geom, params, _ = default_inputs()
    geom = dataclasses.replace(geom, **geom_changes)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        sol, _ = quiet(run, geom, params, None)
    return sol


UPSTREAM = dict(delta=1, rhou=2.0, rhod=0.5, xplume=20.0, NIexpansion=10)  # unchanged upstream inputs


def _mass_flow(sol, I):
    """Mass flow through the C- characteristic stored in column I (0-based)."""
    n = int(sol.LENG_INDI[I])
    x, y, u, v = (np.real(a[:n, I]) for a in (sol.X, sol.Y, sol.U, sol.V))
    _, _, rho, _, _ = (np.real(q[:n, I]) for q in sol.thermo())
    w = 2 * np.pi * y if sol.geom.delta == 1 else np.ones_like(y)
    f = rho * w
    fu = 0.5 * (f[1:] * u[1:] + f[:-1] * u[:-1])
    fv = 0.5 * (f[1:] * v[1:] + f[:-1] * v[:-1])
    return -np.sum(fu * np.diff(y) - fv * np.diff(x))


def _mdot_1d(sol):
    g, p = sol.geom, sol.params
    gm = p.gamma
    Ts = p.T * 2 / (gm + 1)
    ps = p.P * (2 / (gm + 1)) ** (gm / (gm - 1))
    A = g.yt if g.delta == 0 else math.pi * g.yt**2
    return ps / (p.R * Ts) * math.sqrt(gm * p.R * Ts) * A


def _exit_column(sol):
    return int(np.flatnonzero(np.real(sol.X[0]) == sol.geom.xe)[0])


def test_mass_conservation_axisymmetric():
    """Upstream axisymmetric case: mass flow along successive characteristics within 1 % of 1-D choked flow."""
    sol = _solve(**UPSTREAM)
    md = _mdot_1d(sol)
    ex = _exit_column(sol)
    for I in (sol.geom.NI - 1, sol.geom.NI + 5, sol.geom.NI + 20, ex, ex + 10):
        assert _mass_flow(sol, I) / md == pytest.approx(1.0, abs=0.01)


def test_mass_conservation_team_planar():
    """Team planar case (sharp corner, many weak-shock patch deletions): within 7 % of 1-D choked flow."""
    sol = _solve()
    md = _mdot_1d(sol)
    ratios = [_mass_flow(sol, I) / md for I in (sol.geom.NI - 1, sol.geom.NI + 20, _exit_column(sol))]
    assert ratios[0] == pytest.approx(0.977, abs=0.005)  # transonic initial line (2-D throat Cd)
    for r in ratios:
        assert r == pytest.approx(1.0, abs=0.07)


@pytest.mark.parametrize("changes", [{}, UPSTREAM])
def test_exit_lip_mach_vs_isentropic_area_ratio(changes):
    """Mach at the exit lip within 2 % of the 1-D value for the exit area ratio."""
    sol = _solve(**changes)
    ex = _exit_column(sol)
    ar = (np.real(sol.Y[0, ex]) / sol.geom.yt) ** (1 + sol.geom.delta)
    ch = get_Laval_theory(sol.geom, sol.params, [0.0, 0.0], [1.0, ar ** (1 / (1 + sol.geom.delta))])
    Mach, *_ = sol.thermo()
    assert np.real(Mach[0, ex]) == pytest.approx(ch.mach[1], rel=0.02)


@pytest.mark.parametrize("changes", [{}, UPSTREAM])
def test_plume_boundary_is_at_ambient_pressure_and_pm_fan(changes):
    """Jet boundary at p_amb = p_lip / PRatio; lip fan obeys nu(M2) - nu(M1) = deflection."""
    sol = _solve(**changes)
    ex = _exit_column(sol)
    Mach, p, *_ = sol.thermo()
    p_lip = np.real(p[0, ex])
    p_amb = p_lip / sol.params.PRatio
    jet = [I for I in range(ex + 1, sol.indI) if np.real(sol.X[0, I]) > sol.geom.xe and np.imag(sol.X[0, I]) == 0]
    assert len(jet) > 5
    np.testing.assert_allclose(np.real(p[0, jet]), p_amb, rtol=1e-12)
    # Prandtl-Meyer: Mach on the jet boundary is the fully expanded Mach number
    M1 = np.real(Mach[0, ex])
    M2 = np.real(Mach[0, jet[0]])
    nu1, _ = get_prandtl_meyer_function(M1, sol.params.gamma)
    nu2, _ = get_prandtl_meyer_function(M2, sol.params.gamma)
    p2, _, _, _ = MOC_2D_steady_irrotational_get_thermo2(M2, sol.params)
    assert p2 == pytest.approx(p_amb, rel=1e-12)
    assert nu2 - nu1 > 0
    # flow angle on the first jet point ~ lip angle + deflection
    theta_jet = math.degrees(math.atan2(np.real(sol.V[0, jet[0]]), np.real(sol.U[0, jet[0]])))
    assert theta_jet == pytest.approx(sol.geom.te + (nu2 - nu1), abs=1.0)


def test_characteristic_net_geometry():
    """Row 1 is the wall (y = contour), last row of each column is the axis (y = 1e-6)."""
    sol = _solve()
    for I in range(sol.indI):
        n = int(sol.LENG_INDI[I])
        if n > 1:
            assert np.real(sol.Y[n - 1, I]) == pytest.approx(1e-6, abs=1e-9) or np.real(sol.Y[n - 1, I]) > 0
    ex = _exit_column(sol)
    # straight conical wall (ta = te): y = yt + tan(15 deg) * x on the nozzle wall
    xs = np.real(sol.X[0, sol.geom.NI + 15 : ex + 1])
    ys = np.real(sol.Y[0, sol.geom.NI + 15 : ex + 1])
    np.testing.assert_allclose(ys, sol.geom.yt + math.tan(math.radians(15)) * xs, rtol=0, atol=1e-5)


# --------------------------------------------------------------------------- scripts
@pytest.mark.parametrize("M", [1.5, 2.0, 2.4, 3.0])
def test_minimum_length_nozzle_exit_area_ratio(M):
    """Fixed MLN contour: exit area ratio within 0.5 % of isentropic A/A* (planar)."""
    res = minimum_length_nozzle(Mach_out=M, num_rays=30)
    assert res.area_ratio == pytest.approx(area_ratio(M, 1.4), rel=0.005)
    assert res.Wall[-1, 4] == pytest.approx(M, rel=1e-9)  # last wall point at the design Mach
    assert res.theta_max == pytest.approx(0.5 * get_prandtl_meyer_function(M, 1.4)[0], rel=1e-14)


def test_minimum_length_nozzle_default_and_original_bug():
    fixed = minimum_length_nozzle()
    original = minimum_length_nozzle(use_tangent_slopes=False)
    assert fixed.area_ratio == pytest.approx(2.40449, abs=1e-5)  # after fix (theory 2.40310)
    assert original.area_ratio == pytest.approx(2.76005, abs=1e-5)  # MATLAB (radian slopes)


def test_expansion_matches_ambient_pressure():
    res = expansion()
    M2 = get_Mach_from_nu(res.nu1 + res.delta, 1.2)
    p2 = 150.0e5 / (1 + 0.1 * M2**2) ** 6
    assert p2 == pytest.approx(101325, rel=1e-10)


def test_oblique_shock_theta_beta_mach():
    res, _ = quiet(oblique_shock)
    # maximum deflection: independent bounded maximisation of theta(beta) on every grid Mach number
    from scipy.optimize import minimize_scalar

    def theta_of_beta(M, b):
        return math.atan(2 / math.tan(b) * ((M * math.sin(b)) ** 2 - 1) / ((1.4 + math.cos(2 * b)) * M**2 + 2))

    for M, beta_max, tmax in res.machmax:
        opt = minimize_scalar(lambda b: -theta_of_beta(M, b), bounds=(math.asin(1 / M), math.pi / 2), method="bounded",
                              options={"xatol": 1e-10})
        assert tmax == pytest.approx(math.degrees(-opt.fun), abs=1e-7)
        assert beta_max == pytest.approx(math.degrees(opt.x), abs=0.01)
    assert res.machmax[-1, 0] == 5.0 and res.machmax[-1, 2] == pytest.approx(41.12, abs=0.01)  # NACA 1135
    # every (theta, beta_weak, beta_strong) satisfies the theta-beta-M relation
    for j, M in enumerate(res.Mach):
        th, b1, b2 = np.radians(res.betadev[j]).T
        np.testing.assert_allclose(equation_Theta_Beta_Mach(1.4, M, th, b1), 0, atol=1e-12)
        np.testing.assert_allclose(equation_Theta_Beta_Mach(1.4, M, th, b2), 0, atol=1e-12)
        assert np.all(b1 < b2)
    # classical value: M = 2, theta = 10 deg -> beta weak = 39.31 deg
    from MODELING.moc._matlab import fzero

    beta = fzero(lambda b: equation_Theta_Beta_Mach(1.4, 2.0, math.radians(10), b), (math.radians(31), math.radians(60)))
    assert math.degrees(beta) == pytest.approx(39.31, abs=0.01)


# --------------------------------------------------------------------------- RDE toolbox
def test_injection_velocity_satisfies_conservation():
    chamber, inj, _ = rde_default_inputs()
    out, _ = quiet(MoC_Injection_Velocity, inj, chamber)
    A4 = math.pi / 4 * (chamber.outer_chamber_diameter**2 - chamber.inner_chamber_diameter**2) * 1e-6
    A3 = inj.injector_area * inj.N_holes * 1e-6
    V3 = inj.mdot * inj.R * inj.T3 / (inj.P3 * A3)
    V4, P4, T4 = out.V4[0], out.P4[0], out.T4[0]
    assert P4 / (inj.R * T4) * V4 * A4 == pytest.approx(inj.mdot, rel=1e-12)  # continuity
    assert (P4 + inj.mdot / A4 * V4) == pytest.approx(  # momentum
        (inj.P3 - inj.Pa) * A3 / A4 + inj.Pa + inj.mdot / A4 * V3, rel=1e-12
    )
    assert T4 + V4**2 / (2 * inj.Cp) == pytest.approx(inj.T3 + V3**2 / (2 * inj.Cp), rel=1e-12)  # energy
    assert out.M4[0] == pytest.approx(V4 / math.sqrt(inj.gamma1 * inj.R * T4), rel=1e-15)
    assert np.isnan(out.M4[1])  # T4 < 0 root


def test_injection_velocity_no_real_root():
    """Negative discriminant (injector area > annulus area, Pa > P3) raises like the MATLAB error()."""
    from MODELING.moc.rde.structs import ChamberDim, InjV

    A4_mm2 = 36.117592829416715
    chamber = ChamberDim(math.sqrt(4 * A4_mm2 / math.pi + 100.0), 10.0, 20.0)
    inj = InjV(injector_area=370.8322120229682, N_holes=1, P3=45325.55305738061, T3=1011.3885208533596,
               mdot=0.011600297194369542, gamma1=1.3, Pa=239227.1118916968, R=3382.7963075177768,
               Cp=15509.100427631705)
    with pytest.raises(ValueError, match="No real velocity"):
        quiet(MoC_Injection_Velocity, inj, chamber)


def test_triple_point_has_no_root_for_default_inputs():
    """Pressure-matched curve never matches the flow angle (residual >= 0.1485 rad)."""
    from scipy.optimize import brentq

    _, _, t = rde_default_inputs()
    g = t.gamma2
    T1p = t.T2 * (t.P1 / t.P2) ** ((g - 1) / g)
    M1p = t.Vcj / math.sqrt(g * t.R1 * T1p)
    resid = []
    for th in np.linspace(math.asin(1 / M1p) + 1e-6, math.pi / 2 - 1e-6, 400):
        M3 = brentq(lambda M: triple_point_residuals([M, th], t.P1, t.P2, g, g, M1p)[0], 1 + 1e-9, 50)
        resid.append(triple_point_residuals([M3, th], t.P1, t.P2, g, g, M1p)[1])
    assert min(resid) == pytest.approx(0.1485, abs=1e-3)
    assert np.all(np.asarray(resid) > 0)


def test_domain_size():
    chamber, _, _ = rde_default_inputs()
    xMax, yMax = MoC_Calculate_Domain_Size(chamber)
    assert xMax == pytest.approx(math.pi * (52.243 + 37.167) / 2, rel=1e-15)
    assert yMax == 29.21
