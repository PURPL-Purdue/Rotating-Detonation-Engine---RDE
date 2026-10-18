"""Thermochemistry via NASA CEA (cross-platform, no external executable).

This replaces the Windows-only ``FCEA2.exe`` + input-deck + ``.out``-parsing
workflow of the MATLAB wrappers with NASA's official open-source CEA
re-implementation, published on PyPI as ``cea`` (https://github.com/nasa/cea,
wheels for macOS, Linux and Windows). It ships the same ``thermo.lib`` /
``trans.lib`` databases as CEA2 and reproduces the committed CEA2 outputs in
the former ``HADES/SIZING MODEL/CEA`` folder (kept in ``tests/hades_sizing/reference/``)
to every printed digit (see ``tests/hades_sizing/test_rdelib_cea.py``).

Three problem types are wrapped:

* :func:`cj_detonation`  - Chapman-Jouguet detonation (CEA ``det``)
* :func:`rocket`         - rocket performance, infinite-area combustor (``rocket``)
* :func:`equilibrium`    - fixed-state equilibrium (``hp`` or ``tp``)

Units follow CEA: pressure in **bar**, temperature in **K**, enthalpies and
energies in kJ/kg, entropy and cp in kJ/(kg K), velocities in m/s, viscosity
in millipoise, thermal conductivity in mW/(cm K).

Mixtures are given as one fuel and one oxidizer (or lists of them with
relative weights) plus exactly one of ``of`` (O/F weight ratio), ``phi``
(weight-based equivalence ratio, CEA ``phi,eq.ratio``), ``r`` (valence-based
equivalence ratio, CEA ``r,eq.ratio``) or ``pct_fuel`` (weight % fuel).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
from scipy.optimize import nnls

import cea as _cea

from .units import CEA2_PSIA_TO_BAR

__all__ = [
    "DetonationResult",
    "RocketResult",
    "EquilibriumResult",
    "cj_detonation",
    "rocket",
    "equilibrium",
    "reactant_weights",
    "species_elements",
    "cea2_round",
    "CEA2_PSIA_TO_BAR",
]


# ---------------------------------------------------------------------------
# Mixture handling
# ---------------------------------------------------------------------------
def _as_list(x) -> list[str]:
    return [x] if isinstance(x, str) else list(x)


def reactant_weights(
    fuel: str | Sequence[str],
    oxidizer: str | Sequence[str],
    *,
    of: float | None = None,
    phi: float | None = None,
    r: float | None = None,
    pct_fuel: float | None = None,
    fuel_wt: Sequence[float] | None = None,
    ox_wt: Sequence[float] | None = None,
):
    """Build the CEA reactant mixture and the reactant weights.

    Returns ``(names, mixture, weights, of_ratio)`` where ``weights`` are mass
    amounts aligned with ``names`` (fuels first, then oxidizers).
    ``fuel_wt``/``ox_wt`` are the relative weights inside the fuel and
    oxidizer groups (default: equal, like ``wt=100`` on every CEA line).
    """
    fuels = _as_list(fuel)
    oxids = _as_list(oxidizer)
    names = fuels + oxids
    nf, no = len(fuels), len(oxids)
    fw = np.zeros(nf + no)
    ow = np.zeros(nf + no)
    fw[:nf] = 100.0 if fuel_wt is None else np.asarray(fuel_wt, dtype=float)
    ow[nf:] = 100.0 if ox_wt is None else np.asarray(ox_wt, dtype=float)

    given = [v is not None for v in (of, phi, r, pct_fuel)]
    if sum(given) != 1:
        raise ValueError("Specify exactly one mixture type: of, phi, r or pct_fuel.")

    mix = _cea.Mixture(names)
    if of is not None:
        of_ratio = float(of)
    elif phi is not None:
        of_ratio = float(mix.weight_eq_ratio_to_of_ratio(ow, fw, float(phi)))
    elif r is not None:
        of_ratio = float(mix.chem_eq_ratio_to_of_ratio(ow, fw, float(r)))
    else:
        of_ratio = (100.0 - float(pct_fuel)) / float(pct_fuel)
    weights = np.asarray(mix.of_ratio_to_weights(ow, fw, of_ratio), dtype=float)
    return names, mix, weights, of_ratio


def _products(names: list[str], only, omit) -> "_cea.Mixture":
    if only:
        return _cea.Mixture(list(only))
    return _cea.Mixture(names, products_from_reactants=True, omit=list(omit or []))


# --- elemental bookkeeping (pure Python; guards the native solver) ---------
_ELEMENTS = frozenset(
    "H D T He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As "
    "Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd "
    "Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu".split()
)
_PHASE_RE = re.compile(r"\((?:L|cr|gr|a|b|c|d|s|l|g|I|II|III|IV|V|VI)\)$")


def _parse_formula(f: str, i: int, name: str, depth: int = 0):
    counts: dict[str, float] = {}
    while i < len(f):
        c = f[i]
        if c == "(":
            sub, i = _parse_formula(f, i + 1, name, depth + 1)
            mult, i = _parse_number(f, i)
            for el, k in sub.items():
                counts[el] = counts.get(el, 0.0) + k * mult
        elif c == ")":
            if depth == 0:
                raise ValueError(f"cannot parse species formula {name!r}")
            return counts, i + 1
        elif c.isupper():
            j = i + 1
            if j < len(f) and f[j].islower():
                j += 1
            el = f[i:j]
            if el not in _ELEMENTS:
                raise ValueError(f"cannot parse species formula {name!r} (element {el!r})")
            k, i = _parse_number(f, j)
            counts[el] = counts.get(el, 0.0) + k
        else:
            raise ValueError(f"cannot parse species formula {name!r}")
    if depth:
        raise ValueError(f"unbalanced parentheses in species formula {name!r}")
    return counts, i


def _parse_number(f: str, i: int):
    j = i
    while j < len(f) and (f[j].isdigit() or f[j] == "."):
        j += 1
    return (float(f[i:j]) if j > i else 1.0), j


def species_elements(name: str) -> dict[str, float]:
    """Elemental composition of a CEA *product* species, parsed from its name.

    CEA product names are formulas with optional ``,comment`` and phase tags,
    e.g. ``'H2O'``, ``'H2O(L)'``, ``'C(gr)'``, ``'CH3C(CH3)2CH3'``,
    ``'C2H2,acetylene'``. Raises ``ValueError`` for names that are not
    formulas (e.g. reactant blends such as ``'Air'``).
    """
    base = _PHASE_RE.sub("", name.split(",", 1)[0].strip())
    base = base.rstrip("+-")  # ionic charge
    if not base:
        raise ValueError(f"cannot parse species formula {name!r}")
    counts, _ = _parse_formula(base, 0, name)
    return counts


def _reactant_element_moles(names, reac, w) -> dict[str, float]:
    """Element kmol per kg of the reactant mixture.

    Obtained from a standard full-product TP equilibrium (always supported by
    CEA) because elements are conserved; avoids needing formulas for
    reactant blends such as 'Air'.
    """
    prod = _cea.Mixture(names, products_from_reactants=True)
    solver = _cea.EqSolver(prod, reactants=reac)
    sol = _cea.EqSolution(solver)
    solver.solve(sol, _cea.TP, 3000.0, 1.0, w)
    if not sol.converged:
        raise RuntimeError("CEA: could not determine the reactant elemental composition")
    nj = np.asarray(sol.nj, dtype=float)
    cut = 1e-12 * float(np.sum(nj))
    b0: dict[str, float] = {}
    for sp, n in zip(prod.species_names, nj):
        if n <= cut:
            continue
        for el, k in species_elements(sp).items():
            b0[el] = b0.get(el, 0.0) + k * n
    return b0


def _check_product_set(names, reac, w, product_names, rtol: float = 1e-6) -> np.ndarray:
    """Validate a restricted product list *before* any native call.

    Returns non-negative product amounts [kmol/kg] that reproduce the reactant
    elemental composition. Raises ``ValueError`` if a product contains an
    element absent from the reactants, if a reactant element appears in no
    product, or if the product set cannot reproduce the elemental proportions
    (e.g. ``only=['H2O']`` with a non-stoichiometric H2/O2 mixture).
    """
    b0 = _reactant_element_moles(names, reac, w)
    bmax = max(b0.values())
    elems = sorted(el for el, v in b0.items() if v > 1e-12 * bmax)
    comps = [species_elements(p) for p in product_names]
    extra = [p for p, c in zip(product_names, comps) if any(el not in elems for el in c)]
    if extra:
        raise ValueError(f"product species {extra} contain elements not present in the reactants ({elems})")
    missing = [el for el in elems if not any(el in c for c in comps)]
    if missing:
        raise ValueError(f"no product species contains reactant element(s) {missing}")
    A = np.array([[c.get(el, 0.0) for c in comps] for el in elems])
    b = np.array([b0[el] for el in elems])
    n, resid = nnls(A, b)
    if resid > rtol * np.linalg.norm(b):
        raise ValueError(
            "the restricted product set cannot reproduce the reactant elemental proportions "
            f"(relative residual {resid / np.linalg.norm(b):.2e}); e.g. only=['H2O'] requires an "
            "exactly stoichiometric H2/O2 mixture"
        )
    return n


def _guard_products(names, reac, w, prod, only, omit, problem: str) -> np.ndarray | None:
    """Run the product-set checks; returns feasible product amounts or None."""
    if not only and not omit:
        return None
    n = _check_product_set(names, reac, w, prod.species_names)
    if prod.num_species < len(names) and problem != "equilibrium":
        raise NotImplementedError(
            f"CEA {problem}: fewer product species ({prod.num_species}) than reactants ({len(names)}) "
            "is not supported by the CEA 3 solver"
        )
    return n


def _mole_fractions(d) -> dict[str, np.ndarray | float]:
    return {k: (np.array(v, dtype=float) if np.ndim(v) else float(v)) for k, v in d.items()}


# ---------------------------------------------------------------------------
# CJ detonation
# ---------------------------------------------------------------------------
@dataclass
class DetonationResult:
    """Chapman-Jouguet detonation result (CEA ``det`` problem)."""

    # unburned (state 1)
    P1_bar: float
    T1_K: float
    H1_kJkg: float
    M1: float  # molecular weight (1/n) of unburned gas
    gamma1: float
    son_vel1: float  # m/s
    # burned (CJ state)
    P_bar: float
    T_K: float
    rho: float  # kg/m^3
    H_kJkg: float
    U_kJkg: float
    G_kJkg: float
    S_kJkgK: float
    M: float  # molecular weight (1/n)
    MW: float
    gamma_s: float
    son_vel: float  # m/s
    cp_eq: float  # kJ/(kg K)
    cp_frozen: float
    # detonation parameters
    P_P1: float
    T_T1: float
    M_M1: float
    rho_rho1: float
    det_mach: float
    det_vel: float  # CJ velocity, m/s
    # mixture
    of_ratio: float
    # transport (NaN if transport=False)
    visc_mP: float = float("nan")  # millipoise
    k_eq: float = float("nan")  # mW/(cm K)
    k_frozen: float = float("nan")
    Pr_eq: float = float("nan")
    Pr_frozen: float = float("nan")
    mole_fractions: dict = field(default_factory=dict)
    converged: bool = True


def cj_detonation(
    fuel: str | Sequence[str],
    oxidizer: str | Sequence[str],
    P1_bar: float,
    T1_K: float,
    *,
    of: float | None = None,
    phi: float | None = None,
    r: float | None = None,
    pct_fuel: float | None = None,
    fuel_wt: Sequence[float] | None = None,
    ox_wt: Sequence[float] | None = None,
    transport: bool = True,
    only: Sequence[str] | None = None,
    omit: Sequence[str] | None = None,
    frozen: bool = False,
) -> DetonationResult:
    """Chapman-Jouguet detonation of a fuel/oxidizer mixture at (P1, T1).

    Equivalent CEA2 deck::

        problem det p,bar=P1_bar t,k=T1_K phi=...
        reac oxid <oxidizer> wt=100  fuel <fuel> wt=100
        output transport

    Parameters
    ----------
    P1_bar, T1_K : initial (unburned) pressure [bar] and temperature [K]
    of / phi / r / pct_fuel : mixture specification (exactly one)
    transport : also compute viscosity, conductivity and Prandtl numbers
    only, omit : restrict / exclude product species (CEA ``only`` / ``omit``)
    """
    names, reac, w, of_ratio = reactant_weights(
        fuel, oxidizer, of=of, phi=phi, r=r, pct_fuel=pct_fuel, fuel_wt=fuel_wt, ox_wt=ox_wt
    )
    prod = _products(names, only, omit)
    _guard_products(names, reac, w, prod, only, omit, "detonation")
    solver = _cea.DetonationSolver(prod, reactants=reac, transport=transport)
    s = _cea.DetonationSolution(solver)
    solver.solve(s, w, float(T1_K), float(P1_bar), frozen=frozen)
    if not s.converged:
        raise RuntimeError("CEA detonation solution did not converge.")
    nan = float("nan")
    return DetonationResult(
        P1_bar=float(P1_bar),
        T1_K=float(T1_K),
        H1_kJkg=float(s.H1),
        M1=float(s.M1),
        gamma1=float(s.gamma1),
        son_vel1=float(s.sonic_velocity1),
        P_bar=float(s.P),
        T_K=float(s.T),
        rho=float(s.density),
        H_kJkg=float(s.enthalpy),
        U_kJkg=float(s.energy),
        G_kJkg=float(s.gibbs_energy),
        S_kJkgK=float(s.entropy),
        M=float(s.M),
        MW=float(s.MW),
        gamma_s=float(s.gamma_s),
        son_vel=float(s.sonic_velocity),
        cp_eq=float(s.cp_eq),
        cp_frozen=float(s.cp_fr),
        P_P1=float(s.P_P1),
        T_T1=float(s.T_T1),
        M_M1=float(s.M_M1),
        rho_rho1=float(s.rho_rho1),
        det_mach=float(s.Mach),
        det_vel=float(s.velocity),
        of_ratio=of_ratio,
        visc_mP=float(s.viscosity) if transport else nan,
        k_eq=float(s.conductivity_eq) if transport else nan,
        k_frozen=float(s.conductivity_fr) if transport else nan,
        Pr_eq=float(s.Pr_eq) if transport else nan,
        Pr_frozen=float(s.Pr_fr) if transport else nan,
        mole_fractions=_mole_fractions(s.mole_fractions),
        converged=bool(s.converged),
    )


# ---------------------------------------------------------------------------
# Rocket
# ---------------------------------------------------------------------------
@dataclass
class RocketResult:
    """Rocket performance result. Arrays are indexed by station.

    Station order is CEA's: ``[chamber, throat, *pi_p points, *subar points,
    *supar points]`` (``labels`` names them). Performance arrays (``ae_at``,
    ``cstar``, ``cf``, ``isp``, ``ivac``) are 0 at the chamber, as in CEA.
    """

    labels: list[str]
    P_bar: np.ndarray
    T_K: np.ndarray
    rho: np.ndarray
    H_kJkg: np.ndarray
    U_kJkg: np.ndarray
    G_kJkg: np.ndarray
    S_kJkgK: np.ndarray
    M: np.ndarray
    MW: np.ndarray
    cp_eq: np.ndarray
    cp_frozen: np.ndarray
    gamma_s: np.ndarray
    son_vel: np.ndarray
    mach: np.ndarray
    ae_at: np.ndarray
    cstar: np.ndarray  # m/s
    cf: np.ndarray
    ivac: np.ndarray  # m/s
    isp: np.ndarray  # m/s (ambient pressure = exit pressure)
    of_ratio: float
    Pc_bar: float
    mole_fractions: dict = field(default_factory=dict)
    converged: bool = True

    @property
    def pinf_p(self) -> np.ndarray:
        """Pinf/P for every station."""
        return self.Pc_bar / self.P_bar


def rocket(
    fuel: str | Sequence[str],
    oxidizer: str | Sequence[str],
    Pc_bar: float,
    *,
    supar: float | Sequence[float] | None = None,
    subar: float | Sequence[float] | None = None,
    pi_p: float | Sequence[float] | None = None,
    of: float | None = None,
    phi: float | None = None,
    r: float | None = None,
    pct_fuel: float | None = None,
    fuel_wt: Sequence[float] | None = None,
    ox_wt: Sequence[float] | None = None,
    hc: float | None = None,
    T_reac: float | Sequence[float] | None = None,
    n_frz: int | None = None,
    only: Sequence[str] | None = None,
    omit: Sequence[str] | None = None,
    transport: bool = False,
) -> RocketResult:
    """Equilibrium rocket performance (infinite-area combustor).

    Equivalent CEA2 deck::

        problem rocket equilibrium p,bar=Pc_bar supar=... (or pi/p=...) phi=...
        reac oxid <oxidizer> wt=100  fuel <fuel> wt=100

    Chamber enthalpy: pass ``hc`` (h/R, kg-mol K/kg, CEA convention) or the
    reactant temperature(s) ``T_reac`` [K]. If neither is given the reactant
    enthalpy is taken as **zero**, which is what CEA2 did for the HADES decks
    (reactant lines with no ``t=``; the ``.out`` file shows ENERGY = 0, TEMP = 0).

    ``n_frz`` freezes composition at the given station (1 = chamber) for a
    frozen-flow solution; ``None`` gives the equilibrium solution.
    """
    names, reac, w, of_ratio = reactant_weights(
        fuel, oxidizer, of=of, phi=phi, r=r, pct_fuel=pct_fuel, fuel_wt=fuel_wt, ox_wt=ox_wt
    )
    if hc is None:
        if T_reac is not None:
            hc = float(reac.calc_property(_cea.ENTHALPY, w, T_reac) / _cea.R)
        else:
            hc = 0.0
    prod = _products(names, only, omit)
    _guard_products(names, reac, w, prod, only, omit, "rocket")
    solver = _cea.RocketSolver(prod, reactants=reac, transport=transport)
    s = _cea.RocketSolution(solver)

    def _arr(v):
        if v is None:
            return None
        return [float(x) for x in np.atleast_1d(v)]

    solver.solve(
        s, w, float(Pc_bar), pi_p=_arr(pi_p), subar=_arr(subar), supar=_arr(supar),
        iac=True, n_frz=n_frz, hc=float(hc),
    )
    if not s.converged:
        raise RuntimeError("CEA rocket solution did not converge.")
    labels = ["chamber", "throat"]
    labels += [f"pi_p={x:g}" for x in (_arr(pi_p) or [])]
    labels += [f"subar={x:g}" for x in (_arr(subar) or [])]
    labels += [f"supar={x:g}" for x in (_arr(supar) or [])]
    a = lambda v: np.array(v, dtype=float)  # noqa: E731
    return RocketResult(
        labels=labels,
        P_bar=a(s.P),
        T_K=a(s.T),
        rho=a(s.density),
        H_kJkg=a(s.enthalpy),
        U_kJkg=a(s.energy),
        G_kJkg=a(s.gibbs_energy),
        S_kJkgK=a(s.entropy),
        M=a(s.M),
        MW=a(s.MW),
        cp_eq=a(s.cp_eq),
        cp_frozen=a(s.cp_fr),
        gamma_s=a(s.gamma_s),
        son_vel=a(s.sonic_velocity),
        mach=a(s.Mach),
        ae_at=a(s.ae_at),
        cstar=a(s.c_star),
        cf=a(s.coefficient_of_thrust),
        ivac=a(s.Isp_vacuum),
        isp=a(s.Isp),
        of_ratio=of_ratio,
        Pc_bar=float(Pc_bar),
        mole_fractions=_mole_fractions(s.mole_fractions),
        converged=bool(s.converged),
    )


# ---------------------------------------------------------------------------
# Fixed-state equilibrium (hp / tp)
# ---------------------------------------------------------------------------
@dataclass
class EquilibriumResult:
    """Equilibrium state (CEA ``hp`` or ``tp`` problem)."""

    P_bar: float
    T_K: float
    rho: float
    H_kJkg: float
    U_kJkg: float
    G_kJkg: float
    S_kJkgK: float
    M: float
    MW: float
    cp_eq: float
    cp_frozen: float
    gamma_s: float
    son_vel: float
    of_ratio: float
    mole_fractions: dict = field(default_factory=dict)
    converged: bool = True


def equilibrium(
    fuel: str | Sequence[str],
    oxidizer: str | Sequence[str],
    P_bar: float,
    T_K: float,
    *,
    problem: str = "hp",
    of: float | None = None,
    phi: float | None = None,
    r: float | None = None,
    pct_fuel: float | None = None,
    fuel_wt: Sequence[float] | None = None,
    ox_wt: Sequence[float] | None = None,
    only: Sequence[str] | None = None,
    omit: Sequence[str] | None = None,
) -> EquilibriumResult:
    """Equilibrium composition and properties.

    ``problem='hp'``: adiabatic, constant-pressure combustion of reactants at
    temperature ``T_K``. ``problem='tp'``: equilibrium at fixed ``T_K``, ``P_bar``.
    """
    names, reac, w, of_ratio = reactant_weights(
        fuel, oxidizer, of=of, phi=phi, r=r, pct_fuel=pct_fuel, fuel_wt=fuel_wt, ox_wt=ox_wt
    )
    if problem.lower() not in ("hp", "tp"):
        raise ValueError("problem must be 'hp' or 'tp'")
    prod = _products(names, only, omit)
    n_feasible = _guard_products(names, reac, w, prod, only, omit, "equilibrium")
    if n_feasible is not None and prod.num_species < len(names):
        # CEA 3 cannot attach a reactant list longer than the product list
        # (CEA2 reduced the components instead). Equivalent problem: give the
        # amounts as product weights with the same elemental composition and
        # use the reactant enthalpy as the hp target.
        solver = _cea.EqSolver(prod)
        amounts = np.asarray(prod.moles_to_weights(n_feasible), dtype=float)
    else:
        solver = _cea.EqSolver(prod, reactants=reac)
        amounts = w
    s = _cea.EqSolution(solver)
    if problem.lower() == "hp":
        h0 = reac.calc_property(_cea.ENTHALPY, w, float(T_K)) / _cea.R
        solver.solve(s, _cea.HP, float(h0), float(P_bar), amounts)
    else:
        solver.solve(s, _cea.TP, float(T_K), float(P_bar), amounts)
    if not s.converged:
        raise RuntimeError("CEA equilibrium solution did not converge.")
    T = float(s.T)
    MW = float(s.MW)
    gam = float(s.gamma_s)
    son = float(np.sqrt(gam * 8314.46261815324 / float(s.M) * T))
    return EquilibriumResult(
        P_bar=float(P_bar), T_K=T, rho=float(s.density), H_kJkg=float(s.enthalpy),
        U_kJkg=float(s.energy), G_kJkg=float(s.gibbs_energy), S_kJkgK=float(s.entropy),
        M=float(s.M), MW=MW, cp_eq=float(s.cp_eq), cp_frozen=float(s.cp_fr), gamma_s=gam,
        son_vel=son, of_ratio=of_ratio, mole_fractions=_mole_fractions(s.mole_fractions),
        converged=bool(s.converged),
    )


# ---------------------------------------------------------------------------
# CEA2 printed-precision rounding
# ---------------------------------------------------------------------------
def _varfmt_decimals(x: float) -> int:
    """Decimals CEA2's VARFMT routine uses for a value (F9.x formats)."""
    v = abs(x)
    if v >= 1.0e6:
        return 0
    if v >= 1.0e4:
        return 1
    if v >= 100.0:
        return 2
    if v >= 10.0:
        return 3
    if v >= 1.0:
        return 4
    return 5


def cea2_round(x: float, decimals: int | str) -> float:
    """Round ``x`` the way CEA2 printed it in its ``.out`` file.

    The MATLAB wrappers parsed the printed ``.out`` file, so every value they
    returned carried the printed precision. ``decimals`` is either a fixed
    number of decimals (e.g. 1 for velocities) or ``'var'`` for quantities
    printed with CEA2's magnitude-dependent VARFMT format (P, H, U, G,
    viscosity, conductivity, Ae/At).
    """
    if x is None or not np.isfinite(x):
        return x
    d = _varfmt_decimals(x) if decimals == "var" else int(decimals)
    return float(f"{x:.{d}f}")
