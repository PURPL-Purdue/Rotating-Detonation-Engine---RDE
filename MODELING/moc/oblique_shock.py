"""Oblique-shock theta-beta-Mach curves (weak and strong branches).

Python translation of the script ``oblique_shock.m`` (src2/) from
xavierdechamps/Nozzle_MOC_supersonic (GPL-3.0); see MODELING/moc/README.md.

``fminbnd`` -> ``scipy.optimize.fminbound`` (same Brent golden-section /
parabolic algorithm) with MATLAB's default ``TolX = 1e-4``; ``fzero`` uses the
MATLAB algorithm port.  Porting fix (see PORTING_NOTES.md): the final axis
labels were swapped in the MATLAB (x is the deflection angle theta).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import fminbound

from ._matlab import colon, fzero, num2str


def equation_max_theta(gamma, Mach, beta):
    """cot(theta) from the theta-beta-M relation (minimum -> maximum deflection)."""
    return np.tan(beta) * (0.5 * ((gamma + np.cos(2 * beta)) * Mach**2 + 2)) / ((Mach * np.sin(beta)) ** 2 - 1)


def equation_Theta_Beta_Mach(gamma, Mach, theta, beta):
    """Residual of the theta-beta-M relation: tan(theta(beta)) - tan(theta)."""
    return 2 * (1.0 / np.tan(beta)) * ((Mach * np.sin(beta)) ** 2 - 1) / ((gamma + np.cos(2 * beta)) * Mach**2 + 2) - np.tan(theta)


def _acot(x):
    return math.atan(1.0 / x)


@dataclass
class ObliqueShockResult:
    """``machmax`` rows: [Mach, beta_max (deg), theta_max (deg)];
    ``betadev[j, i, :]`` = [theta, beta_weak, beta_strong] (deg) for ``Mach[j]``;
    ``n_valid[j]`` = number of filled ``i`` entries (``i-1`` in MATLAB)."""

    Mach: np.ndarray
    machmax: np.ndarray
    betadev: np.ndarray
    n_valid: list
    theta_max: float  # [rad] maximum deflection of the last Mach number


def oblique_shock(Mach=(2, 3, 5), gamma=1.4, theta_min=1 * math.pi / 180, theta_steps=100, plot=False):
    """Compute (and optionally plot, like the MATLAB script) the theta-beta-M curves."""
    Mach = np.asarray(Mach, dtype=float)
    Mach_min = 1.1
    Mach_max = Mach[-1]
    Mach_theta_max = colon(Mach_min, (Mach_max - Mach_min) / 100, Mach_max)
    machmax = np.zeros((Mach_theta_max.size, 3))
    for k, mach_j in enumerate(Mach_theta_max):
        beta_max = fminbound(lambda beta: equation_max_theta(gamma, mach_j, beta), math.pi / 3, math.pi / 2 - 0.0001, xtol=1e-4)
        theta_max = _acot(equation_max_theta(gamma, mach_j, beta_max))
        machmax[k, :] = [mach_j, beta_max * 180 / math.pi, theta_max * 180 / math.pi]

    if plot:
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(111)
    hold = False

    def _plot(*args, **kw):
        if not hold:
            ax.cla()
        ax.plot(*args, **kw)
        ax.grid(True)

    betadev_rows = []
    n_valid = []
    for mach_j in Mach:
        # theta_max and beta (rad) corresponding to theta_max
        beta_max = fminbound(lambda beta: equation_max_theta(gamma, mach_j, beta), math.pi / 4, math.pi / 2 - 0.001, xtol=1e-4)
        theta_max = _acot(equation_max_theta(gamma, mach_j, beta_max))

        theta_max_search = theta_max - 0.00001
        theta = colon(theta_min, (theta_max_search - theta_min) / theta_steps, theta_max_search)
        rows = np.zeros((theta.size, 3))
        i = 1
        for theta_i in theta:
            try:
                guessBeta = (math.pi / 180, beta_max)
                if plot:
                    beta_r = colon(guessBeta[0], math.pi / 100, guessBeta[1])
                    _plot(beta_r, equation_Theta_Beta_Mach(gamma, mach_j, theta_i, beta_r))
                betaEval_1 = fzero(lambda beta: equation_Theta_Beta_Mach(gamma, mach_j, theta_i, beta), guessBeta)
                guessBeta = (betaEval_1 + 0.0001, math.pi / 2)
                betaEval_2 = fzero(lambda beta: equation_Theta_Beta_Mach(gamma, mach_j, theta_i, beta), guessBeta)
                rows[i - 1, :] = np.array([theta_i, betaEval_1, betaEval_2]) * 180 / math.pi
                i += 1
            except ValueError:
                print("fzero failed")
                break
        betadev_rows.append(rows)
        n_valid.append(i - 1)
        if plot:
            _plot(rows[: i - 1, 0], rows[: i - 1, 1], "r")
            ax.plot(rows[: i - 1, 0], rows[: i - 1, 2], "b")
            hold = True

    width = max(r.shape[0] for r in betadev_rows)
    betadev = np.zeros((Mach.size, width, 3))
    for j, r in enumerate(betadev_rows):
        betadev[j, : r.shape[0], :] = r

    if plot:
        ax.plot(machmax[:, 2], machmax[:, 1], "k")
        ax.set_title("Mach number = " + num2str(Mach) + ", max theta = " + num2str(theta_max * 180 / math.pi))
        ax.set_xlabel("Deflection angle theta (deg)")
        ax.set_ylabel("Shock wave angle beta (deg)")
    return ObliqueShockResult(Mach, machmax, betadev, n_valid, float(theta_max))


def main():
    """Equivalent of running ``oblique_shock.m``."""
    import matplotlib.pyplot as plt

    res = oblique_shock(plot=True)
    plt.show()
    return res


if __name__ == "__main__":
    main()
