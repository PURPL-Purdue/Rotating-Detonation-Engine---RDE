"""Minimum-length (sharp-throat) planar nozzle by the method of characteristics.

Python translation of the script ``minimum_length_nozzle.m`` (src2/) from
xavierdechamps/Nozzle_MOC_supersonic (GPL-3.0); see MODELING/moc/README.md.

Porting fix (see PORTING_NOTES.md): the MATLAB used the mean characteristic /
wall *angle in radians* as the slope of the straight segments instead of
``tan(angle)``, which overestimates the exit area ratio (A/A* = 2.760 instead of
the isentropic 2.403 for M = 2.4, gamma = 1.4).  ``use_tangent_slopes=False``
reproduces the original MATLAB numbers exactly.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from ._matlab import colon, tand
from .get_Mach_from_nu import get_Mach_from_nu
from .get_prandtl_meyer_function import get_prandtl_meyer_function


@dataclass
class MinimumLengthNozzle:
    """Workspace of ``minimum_length_nozzle.m`` (MATLAB variable names; ``abs`` -> ``abs_``)."""

    theta_max: float
    nu_out: float
    theta: np.ndarray
    nu: np.ndarray
    Mach: np.ndarray
    mu: np.ndarray
    Kminus: np.ndarray
    Kplus: np.ndarray
    abs_: np.ndarray  # abscissae of the characteristic intersections (MATLAB ``abs``)
    ord: np.ndarray  # ordinates of the characteristic intersections
    Wall: np.ndarray  # [K- K+ theta nu Mach mu] at each wall point
    absWall: np.ndarray
    ordWall: np.ndarray
    absp: np.ndarray
    ordp: np.ndarray
    Machp: np.ndarray
    area_ratio: float


def minimum_length_nozzle(gamma=1.4, Mach_out=2.4, num_rays=7, throat_rad=1.0, use_tangent_slopes=True):
    """Compute the minimum-length nozzle contour (defaults = values of the MATLAB script)."""
    rad2deg = 180 / math.pi

    def slope(angle_deg):
        # MATLAB original: angle_deg / rad2deg (angle in radians used as a slope)
        return tand(angle_deg) if use_tangent_slopes else angle_deg / rad2deg

    nu_out, mu_out = get_prandtl_meyer_function(Mach_out, gamma)
    theta_max = 0.5 * nu_out

    n = num_rays
    theta = np.zeros((n, n))
    nu = np.zeros((n, n))
    Mach = np.zeros((n, n))
    mu = np.zeros((n, n))
    Kminus = np.zeros((n, n))
    Kplus = np.zeros((n, n))
    abs_ = np.zeros((n, n))
    ord_ = np.zeros((n, n))
    Wall = np.zeros((n, 6))
    absWall = np.zeros(n)
    ordWall = np.zeros(n)

    theta_init = theta_max - math.floor(theta_max)
    theta[0, :] = colon(theta_init, (theta_max - theta_init) / (num_rays - 1), theta_max)
    nu[0, :] = theta[0, :]
    Mach[0, :] = get_Mach_from_nu(nu[0, :], gamma)
    _, mu[0, :] = get_prandtl_meyer_function(Mach[0, :], gamma)

    Kminus[0, :] = theta[0, :] + nu[0, :]
    Kplus[0, :] = theta[0, :] - nu[0, :]

    # 1-based (j, k) of the MATLAB script -> [j-1, k-1]
    abs_[0, 0] = throat_rad / tand(abs(theta[0, 0] - mu[0, 0]))
    ord_[0, 0] = 0.0
    for k in range(2, num_rays + 1):
        meanSlopeplus = slope(0.5 * (theta[0, k - 2] + theta[0, k - 1] + mu[0, k - 2] + mu[0, k - 1]))
        meanSlopeminus = slope(theta[0, k - 1] - mu[0, k - 1])
        abs_[0, k - 1] = (throat_rad - ord_[0, k - 2] + abs_[0, k - 2] * meanSlopeplus) / (
            meanSlopeplus - meanSlopeminus
        )
        ord_[0, k - 1] = throat_rad + meanSlopeminus * abs_[0, k - 1]

    k = num_rays
    Wall[0, :] = [Kminus[0, k - 1], Kplus[0, k - 1], theta[0, k - 1], nu[0, k - 1], Mach[0, k - 1], mu[0, k - 1]]
    meanSlopeplus = slope(0.5 * (theta[0, k - 1] + Wall[0, 2] + mu[0, k - 1] + Wall[0, 5]))
    meanSlopeWall = slope(0.5 * (theta_max + Wall[0, 2]))
    absWall[0] = (throat_rad - ord_[0, k - 1] + abs_[0, k - 1] * meanSlopeplus) / (meanSlopeplus - meanSlopeWall)
    ordWall[0] = throat_rad + meanSlopeWall * absWall[0]

    for j in range(2, num_rays + 1):
        lk = num_rays - j + 1
        for k in range(1, lk + 1):
            J, K = j - 1, k - 1
            if k == 1:
                theta[J, K] = 0.0
                Kminus[J, K] = Kminus[J - 1, K + 1]
                nu[J, K] = Kminus[J, K] - theta[J, K]
                Kplus[J, K] = theta[J, K] - nu[J, K]
                Mach[J, K] = get_Mach_from_nu(nu[J, K], gamma)
                _, mu[J, K] = get_prandtl_meyer_function(Mach[J, K], gamma)

                meanSlopeminus = slope(0.5 * (theta[J - 1, K + 1] + theta[J, K] - mu[J - 1, K + 1] - mu[J, K]))
                abs_[J, K] = abs_[J - 1, K + 1] - ord_[J - 1, K + 1] / meanSlopeminus
                ord_[J, K] = 0.0
            else:
                Kplus[J, K] = Kplus[J, K - 1]
                Kminus[J, K] = Kminus[J - 1, K + 1]
                theta[J, K] = 0.5 * (Kminus[J, K] + Kplus[J, K])
                nu[J, K] = 0.5 * (Kminus[J, K] - Kplus[J, K])
                Mach[J, K] = get_Mach_from_nu(nu[J, K], gamma)
                _, mu[J, K] = get_prandtl_meyer_function(Mach[J, K], gamma)

                meanSlopeplus = slope(0.5 * (theta[J, K - 1] + theta[J, K] + mu[J, K - 1] + mu[J, K]))
                meanSlopeminus = slope(0.5 * (theta[J - 1, K + 1] + theta[J, K] - mu[J - 1, K + 1] - mu[J, K]))
                abs_[J, K] = (
                    ord_[J - 1, K + 1] - ord_[J, K - 1] + abs_[J, K - 1] * meanSlopeplus
                    - abs_[J - 1, K + 1] * meanSlopeminus
                ) / (meanSlopeplus - meanSlopeminus)
                ord_[J, K] = ord_[J - 1, K + 1] + meanSlopeminus * (abs_[J, K] - abs_[J - 1, K + 1])

        L = lk - 1
        J = j - 1
        Wall[J, :] = [Kminus[J, L], Kplus[J, L], theta[J, L], nu[J, L], Mach[J, L], mu[J, L]]

        meanSlopeplus = slope(0.5 * (theta[J, L] + Wall[J, 2] + mu[J, L] + Wall[J, 5]))
        meanSlopeWall = slope(0.5 * (Wall[J - 1, 2] + Wall[J, 2]))
        absWall[J] = (
            ordWall[J - 1] - ord_[J, L] + abs_[J, L] * meanSlopeplus - absWall[J - 1] * meanSlopeWall
        ) / (meanSlopeplus - meanSlopeWall)
        ordWall[J] = ordWall[J - 1] + meanSlopeWall * (absWall[J] - absWall[J - 1])

    absp, ordp, Machp = [], [], []
    length = num_rays
    for j in range(num_rays):
        absp.extend(abs_[j, :length])
        ordp.extend(ord_[j, :length])
        Machp.extend(Mach[j, :length])
        length -= 1

    area_ratio = ordWall[-1] / throat_rad
    return MinimumLengthNozzle(
        theta_max=float(theta_max), nu_out=float(nu_out), theta=theta, nu=nu, Mach=Mach, mu=mu,
        Kminus=Kminus, Kplus=Kplus, abs_=abs_, ord=ord_, Wall=Wall, absWall=absWall, ordWall=ordWall,
        absp=np.asarray(absp), ordp=np.asarray(ordp), Machp=np.asarray(Machp), area_ratio=float(area_ratio),
    )


def plot_minimum_length_nozzle(res: MinimumLengthNozzle, throat_rad=1.0):
    """Scatter of the characteristic nodes coloured by Mach number + wall contour."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    sc = ax.scatter(res.absp, res.ordp, 25, res.Machp)
    ax.grid(True)
    fig.colorbar(sc, ax=ax)
    ax.set_ylim(0, res.ordWall.max())
    ax.plot(np.concatenate([[0.0], res.absWall]), np.concatenate([[throat_rad], res.ordWall]), "k")
    return fig


def main():
    """Equivalent of running ``minimum_length_nozzle.m`` (with the slope fix)."""
    import matplotlib.pyplot as plt

    res = minimum_length_nozzle()
    plot_minimum_length_nozzle(res)
    print(f"area_ratio = {res.area_ratio:.15e}")
    plt.show()
    return res


if __name__ == "__main__":
    main()
