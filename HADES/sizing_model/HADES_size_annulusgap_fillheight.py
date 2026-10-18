"""Port of ``HADES_size_annulusgap_fillheight.m`` (cell width, annulus gap, fill height).

Author of the MATLAB: anjali. Cell-size data from the Caltech detonation
database (https://shepherd.caltech.edu/detn_db/html/H2-Air2.html), stored in
``HADES_size_cellsizedata.xlsx`` next to this module.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

__all__ = [
    "HADES_size_annulusgap_fillheight",
    "read_cellsize_data",
    "cellwidth_from_caltech_fit",
    "kawashima_fill_height",
]

CELLSIZE_XLSX = Path(__file__).with_name("HADES_size_cellsizedata.xlsx")


def read_cellsize_data(path: str | Path | None = None) -> np.ndarray:
    """``readmatrix("HADES_size_cellsizedata.xlsx")``: numeric rows only.

    Returns an (n, 2) array: equivalence ratio, detonation cell width [mm].
    """
    import openpyxl

    wb = openpyxl.load_workbook(path or CELLSIZE_XLSX, data_only=True, read_only=True)
    ws = wb.worksheets[0]
    rows = []
    for row in ws.iter_rows(values_only=True):
        vals = [v for v in row[:2]]
        if all(isinstance(v, (int, float)) for v in vals):
            rows.append([float(v) for v in vals])
    wb.close()
    return np.array(rows)


def _log_fit(data: np.ndarray):
    log_phis = np.log10(data[:, 0])
    log_cell_width = np.log10(data[:, 1])
    coeffs = np.polyfit(log_phis, log_cell_width, 3)
    return log_phis, log_cell_width, coeffs


def cellwidth_from_caltech_fit(phi, path: str | Path | None = None):
    """Cell width [mm] from the cubic log-log fit of the Caltech data.

    This is the calculation of the superseded script ``cell_size2.m``
    (``10^polyval(coeffs, log10(phi))``), kept so no logic is lost.
    """
    _, _, coeffs = _log_fit(read_cellsize_data(path))
    return 10 ** np.polyval(coeffs, np.log10(phi))


def kawashima_fill_height(L_theta, N_det, D_CJ, a2, P2, Pc, u0):
    """Detonation (fill) height [m], Kawashima Eq. (15) with Eq. (9).

    L_theta azimuthal annulus length [m], N_det wave number [-], D_CJ CJ speed
    [m/s], a2 burned-gas sound speed [m/s], P2 burned-gas pressure and Pc
    injector critical pressure (same units), u0 unburned axial velocity [m/s].
    """
    x = np.asarray(P2) / Pc
    f_val = (
        (231 / 2**10) * x ** (1 / 15)
        + (63 / 2**9) * x ** (3 / 15)
        + (105 / 2**10) * x ** (5 / 15)
        + (25 / 2**8) * x ** (7 / 15)
        + (105 / 2**10) * x ** (9 / 15)
        + (63 / 2**9) * x ** (11 / 15)
        + (231 / 2**10) * x ** (13 / 15)
    )
    bracket = (1 / a2) * f_val + 1 / u0
    return (L_theta / (np.asarray(N_det) * D_CJ)) * (1 / bracket)


def HADES_size_annulusgap_fillheight(chosen_phi, L_theta, N_det, D_CJ, a2, P2, Pc, u0, *,
                                     plot: bool = True, data_file: str | Path | None = None):
    """Cell width, chosen annulus gap and fill height.

    Parameters
    ----------
    chosen_phi : equivalence ratio [-]
    L_theta : azimuthal annulus length [m]
    N_det : detonation wave number [-]
    D_CJ : Chapman-Jouguet detonation speed [m/s]
    a2 : sound speed in burned gas [m/s]
    P2 : pressure of burned gas region [Pa or bar, consistent with Pc]
    Pc : injector critical pressure [same units as P2]; also used (as bar,
        converted to atm) in the cell-width correlation
    u0 : axial velocity of unburned propellant [m/s]

    Returns
    -------
    (chosen_cellwidth [mm], chosen_annulus_gap [mm], chosen_fillheight [m])

    ``chosen_annulus_gap`` is the HADES design choice, 3/8 in = 9.525 mm, not
    the 1.4 x cell-width correlation value (which is only printed).
    """
    data = read_cellsize_data(data_file)
    log_phis, log_cell_width, coeffs = _log_fit(data)
    log_phis_fit = np.linspace(log_phis.min(), log_phis.max(), 100)
    log_cw_fit = np.polyval(coeffs, log_phis_fit)

    if plot:
        import matplotlib.pyplot as plt

        plt.figure()
        plt.scatter(log_phis, log_cell_width)
        plt.plot(log_phis_fit, log_cw_fit)
        plt.xlabel("log(equivalence ratio)")
        plt.ylabel("log(detonation cell width)")
        plt.title("detonation cell width vs equivalence ratio")
        plt.grid(True)

    # cell width correlation (pressure converted bar -> atm)
    chosen_cellwidth = -0.435 + 8.33 * (Pc * 0.98692327) ** -0.4 + 0.219 * chosen_phi**-7

    # annulus gap correlation from Karasu (continuous spin detonation)
    annulus_gap = 1.4 * chosen_cellwidth
    chosen_annulus_gap = 9.525  # mm, HADES chose 3/8 in

    print("Calculated Cell Width: %.3f mm " % chosen_cellwidth)
    print("Calculated Annulus Gap: %.3f mm " % annulus_gap)

    # fill height (Kawashima Eq. 15)
    chosen_fillheight = float(kawashima_fill_height(L_theta, N_det, D_CJ, a2, P2, Pc, u0))

    if plot:
        import matplotlib.pyplot as plt

        N_plot = np.array([1, 2, 3, 4, 5, 6], dtype=float)
        h_plot = kawashima_fill_height(L_theta, N_plot, D_CJ, a2, P2, Pc, u0)
        # fit(N, h, fittype('a./x')) is linear in a: least-squares a = sum(h/N)/sum(1/N^2)
        a_fit = np.sum(h_plot / N_plot) / np.sum(1 / N_plot**2)
        plt.figure()
        plt.scatter(N_plot, h_plot)
        xs = np.linspace(N_plot.min(), N_plot.max(), 100)
        plt.plot(xs, a_fit / xs)
        plt.plot([N_det], [chosen_fillheight])
        plt.xlabel("Wave #")
        plt.ylabel("Fill Height")
        plt.title("Fill height as function of wave #")
        plt.grid(True)

    print("Fill Height: %.3f m (%.3f mm). " % (chosen_fillheight, chosen_fillheight * 1000))
    return chosen_cellwidth, chosen_annulus_gap, chosen_fillheight
