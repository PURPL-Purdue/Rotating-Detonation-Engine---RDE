"""Helpers shared by the MODELING tests (Octave reference loading)."""

from __future__ import annotations

import contextlib
import io
from pathlib import Path

import numpy as np
import scipy.io as sio

from MODELING.moc.structs import Geom, Params

DATA = Path(__file__).parent / "data"

NOZZLE_CASES = [
    "moc_team",
    "moc_team_te5",
    "moc_team_strongshock",
    "moc_upstream",
    "moc_upstream_te0",
    "moc_upstream_weakshock",
    "moc_upstream_NI21",
]


def load(name: str, structs: bool = False) -> dict:
    """Load an Octave reference file (``structs=True`` keeps MATLAB structs as objects)."""
    if structs:
        return sio.loadmat(DATA / f"{name}.mat", squeeze_me=True, struct_as_record=False)
    return sio.loadmat(DATA / f"{name}.mat")


def inputs_from_reference(name: str):
    """``(geom, params)`` exactly as used by the Octave run of case ``name``."""
    r = load(name, structs=True)
    g, p = r["case_geom"], r["case_params"]
    geom = Geom(
        delta=int(g.delta), yt=float(g.yt), rhou=float(g.rhou), rhod=float(g.rhod), xe=float(g.xe),
        xplume=float(g.xplume), ta=float(g.ta), te=float(g.te), NI=int(g.NI),
        circdownTheta=np.atleast_1d(g.circdownTheta).astype(float), NIexpansion=int(g.NIexpansion),
    )
    params = Params(gamma=float(p.gamma), R=float(p.R), P=float(p.P), PRatio=float(p.PRatio), T=float(p.T))
    return geom, params


def quiet(fn, *args, **kwargs):
    """Call ``fn`` with stdout captured; returns ``(result, printed_text)``."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        out = fn(*args, **kwargs)
    return out, buf.getvalue()


def rel_err(a, b) -> float:
    """max|a-b| / max|b| (arrays of identical shape)."""
    a = np.asarray(a)
    b = np.asarray(b)
    assert a.shape == b.shape, (a.shape, b.shape)
    scale = np.max(np.abs(b))
    return float(np.max(np.abs(a - b)) / (scale if scale > 0 else 1.0))
